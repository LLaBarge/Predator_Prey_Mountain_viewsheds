library(inlabru)
library(INLA)
library(dplyr)
library(amt)
library(raster)
library(sf)

# Load data and rasters
dem <- raster("dem.tif")
canopy_height_raster <- crop(raster("meta_tree_height.tif"), dem)
cc_raster <- raster("canopy_cover.tif")
vervet_clean <- readRDS("vervet_clean.rds")
baboon_clean <- readRDS("baboon_clean.rds")

if (!("x_utm" %in% names(vervet_clean))) {
  vervet_sf <- st_as_sf(vervet_clean, coords = c("longitude", "latitude"), crs = 4326)
  vervet_utm <- st_transform(vervet_sf, crs = 32735)
  coords <- st_coordinates(vervet_utm)
  vervet_clean$x_utm <- coords[, 1]
  vervet_clean$y_utm <- coords[, 2]
}

if (!("x_utm" %in% names(baboon_clean))) {
  baboon_sf <- st_as_sf(baboon_clean, coords = c("longitude", "latitude"), crs = 4326)
  baboon_utm <- st_transform(baboon_sf, crs = 32735)
  coords <- st_coordinates(baboon_utm)
  baboon_clean$x_utm <- coords[, 1]
  baboon_clean$y_utm <- coords[, 2]
}

if ("New_Timestamp" %in% names(vervet_clean)) {
  vervet_clean$timestamp <- vervet_clean$New_Timestamp
  baboon_clean$timestamp <- baboon_clean$New_Timestamp
}

gps_dates <- unique(as.Date(vervet_clean$timestamp, tz = "UTC"))

# Load and parse scan data for observer height
scans <- read.csv("scans_com.csv", stringsAsFactors = FALSE, quote = "\"")

scan_date_col <- names(scans)[grepl("Date", names(scans), ignore.case = TRUE)][1]
scan_time_col <- names(scans)[grepl("Time", names(scans), ignore.case = TRUE)][1]
scan_height_col <- names(scans)[grepl("Height", names(scans), ignore.case = TRUE)][1]

clean_height <- function(h) {
  h <- trimws(as.character(h))
  if (toupper(h) == "UN" || h == "" || is.na(h)) return(NA)
  if (grepl("^>\\s*10", h, ignore.case = TRUE)) return(11)
  h_clean <- gsub("[^0-9.]", "", h)
  height <- as.numeric(h_clean)
  if (!is.na(height) && height > 10) return(11)
  return(height)
}

scans$height_clean <- sapply(scans[[scan_height_col]], clean_height)

parse_date_gps_informed <- function(date_str, gps_dates_ref) {
  date_str <- trimws(as.character(date_str))
  if (grepl("^\\d+$", date_str)) {
    r_numeric <- as.numeric(date_str)
    date <- as.Date(r_numeric, origin = "1970-01-01")
    year <- as.integer(format(date, "%Y"))
    if (!is.na(date) && year >= 2014 && year <= 2019) return(date)
    return(NA)
  }
  if (grepl("/", date_str)) {
    parts <- strsplit(date_str, "/")[[1]]
    if (length(parts) != 3) return(NA)
    part1 <- as.integer(parts[1])
    part2 <- as.integer(parts[2])
    year <- as.integer(parts[3])
    if (part1 > 12) return(as.Date(paste(year, part2, part1, sep = "-")))
    if (part2 > 12) return(as.Date(paste(year, part1, part2, sep = "-")))
    uk_date <- tryCatch(as.Date(paste(year, part2, part1, sep = "-")), error = function(e) NA)
    us_date <- tryCatch(as.Date(paste(year, part1, part2, sep = "-")), error = function(e) NA)
    uk_in_gps <- !is.na(uk_date) && (uk_date %in% gps_dates_ref)
    us_in_gps <- !is.na(us_date) && (us_date %in% gps_dates_ref)
    if (uk_in_gps && !us_in_gps) return(uk_date)
    if (us_in_gps && !uk_in_gps) return(us_date)
    if (!is.na(uk_date)) return(uk_date)
    return(us_date)
  }
  return(NA)
}

scans$date_clean <- sapply(scans[[scan_date_col]], parse_date_gps_informed,
                           gps_dates_ref = gps_dates)
scans$date_clean <- as.Date(scans$date_clean, origin = "1970-01-01")

parse_time_robust <- function(time_str) {
  time_str <- trimws(as.character(time_str))
  if (grepl("^\\d{1,2}:\\d{2}:\\d{2}$", time_str)) {
    parts <- strsplit(time_str, ":")[[1]]
    return(sprintf("%02d:%02d:%02d", as.integer(parts[1]),
                   as.integer(parts[2]), as.integer(parts[3])))
  }
  if (grepl("^\\d{1,2}:\\d{2}$", time_str)) {
    parts <- strsplit(time_str, ":")[[1]]
    return(sprintf("%02d:%02d:00", as.integer(parts[1]), as.integer(parts[2])))
  }
  return(NA)
}

scans$time_clean <- sapply(scans[[scan_time_col]], parse_time_robust)
scans$datetime <- as.POSIXct(paste(scans$date_clean, scans$time_clean),
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

scan_samples <- scans %>%
  filter(!is.na(datetime)) %>%
  group_by(date_clean, datetime) %>%
  summarize(
    n_individuals = n(),
    n_valid_heights = sum(!is.na(height_clean)),
    max_height_m = ifelse(n_valid_heights > 0,
                          max(height_clean, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(max_height_m)) %>%
  arrange(datetime)

# Build SSF steps
ver_track <- vervet_clean %>%
  make_track(longitude, latitude, timestamp, crs = 4326, all_cols = TRUE) %>%
  transform_coords(32735)

ver_steps <- ver_track %>%
  steps() %>%
  mutate(
    dt_hours = as.numeric(difftime(t2_, t1_, units = "hours")),
    date1 = as.Date(t1_, tz = "UTC"),
    date2 = as.Date(t2_, tz = "UTC"),
    same_date = (date1 == date2)
  ) %>%
  filter(same_date == TRUE, dt_hours <= 2.0) %>%
  dplyr::select(-date1, -date2, -same_date)

ver_ssf_data <- ver_steps %>% random_steps(n_control = 10)

# Match observer heights from scan data
strata_info <- ver_ssf_data %>%
  filter(case_ == TRUE) %>%
  dplyr::select(step_id_, t2_)

sentinel_heights <- numeric(nrow(strata_info))
height_source <- character(nrow(strata_info))
time_diff_hours <- numeric(nrow(strata_info))

for (i in 1:nrow(strata_info)) {
  step_time <- strata_info$t2_[i]
  step_date <- as.Date(step_time, tz = "UTC")
  same_day_scans <- scan_samples %>% filter(date_clean == step_date)
  if (nrow(same_day_scans) > 0) {
    time_diffs <- abs(difftime(same_day_scans$datetime, step_time, units = "hours"))
    nearest_idx <- which.min(time_diffs)
    sentinel_heights[i] <- same_day_scans$max_height_m[nearest_idx]
    time_diff_hours[i] <- as.numeric(time_diffs[nearest_idx])
    height_source[i] <- "scan_same_day"
  } else {
    sentinel_heights[i] <- NA
    time_diff_hours[i] <- NA
    height_source[i] <- "no_scan_on_date"
  }
}

strata_info$sentinel_height_m <- sentinel_heights
strata_info$height_source <- height_source
strata_info$time_diff_hours <- time_diff_hours

ver_ssf_data <- ver_ssf_data %>%
  left_join(strata_info %>% dplyr::select(step_id_, sentinel_height_m,
                                          height_source, time_diff_hours),
            by = "step_id_")

valid_strata <- strata_info %>%
  filter(!is.na(sentinel_height_m)) %>%
  pull(step_id_)

ver_ssf_data <- ver_ssf_data %>%
  filter(step_id_ %in% valid_strata)

# Pre-match baboon fixes with vervets
bab_times_num <- as.numeric(baboon_clean$timestamp)
step_times_num <- as.numeric(ver_ssf_data$t2_)

bab_order <- order(bab_times_num)
bab_sorted <- bab_times_num[bab_order]

nearest_bab <- function(st) {
  idx <- findInterval(st, bab_sorted)
  idx <- max(1, min(idx, length(bab_sorted)))
  candidates <- unique(pmax(1, pmin(c(idx, idx + 1), length(bab_sorted))))
  diffs <- abs(st - bab_sorted[candidates])
  best <- candidates[which.min(diffs)]
  bab_order[best]
}

bab_match_idx <- sapply(step_times_num, nearest_bab)
bab_time_gap  <- abs(step_times_num - bab_times_num[bab_match_idx]) / 60

ver_ssf_data$baboon_x <- baboon_clean$x_utm[bab_match_idx]
ver_ssf_data$baboon_y <- baboon_clean$y_utm[bab_match_idx]
ver_ssf_data$bab_time_gap_min <- bab_time_gap
ver_ssf_data$bab_valid <- ver_ssf_data$bab_time_gap_min <= 10

ver_ssf_data$dist_to_baboon <- sqrt(
  (ver_ssf_data$x2_ - ver_ssf_data$baboon_x)^2 +
    (ver_ssf_data$y2_ - ver_ssf_data$baboon_y)^2
)
ver_ssf_data$dist_to_baboon[!ver_ssf_data$bab_valid] <- NA

# Sector viewshed with canopy permeability
# 9 rays across a 45-degree cone centred on the bearing to the baboon group
# prop_visible = fraction of rays with clearance >= 0 (primary metric)
# Sights/sounds assumed mostly blocked beyond 1km .

calculate_sector_viewshed <- function(observer_x, observer_y, observer_height_agl,
                                      target_x, target_y,
                                      dem, canopy_height_raster, cc_raster,
                                      target_height_agl = 1.5,
                                      near_exclusion_m = 20,
                                      sector_degrees = 45,
                                      n_rays = 9) {
  
  na_result <- list(prop_visible = NA, max_clearance = NA,
                    mean_clearance = NA, distance_m = NA)
  
  if (any(is.na(c(observer_x, observer_y, observer_height_agl, target_x, target_y)))) {
    return(na_result)
  }
  
  distance <- sqrt((target_x - observer_x)^2 + (target_y - observer_y)^2)
  
  if (is.na(distance) || distance < 1) {
    return(list(prop_visible = 1, max_clearance = observer_height_agl,
                mean_clearance = observer_height_agl, distance_m = distance))
  }
  
  obs_ground <- raster::extract(dem, matrix(c(observer_x, observer_y), ncol = 2))
  if (is.na(obs_ground)) return(na_result)
  observer_elev <- obs_ground + observer_height_agl
  
  dx <- target_x - observer_x
  dy <- target_y - observer_y
  central_bearing <- atan2(dx, dy)
  
  half_sector <- (sector_degrees / 2) * pi / 180
  angles <- seq(-half_sector, half_sector, length.out = n_rays)
  bearings <- central_bearing + angles
  
  ray_clearances <- numeric(n_rays)
  
  for (r in seq_along(bearings)) {
    ray_x <- observer_x + distance * sin(bearings[r])
    ray_y <- observer_y + distance * cos(bearings[r])
    
    ray_ground <- raster::extract(dem, matrix(c(ray_x, ray_y), ncol = 2))
    if (is.na(ray_ground)) {
      ray_clearances[r] <- NA
      next
    }
    ray_elev <- ray_ground + target_height_agl
    
    n_points <- max(5, ceiling(distance / 50))
    x_points <- seq(observer_x, ray_x, length.out = n_points)
    y_points <- seq(observer_y, ray_y, length.out = n_points)
    
    terrain_elevs  <- raster::extract(dem, cbind(x_points, y_points))
    canopy_heights <- raster::extract(canopy_height_raster, cbind(x_points, y_points))
    canopy_covers  <- raster::extract(cc_raster, cbind(x_points, y_points))
    
    canopy_heights[is.na(canopy_heights)] <- 0
    canopy_covers[is.na(canopy_covers)]   <- 0
    
    cover_fraction <- canopy_covers / 100
    obstacle_elevs <- terrain_elevs + canopy_heights * cover_fraction
    
    fractions <- seq(0, 1, length.out = n_points)
    los_elevs <- observer_elev + fractions * (ray_elev - observer_elev)
    clearances <- los_elevs - obstacle_elevs
    
    point_distances <- fractions * distance
    valid <- point_distances > near_exclusion_m &
      point_distances < (distance - 5)
    
    if (sum(valid) == 0) {
      ray_clearances[r] <- observer_height_agl
    } else {
      ray_clearances[r] <- min(clearances[valid], na.rm = TRUE)
    }
  }
  
  valid_rays <- ray_clearances[!is.na(ray_clearances)]
  if (length(valid_rays) == 0) return(na_result)
  
  return(list(
    prop_visible   = sum(valid_rays >= 0) / length(valid_rays),
    max_clearance  = max(valid_rays),
    mean_clearance = mean(valid_rays),
    distance_m     = distance
  ))
}

# Compute sector viewshed covariates
# Pairs within 1km get full sector calculation, >1km assumed fully blocked

compute_idx <- which(ver_ssf_data$bab_valid)

ver_ssf_data$viewshed_prop_visible   <- NA
ver_ssf_data$viewshed_max_clearance  <- NA
ver_ssf_data$viewshed_mean_clearance <- NA

start_time <- Sys.time()

for (k in seq_along(compute_idx)) {
  i   <- compute_idx[k]
  row <- ver_ssf_data[i, ]
  dist <- row$dist_to_baboon
  
  if (!is.na(dist) && dist > 1000) {
    ver_ssf_data$viewshed_prop_visible[i]   <- 0
    ver_ssf_data$viewshed_max_clearance[i]  <- -99
    ver_ssf_data$viewshed_mean_clearance[i] <- -99
    next
  }
  
  vs <- calculate_sector_viewshed(
    row$x2_, row$y2_, row$sentinel_height_m,
    row$baboon_x, row$baboon_y,
    dem, canopy_height_raster, cc_raster
  )
  
  ver_ssf_data$viewshed_prop_visible[i]   <- vs$prop_visible
  ver_ssf_data$viewshed_max_clearance[i]  <- vs$max_clearance
  ver_ssf_data$viewshed_mean_clearance[i] <- vs$mean_clearance
  
  if (k %% 2000 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- k / elapsed
    remaining <- (length(compute_idx) - k) / rate
    cat(sprintf("  %d/%d (%.0f%%) - %.0f rows/sec - ETA: %.0f min\n",
                k, length(compute_idx), 100 * k / length(compute_idx),
                rate, remaining / 60))
  }
}

# Extract horizontal visibility
hv_raster <- raster("horizontal_visibility.tif")
coords_matrix <- cbind(ver_ssf_data$x2_, ver_ssf_data$y2_)
ver_ssf_data$horizontal_visibility_m <- raster::extract(hv_raster, coords_matrix)
ver_ssf_data$canopy_cover_pct        <- raster::extract(cc_raster, coords_matrix)

# Primary viewshed metric: proportion of sector visible toward baboon
ver_ssf_data$viewshed_clearance_m <- ver_ssf_data$viewshed_prop_visible

write.csv(ver_ssf_data, "ver_ssf_data_viewshed_sector.csv", row.names = FALSE)

