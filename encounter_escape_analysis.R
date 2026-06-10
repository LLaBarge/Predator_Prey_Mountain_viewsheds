# Encounter escape analysis
# Do vervets with clearer views of baboons escape perpendicular (monitoring
# baboon trajectory) while those relying on auditory cues flee directly away?
#
# Predictions
#   High viewshed → escape angle ~90° (perpendicular, can see baboon trajectory)
#   Low viewshed  → escape angle ~0° (fleeing ahead, auditory cues only)
#
# Also tests whether vervets move into denser cover after detection,
# and whether this differs with viewshed clarity.

library(dplyr)
library(sf)
library(raster)
library(ggplot2)

library(brms)
library(bayestestR)
library(tidyr)

dem <- raster("dem.tif")
canopy_height_raster <- crop(raster("meta_tree_height.tif"), dem)
cc_raster <- raster("canopy_cover.tif")
hv_raster <- raster("horizontal_visibility.tif")

# Load and prepare GPS data
vervet_clean <- readRDS("vervet_clean.rds")
baboon_clean <- readRDS("baboon_clean.rds")

vervet_clean$timestamp <- as.POSIXct(vervet_clean$New_Timestamp, tz = "UTC")
baboon_clean$timestamp <- as.POSIXct(baboon_clean$New_Timestamp, tz = "UTC")

vervet_sf <- st_as_sf(vervet_clean, coords = c("longitude", "latitude"), crs = 4326)
vervet_utm <- st_transform(vervet_sf, crs = 32735)
coords_v <- st_coordinates(vervet_utm)
vervet_clean$x_utm <- coords_v[, 1]
vervet_clean$y_utm <- coords_v[, 2]

baboon_sf <- st_as_sf(baboon_clean, coords = c("longitude", "latitude"), crs = 4326)
baboon_utm <- st_transform(baboon_sf, crs = 32735)
coords_b <- st_coordinates(baboon_utm)
baboon_clean$x_utm <- coords_b[, 1]
baboon_clean$y_utm <- coords_b[, 2]

vervet_clean <- vervet_clean %>% arrange(timestamp)
baboon_clean <- baboon_clean %>% arrange(timestamp)

# Load encounter data for vervet heights
encounters_ref <- read.csv("encounter_behavioral_analysis.csv", stringsAsFactors = FALSE)
encounters_ref$start_time <- as.POSIXct(encounters_ref$start_time, tz = "UTC")

# Match each vervet fix to nearest baboon fix within 30 minutes
proximity_threshold <- 350

vervet_clean$bab_idx <- NA
vervet_clean$dist_to_bab <- NA

bab_ts <- as.numeric(baboon_clean$timestamp)

for (i in 1:nrow(vervet_clean)) {
  ver_t <- as.numeric(vervet_clean$timestamp[i])
  diffs <- abs(ver_t - bab_ts)
  best <- which.min(diffs)
  gap_min <- diffs[best] / 60
  
  if (gap_min <= 30) {
    d <- sqrt((vervet_clean$x_utm[i] - baboon_clean$x_utm[best])^2 +
                (vervet_clean$y_utm[i] - baboon_clean$y_utm[best])^2)
    vervet_clean$bab_idx[i] <- best
    vervet_clean$dist_to_bab[i] <- d
  }
}

print(paste("Fixes with baboon match:", sum(!is.na(vervet_clean$dist_to_bab))))
print(paste("Fixes within 350m:", sum(vervet_clean$dist_to_bab <= proximity_threshold, na.rm = TRUE)))

# Identify encounter bouts
vervet_clean$date <- as.Date(vervet_clean$timestamp, tz = "UTC")
vervet_clean$in_proximity <- !is.na(vervet_clean$dist_to_bab) &
  vervet_clean$dist_to_bab <= proximity_threshold

prox_fixes <- which(vervet_clean$in_proximity)

vervet_clean$encounter_id <- NA
enc_id <- 0

if (length(prox_fixes) > 0) {
  enc_id <- 1
  vervet_clean$encounter_id[prox_fixes[1]] <- enc_id
  
  for (j in 2:length(prox_fixes)) {
    prev <- prox_fixes[j - 1]
    curr <- prox_fixes[j]
    
    same_day <- vervet_clean$date[curr] == vervet_clean$date[prev]
    time_gap <- as.numeric(difftime(vervet_clean$timestamp[curr],
                                    vervet_clean$timestamp[prev], units = "hours"))
    
    if (!same_day || is.na(time_gap) || time_gap > 2) {
      enc_id <- enc_id + 1
    }
    vervet_clean$encounter_id[curr] <- enc_id
  }
}

print(paste("Total proximity encounters:", enc_id))

# Sector viewshed with canopy permeability (45 degrees, 9 rays)
# Canopy scaled by cover fraction: dense forest blocks, open woodland transparent
viewshed_sector <- function(obs_x, obs_y, obs_h, tgt_x, tgt_y, tgt_h,
                            dem, canopy_height_raster, cc_raster,
                            sector_deg = 45, n_rays = 9) {
  if (any(is.na(c(obs_x, obs_y, tgt_x, tgt_y)))) return(NA)
  obs_g <- raster::extract(dem, matrix(c(obs_x, obs_y), ncol = 2))
  if (is.na(obs_g)) return(NA)
  d <- sqrt((tgt_x - obs_x)^2 + (tgt_y - obs_y)^2)
  if (d < 1) return(1)
  obs_elev <- obs_g + obs_h
  central <- atan2(tgt_x - obs_x, tgt_y - obs_y)
  half <- (sector_deg / 2) * pi / 180
  bearings <- central + seq(-half, half, length.out = n_rays)
  clear <- 0; total <- 0
  for (r in seq_along(bearings)) {
    rx <- obs_x + d * sin(bearings[r])
    ry <- obs_y + d * cos(bearings[r])
    rg <- raster::extract(dem, matrix(c(rx, ry), ncol = 2))
    if (is.na(rg)) next
    n_pts <- max(5, ceiling(d / 30))
    xs <- seq(obs_x, rx, length.out = n_pts)
    ys <- seq(obs_y, ry, length.out = n_pts)
    ter <- raster::extract(dem, cbind(xs, ys))
    ch  <- raster::extract(canopy_height_raster, cbind(xs, ys))
    cc  <- raster::extract(cc_raster, cbind(xs, ys))
    ch[is.na(ch)] <- 0; cc[is.na(cc)] <- 0
    obstacle <- ter + ch * (cc / 100)
    frac <- seq(0, 1, length.out = n_pts)
    los <- obs_elev + frac * ((rg + tgt_h) - obs_elev)
    cl <- los[2:(n_pts-1)] - obstacle[2:(n_pts-1)]
    total <- total + 1
    if (min(cl, na.rm = TRUE) >= 0) clear <- clear + 1
  }
  if (total == 0) return(NA)
  clear / total
}

# Get vervet height from encounter reference data
get_vervet_height <- function(fix_timestamp, fix_date, encounters_ref) {
  same_day <- encounters_ref %>% filter(as.Date(start_time, tz = "UTC") == fix_date)
  if (nrow(same_day) == 0) return(6)
  time_diffs <- abs(difftime(same_day$start_time, fix_timestamp, units = "mins"))
  nearest <- which.min(time_diffs)
  h <- same_day$height_before_m[nearest]
  if (is.na(h)) h <- same_day$height_after_m[nearest]
  if (is.na(h)) h <- 6
  return(h)
}

# For each encounter, find closest approach and analyse escape movement
results <- data.frame()

for (eid in 1:enc_id) {
  fixes <- vervet_clean %>% filter(encounter_id == eid) %>% arrange(timestamp)
  if (nrow(fixes) < 2) next
  
  closest_idx <- which.min(fixes$dist_to_bab)
  if (is.na(closest_idx)) next
  if (closest_idx >= nrow(fixes)) next
  
  bab_idx_at_closest <- fixes$bab_idx[closest_idx]
  if (is.na(bab_idx_at_closest) || bab_idx_at_closest <= 1) next
  
  # Positions at closest approach
  vx0 <- fixes$x_utm[closest_idx]
  vy0 <- fixes$y_utm[closest_idx]
  bx0 <- baboon_clean$x_utm[bab_idx_at_closest]
  by0 <- baboon_clean$y_utm[bab_idx_at_closest]
  
  # Vervet position one step after closest approach
  vx1 <- fixes$x_utm[closest_idx + 1]
  vy1 <- fixes$y_utm[closest_idx + 1]
  
  # Baboon position one step before closest approach
  bx_prev <- baboon_clean$x_utm[bab_idx_at_closest - 1]
  by_prev <- baboon_clean$y_utm[bab_idx_at_closest - 1]
  
  if (any(is.na(c(vx0, vy0, vx1, vy1, bx0, by0, bx_prev, by_prev)))) next
  
  # Baboon approach bearing
  bab_approach <- atan2(bx0 - bx_prev, by0 - by_prev) * 180 / pi
  
  # Vervet escape bearing
  ver_escape <- atan2(vx1 - vx0, vy1 - vy0) * 180 / pi
  
  # Relative angle
  # 0 = fleeing ahead of baboons, 90 = perpendicular, 180 = toward baboons
  rel_angle <- abs(ver_escape - bab_approach) %% 360
  if (rel_angle > 180) rel_angle <- 360 - rel_angle
  
  # Vervet step length
  ver_step <- sqrt((vx1 - vx0)^2 + (vy1 - vy0)^2)
  
  # Habitat at origin and destination
  cc_origin <- raster::extract(cc_raster, matrix(c(vx0, vy0), ncol = 2))
  cc_dest   <- raster::extract(cc_raster, matrix(c(vx1, vy1), ncol = 2))
  hv_origin <- raster::extract(hv_raster, matrix(c(vx0, vy0), ncol = 2))
  hv_dest   <- raster::extract(hv_raster, matrix(c(vx1, vy1), ncol = 2))
  
  # Vervet height from encounter reference data
  obs_height <- get_vervet_height(fixes$timestamp[closest_idx],
                                  fixes$date[closest_idx],
                                  encounters_ref)
  
  # Viewshed at closest approach (continuous proportion visible)
  vs <- viewshed_sector(vx0, vy0, obs_height, bx0, by0, 1.5,
                        dem, canopy_height_raster, cc_raster)
  
  results <- rbind(results, data.frame(
    encounter_id   = eid,
    date           = as.character(fixes$date[closest_idx]),
    min_distance   = fixes$dist_to_bab[closest_idx],
    obs_height_m   = obs_height,
    rel_angle      = rel_angle,
    ver_step_m     = ver_step,
    cc_origin      = cc_origin,
    cc_dest        = cc_dest,
    cc_change      = cc_dest - cc_origin,
    hv_origin      = hv_origin,
    hv_dest        = hv_dest,
    hv_change      = hv_dest - hv_origin,
    viewshed       = vs
  ))
}

print(paste("Encounters with escape data:", nrow(results)))

# Filter to encounters where vervets actually moved (step > 10m)
movers <- results %>% filter(ver_step_m > 10)
print(paste("Encounters with movement >10m:", nrow(movers)))

# Viewshed tertiles for grouping
movers$viewshed_group <- cut(movers$viewshed,
                             breaks = c(-0.05, 0.1, 0.5, 1.01),
                             labels = c("Low visibility", "Moderate visibility", "High visibility"))

# Escape angle by viewshed group
print("Escape angle by viewshed tertile")
print("Prediction: high visibility → ~90° (perpendicular), low → ~0° (fleeing ahead)")
movers %>%
  group_by(viewshed_group) %>%
  summarize(
    n = n(),
    mean_viewshed = mean(viewshed),
    mean_angle = mean(rel_angle),
    median_angle = median(rel_angle),
    sd_angle = sd(rel_angle),
    .groups = "drop"
  ) %>% print()

# Habitat change by viewshed group
print("Habitat change by viewshed tertile")
movers %>%
  group_by(viewshed_group) %>%
  summarize(
    n = n(),
    mean_cc_change = mean(cc_change, na.rm = TRUE),
    mean_hv_change = mean(hv_change, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

# Escape angle vs viewshed (continuous)
p_scatter <- ggplot(movers, aes(x = viewshed, y = rel_angle)) +
  geom_point(size = 2.5, alpha = 0.5, colour = "#2E86AB") +
  geom_smooth(method = "lm", colour = "#D32F2F", fill = "#D32F2F", alpha = 0.2) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "grey30") +
  annotate("text", x = 0, y = 93, label = "Perpendicular", hjust = 0,
           size = 3, colour = "grey30") +
  labs(x = "Sector viewshed (proportion visible toward baboons)",
       y = "Escape angle relative to baboon approach (°)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("escape_angle_vs_viewshed.png", p_scatter, width = 8, height = 5, dpi = 300)
p_scatter


write.csv(results, "encounter_escape_analysis.csv", row.names = FALSE)


# Reverse viewshed analysis: do vervets move to positions that are
# (a) more concealed FROM baboons, and/or
# (b) still allow monitoring OF baboons?
#
# For each encounter, calculate:
#   baboon → vervet viewshed BEFORE escape (how visible was the vervet?)
#   baboon → vervet viewshed AFTER escape (how visible is the vervet now?)
#   vervet → baboon viewshed AFTER escape (can the vervet still see baboons?)

results$vs_bab_to_ver_before <- NA  # baboon's view of vervet at t0
results$vs_bab_to_ver_after  <- NA  # baboon's view of vervet at t1
results$vs_ver_to_bab_after  <- NA  # vervet's view of baboon from new position

for (eid in 1:enc_id) {
  row_idx <- which(results$encounter_id == eid)
  if (length(row_idx) == 0) next
  
  fixes <- vervet_clean %>% filter(encounter_id == eid) %>% arrange(timestamp)
  if (nrow(fixes) < 2) next
  
  closest_idx <- which.min(fixes$dist_to_bab)
  if (is.na(closest_idx) || closest_idx >= nrow(fixes)) next
  
  bab_idx_at_closest <- fixes$bab_idx[closest_idx]
  if (is.na(bab_idx_at_closest)) next
  
  vx0 <- fixes$x_utm[closest_idx]
  vy0 <- fixes$y_utm[closest_idx]
  vx1 <- fixes$x_utm[closest_idx + 1]
  vy1 <- fixes$y_utm[closest_idx + 1]
  bx0 <- baboon_clean$x_utm[bab_idx_at_closest]
  by0 <- baboon_clean$y_utm[bab_idx_at_closest]
  
  if (any(is.na(c(vx0, vy0, vx1, vy1, bx0, by0)))) next
  
  obs_height <- results$obs_height_m[row_idx]
  
  # Baboon's view of vervet BEFORE escape (baboon at 1.5m looking at vervet at tree height)
  results$vs_bab_to_ver_before[row_idx] <- viewshed_sector(
    bx0, by0, 1.5, vx0, vy0, obs_height,
    dem, canopy_height_raster, cc_raster
  )
  
  # Baboon's view of vervet AFTER escape (baboon at 1.5m looking at vervet's new position)
  results$vs_bab_to_ver_after[row_idx] <- viewshed_sector(
    bx0, by0, 1.5, vx1, vy1, obs_height,
    dem, canopy_height_raster, cc_raster
  )
  
  # Vervet's view of baboon FROM new position (vervet at tree height looking back)
  results$vs_ver_to_bab_after[row_idx] <- viewshed_sector(
    vx1, vy1, obs_height, bx0, by0, 1.5,
    dem, canopy_height_raster, cc_raster
  )
}

# Change in concealment and surveillance
results$concealment_change <- results$vs_bab_to_ver_before - results$vs_bab_to_ver_after
results$surveillance_after <- results$vs_ver_to_bab_after

# Positive concealment_change = vervet moved to where baboons see it LESS
# High surveillance_after = vervet can still see baboons from new position

movers <- results %>% filter(ver_step_m > 10) # 10 chosen b/c gps error might result in 5m

# Summaries by viewshed tertile
movers$viewshed_group <- dplyr::ntile(movers$viewshed, 3)
movers$viewshed_group <- factor(movers$viewshed_group,
                                labels = c("Low visibility", "Moderate visibility", "High visibility"))

print("Concealment and surveillance by viewshed at encounter")
movers %>%
  group_by(viewshed_group) %>%
  summarize(
    n = n(),
    mean_conceal_change = mean(concealment_change, na.rm = TRUE),
    mean_surveillance   = mean(surveillance_after, na.rm = TRUE),
    mean_bab_view_before = mean(vs_bab_to_ver_before, na.rm = TRUE),
    mean_bab_view_after  = mean(vs_bab_to_ver_after, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

# Overall: did vervets increase concealment?
print("Overall concealment change (positive = more hidden after escape)")
print(paste("Mean:", round(mean(movers$concealment_change, na.rm = TRUE), 3)))

# Overall: can vervets still see baboons from new position?
print("Surveillance from new position")
print(paste("Mean prop visible:", round(mean(movers$surveillance_after, na.rm = TRUE), 3)))

# Plot: concealment vs surveillance (the strategic tradeoff)
p_tradeoff <- ggplot(movers, aes(x = surveillance_after, y = concealment_change,
                                 colour = viewshed_group)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Low visibility" = "#D32F2F",
                                 "Moderate visibility" = "#FF8F00",
                                 "High visibility" = "#1565C0"),
                      name = "Viewshed at encounter") +
  annotate("text", x = 0.75, y = max(movers$concealment_change, na.rm = TRUE),
           label = "Hidden but\ncan monitor", size = 3, colour = "grey30") +
  annotate("text", x = 0.25, y = max(movers$concealment_change, na.rm = TRUE),
           label = "Hidden and\nblind", size = 3, colour = "grey30") +
  labs(x = "Surveillance after escape (vervet's view of baboons)",
       y = "Concealment change (positive = more hidden)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("concealment_surveillance_tradeoff.png", p_tradeoff, width = 8, height = 6, dpi = 300)
p_tradeoff

# end exploratory plots--------------------------------------------------------
# Encounter escape analysis, Bayesian paired Beta model.
#
# Before: vervet's view of baboon from vervet's pre-encounter position
# After:  vervet's view of baboon from vervet's post-encounter position


# Load and prepare GPS data
vervet_clean <- readRDS("vervet_clean.rds")
baboon_clean <- readRDS("baboon_clean.rds")

vervet_clean$timestamp <- as.POSIXct(vervet_clean$New_Timestamp, tz = "UTC")
baboon_clean$timestamp <- as.POSIXct(baboon_clean$New_Timestamp, tz = "UTC")

vervet_sf <- st_as_sf(vervet_clean, coords = c("longitude", "latitude"), crs = 4326)
vervet_utm <- st_transform(vervet_sf, crs = 32735)
coords_v <- st_coordinates(vervet_utm)
vervet_clean$x_utm <- coords_v[, 1]
vervet_clean$y_utm <- coords_v[, 2]

baboon_sf <- st_as_sf(baboon_clean, coords = c("longitude", "latitude"), crs = 4326)
baboon_utm <- st_transform(baboon_sf, crs = 32735)
coords_b <- st_coordinates(baboon_utm)
baboon_clean$x_utm <- coords_b[, 1]
baboon_clean$y_utm <- coords_b[, 2]

vervet_clean <- vervet_clean %>% arrange(timestamp)
baboon_clean <- baboon_clean %>% arrange(timestamp)
vervet_clean$date <- as.Date(vervet_clean$timestamp, tz = "UTC")

# Load scan data for vervet heights
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

gps_dates <- unique(vervet_clean$date)

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
    n_valid_heights = sum(!is.na(height_clean)),
    max_height_m = ifelse(n_valid_heights > 0,
                          max(height_clean, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(max_height_m)) %>%
  arrange(datetime)

get_height_at_time <- function(ts, ts_date) {
  same_day <- scan_samples %>% filter(date_clean == ts_date)
  if (nrow(same_day) == 0) return(NA)
  time_diffs <- abs(difftime(same_day$datetime, ts, units = "hours"))
  nearest <- which.min(time_diffs)
  return(same_day$max_height_m[nearest])
}

# Match each vervet fix to nearest baboon fix within 30 minutes
proximity_threshold <- 350

vervet_clean$bab_idx <- NA
vervet_clean$dist_to_bab <- NA

bab_ts <- as.numeric(baboon_clean$timestamp)

for (i in 1:nrow(vervet_clean)) {
  ver_t <- as.numeric(vervet_clean$timestamp[i])
  diffs <- abs(ver_t - bab_ts)
  best <- which.min(diffs)
  gap_min <- diffs[best] / 60
  
  if (gap_min <= 30) {
    d <- sqrt((vervet_clean$x_utm[i] - baboon_clean$x_utm[best])^2 +
                (vervet_clean$y_utm[i] - baboon_clean$y_utm[best])^2)
    vervet_clean$bab_idx[i] <- best
    vervet_clean$dist_to_bab[i] <- d
  }
}

print(paste("Fixes with baboon match:", sum(!is.na(vervet_clean$dist_to_bab))))
print(paste("Fixes within 350m:", sum(vervet_clean$dist_to_bab <= proximity_threshold, na.rm = TRUE)))

# Identify encounter bouts
vervet_clean$in_proximity <- !is.na(vervet_clean$dist_to_bab) &
  vervet_clean$dist_to_bab <= proximity_threshold

prox_fixes <- which(vervet_clean$in_proximity)

vervet_clean$encounter_id <- NA
enc_id <- 0

if (length(prox_fixes) > 0) {
  enc_id <- 1
  vervet_clean$encounter_id[prox_fixes[1]] <- enc_id
  
  for (j in 2:length(prox_fixes)) {
    prev <- prox_fixes[j - 1]
    curr <- prox_fixes[j]
    
    same_day <- vervet_clean$date[curr] == vervet_clean$date[prev]
    time_gap <- as.numeric(difftime(vervet_clean$timestamp[curr],
                                    vervet_clean$timestamp[prev], units = "hours"))
    
    if (!same_day || is.na(time_gap) || time_gap > 2) {
      enc_id <- enc_id + 1
    }
    vervet_clean$encounter_id[curr] <- enc_id
  }
}

print(paste("Total proximity encounters:", enc_id))

# Sector viewshed with canopy permeability (45 degrees, 9 rays)
viewshed_sector <- function(obs_x, obs_y, obs_h, tgt_x, tgt_y, tgt_h,
                            dem, canopy_height_raster, cc_raster,
                            sector_deg = 45, n_rays = 9) {
  if (any(is.na(c(obs_x, obs_y, tgt_x, tgt_y)))) return(NA)
  obs_g <- raster::extract(dem, matrix(c(obs_x, obs_y), ncol = 2))
  if (is.na(obs_g)) return(NA)
  d <- sqrt((tgt_x - obs_x)^2 + (tgt_y - obs_y)^2)
  if (d < 1) return(1)
  obs_elev <- obs_g + obs_h
  central <- atan2(tgt_x - obs_x, tgt_y - obs_y)
  half <- (sector_deg / 2) * pi / 180
  bearings <- central + seq(-half, half, length.out = n_rays)
  clear <- 0; total <- 0
  for (r in seq_along(bearings)) {
    rx <- obs_x + d * sin(bearings[r])
    ry <- obs_y + d * cos(bearings[r])
    rg <- raster::extract(dem, matrix(c(rx, ry), ncol = 2))
    if (is.na(rg)) next
    n_pts <- max(5, ceiling(d / 30))
    xs <- seq(obs_x, rx, length.out = n_pts)
    ys <- seq(obs_y, ry, length.out = n_pts)
    ter <- raster::extract(dem, cbind(xs, ys))
    ch  <- raster::extract(canopy_height_raster, cbind(xs, ys))
    cc  <- raster::extract(cc_raster, cbind(xs, ys))
    ch[is.na(ch)] <- 0; cc[is.na(cc)] <- 0
    obstacle <- ter + ch * (cc / 100)
    frac <- seq(0, 1, length.out = n_pts)
    los <- obs_elev + frac * ((rg + tgt_h) - obs_elev)
    cl <- los[2:(n_pts-1)] - obstacle[2:(n_pts-1)]
    total <- total + 1
    if (min(cl, na.rm = TRUE) >= 0) clear <- clear + 1
  }
  if (total == 0) return(NA)
  clear / total
}

# For each encounter, find before and after fixes and compute viewsheds
# Viewshed: vervet looking toward baboon (vervet at scan height, baboon at 1.5m)
results <- data.frame()

for (eid in 1:enc_id) {
  enc_rows <- which(vervet_clean$encounter_id == eid)
  if (length(enc_rows) == 0) next
  
  first_prox <- min(enc_rows)
  last_prox  <- max(enc_rows)
  enc_date   <- vervet_clean$date[first_prox]
  
  # Before fix: row immediately before the first proximity fix, same day
  if (first_prox <= 1) next
  before_idx <- first_prox - 1
  if (vervet_clean$date[before_idx] != enc_date) next
  
  # After fix: row after the last proximity fix, within 40 min same day
  after_idx <- NA
  if (last_prox < nrow(vervet_clean)) {
    candidate <- last_prox + 1
    if (vervet_clean$date[candidate] == enc_date) {
      gap <- as.numeric(difftime(vervet_clean$timestamp[candidate],
                                 vervet_clean$timestamp[last_prox], units = "mins"))
      if (!is.na(gap) && gap <= 40) after_idx <- candidate
    }
    if (is.na(after_idx) && (last_prox + 2) <= nrow(vervet_clean)) {
      candidate2 <- last_prox + 2
      if (vervet_clean$date[candidate2] == enc_date) {
        gap2 <- as.numeric(difftime(vervet_clean$timestamp[candidate2],
                                    vervet_clean$timestamp[last_prox], units = "mins"))
        if (!is.na(gap2) && gap2 <= 40) after_idx <- candidate2
      }
    }
  }
  if (is.na(after_idx)) next
  
  # Vervet positions
  vx_before <- vervet_clean$x_utm[before_idx]
  vy_before <- vervet_clean$y_utm[before_idx]
  vx_after  <- vervet_clean$x_utm[after_idx]
  vy_after  <- vervet_clean$y_utm[after_idx]
  
  # Baboon positions at before and after timestamps
  bab_before_idx <- vervet_clean$bab_idx[before_idx]
  bab_after_idx  <- vervet_clean$bab_idx[after_idx]
  if (is.na(bab_before_idx) || is.na(bab_after_idx)) next
  
  bx_before <- baboon_clean$x_utm[bab_before_idx]
  by_before <- baboon_clean$y_utm[bab_before_idx]
  bx_after  <- baboon_clean$x_utm[bab_after_idx]
  by_after  <- baboon_clean$y_utm[bab_after_idx]
  
  if (any(is.na(c(vx_before, vy_before, vx_after, vy_after,
                  bx_before, by_before, bx_after, by_after)))) next
  
  # Vervet heights from scan data
  h_before <- get_height_at_time(vervet_clean$timestamp[before_idx], enc_date)
  h_after  <- get_height_at_time(vervet_clean$timestamp[after_idx], enc_date)
  if (is.na(h_before)) h_before <- 6
  if (is.na(h_after))  h_after  <- 6
  
  # Distances
  dist_before <- sqrt((vx_before - bx_before)^2 + (vy_before - by_before)^2)
  dist_after  <- sqrt((vx_after - bx_after)^2 + (vy_after - by_after)^2)
  min_dist    <- min(vervet_clean$dist_to_bab[enc_rows], na.rm = TRUE)
  
  # Vervet's view of baboon BEFORE encounter
  # Vervet at scan height looking toward baboon at 1.5m
  vs_before <- viewshed_sector(vx_before, vy_before, h_before,
                               bx_before, by_before, 1.5,
                               dem, canopy_height_raster, cc_raster)
  
  # Vervet's view of baboon AFTER encounter
  vs_after <- viewshed_sector(vx_after, vy_after, h_after,
                              bx_after, by_after, 1.5,
                              dem, canopy_height_raster, cc_raster)
  
  # Habitat at before and after positions
  cc_before <- raster::extract(cc_raster, matrix(c(vx_before, vy_before), ncol = 2))
  cc_after  <- raster::extract(cc_raster, matrix(c(vx_after, vy_after), ncol = 2))
  hv_before <- raster::extract(hv_raster, matrix(c(vx_before, vy_before), ncol = 2))
  hv_after  <- raster::extract(hv_raster, matrix(c(vx_after, vy_after), ncol = 2))
  
  results <- rbind(results, data.frame(
    encounter_id     = eid,
    date             = as.character(enc_date),
    n_prox_fixes     = length(enc_rows),
    min_distance     = min_dist,
    dist_before      = dist_before,
    dist_after       = dist_after,
    h_before         = h_before,
    h_after          = h_after,
    vs_ver_before    = vs_before,
    vs_ver_after     = vs_after,
    cc_before        = cc_before,
    cc_after         = cc_after,
    hv_before        = hv_before,
    hv_after         = hv_after
  ))
}

print(paste("Encounters with before/after data:", nrow(results)))
print(paste("  Including single-fix encounters:", sum(results$n_prox_fixes == 1)))

# Descriptive summaries
print("Vervet's view of baboon: before vs after encounter")
print(paste("Mean before:", round(mean(results$vs_ver_before, na.rm = TRUE), 3)))
print(paste("Mean after:", round(mean(results$vs_ver_after, na.rm = TRUE), 3)))

print("Vervet height: before vs after")
print(paste("Mean before:", round(mean(results$h_before, na.rm = TRUE), 2), "m"))
print(paste("Mean after:", round(mean(results$h_after, na.rm = TRUE), 2), "m"))

# Bayesian paired Beta model
results_long <- results %>%
  dplyr::select(encounter_id, vs_ver_before, vs_ver_after) %>%
  filter(!is.na(vs_ver_before), !is.na(vs_ver_after)) %>%
  tidyr::pivot_longer(cols = c(vs_ver_before, vs_ver_after),
               names_to = "timing",
               values_to = "ver_view") %>%
  mutate(timing = ifelse(timing == "vs_ver_before", "Before", "After"),
         timing = factor(timing, levels = c("Before", "After")))

# Beta requires (0,1) exclusive
eps <- 0.001
results_long$ver_view <- pmin(pmax(results_long$ver_view, eps), 1 - eps)

print(paste("Paired observations:", nrow(results_long) / 2))

# Priors
# Intercept: no strong expectation about baseline visibility
# timing (After): centred at 0, vervets could move to positions with
#   better view (monitoring) or worse view (breaking line of sight)
# Random intercept per encounter: paired structure
priors_beta <- c(
  prior(normal(0, 1.5),  class = "Intercept"),
  prior(normal(0, 1),    class = "b", coef = "timingAfter"),
  prior(gamma(2, 0.1),   class = "phi"),
  prior(exponential(1),  class = "sd")
)

m_paired <- brm(
  ver_view ~ timing + (1 | encounter_id),
  family = Beta(),
  data = results_long,
  prior = priors_beta,
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2,
  control = list(adapt_delta = 0.95)
)

print(summary(m_paired))

# Posterior of the timing effect
post_timing <- as_draws_df(m_paired)$b_timingAfter

# HDI and ROPE
rope_range <- c(-0.1, 0.1)

print("=== Vervet's view of baboon: Before vs After encounter ===")
print(paste("Posterior mean (logit):", round(mean(post_timing), 4)))
print(paste("Posterior median (logit):", round(median(post_timing), 4)))
print(hdi(post_timing, ci = 0.95))
print(rope(post_timing, range = rope_range, ci = 1))
print(paste("P(better view after):", round(mean(post_timing > 0), 3)))
print(paste("P(worse view after):", round(mean(post_timing < 0), 3)))

# Conditional effects
p_ce <- plot(conditional_effects(m_paired, effects = "timing"), plot = FALSE)
ggsave("vervet_view_before_after.png", p_ce[[1]], width = 6, height = 5, dpi = 300)

# Posterior density
p_post <- ggplot(data.frame(x = post_timing), aes(x = x)) +
  geom_density(fill = "#2E86AB", alpha = 0.4, colour = "#2E86AB") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = mean(post_timing), colour = "#D32F2F", linewidth = 0.8) +
  annotate("rect", xmin = rope_range[1], xmax = rope_range[2],
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey50") +
  annotate("text", x = 0, y = Inf, label = "ROPE", vjust = 2,
           size = 3, colour = "grey30") +
  labs(x = "Timing effect (logit scale)\n\u2190 Concealment                        Greater View \u2192",
       y = "Posterior density") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("posterior_viewshed_change.png", p_post, width = 7, height = 4.5, dpi = 300)
p_post
write.csv(results, "encounter_before_after.csv", row.names = FALSE)