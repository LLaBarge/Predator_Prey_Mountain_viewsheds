library(readxl)
library(dplyr)
library(tidyr)
library(sf)

# ==============================================================================
# Utility functions
# ==============================================================================

# Parse checkerboard value: must be 0-225 integer, else NA
parse_cb <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.na(x) & (x > 225 | x < 0)] <- NA
  x
}

# Remove prefixes like (0), (00), [0], E(0), S from coordinate strings
strip_coord_prefix <- function(s) {
  s <- gsub("^[ESes]\\(0\\)", "", s)
  s <- gsub("^\\(0+\\)", "", s)
  s <- gsub("^\\(0+\\(", "", s)
  s <- gsub("^\\[0\\]", "", s)
  s <- gsub(",", ".", s)
  s
}

# Remove extra decimal points safely (max 3 attempts, no infinite loop)
fix_double_dots <- function(s) {
  for (i in 1:3) {
    if (!grepl("\\d\\.\\d+\\.", s)) break
    s <- sub("(\\d+\\.\\d+)\\.(\\d+)", "\\1\\2", s)
  }
  s
}

# Clean a single coordinate string -> numeric or NA
clean_coord <- function(s) {
  if (is.na(s) || s == "NULL" || s == "") return(NA_real_)
  s <- trimws(s)
  if (grepl("cliff|could not|missing|image|marsh|thick", s, ignore.case = TRUE)) return(NA_real_)
  s <- strip_coord_prefix(s)
  s <- fix_double_dots(s)
  suppressWarnings(as.numeric(trimws(s)))
}

# Vectorised wrapper
clean_coords <- function(x) sapply(x, clean_coord, USE.NAMES = FALSE)

# Assign lon/lat from two cleaned numeric vectors based on value ranges
# lon ~ 29, lat ~ 23 at this site
assign_lonlat <- function(x_clean, y_clean) {
  lon <- rep(NA_real_, length(x_clean))
  lat <- rep(NA_real_, length(x_clean))
  
  # Case 1: x=lat(~23), y=lon(~29) — most common
  normal <- !is.na(x_clean) & !is.na(y_clean) &
    x_clean > 22 & x_clean < 24 & y_clean > 28 & y_clean < 31
  lat[normal] <- x_clean[normal]
  lon[normal] <- y_clean[normal]
  
  # Case 2: x=lon(~29), y=lat(~23) — swapped (J-series etc)
  swapped <- !is.na(x_clean) & !is.na(y_clean) &
    x_clean > 28 & x_clean < 31 & y_clean > 22 & y_clean < 24
  lat[swapped] <- y_clean[swapped]
  lon[swapped] <- x_clean[swapped]
  
  # Latitude should be negative (Southern Hemisphere)
  lat <- -abs(lat)
  
  data.frame(longitude = lon, latitude = lat)
}

# Check coordinates are within study area
valid_coords <- function(lon, lat) {
  !is.na(lon) & !is.na(lat) &
    lon > 29 & lon < 30 &
    lat < -22.5 & lat > -23.5
}


# ==============================================================================
# Standard checkerboard column definitions
# ==============================================================================

cb_defs <- data.frame(
  pattern = c(
    "North.*5 m.*0 m",  "North.*5 m.*1\\.25",
    "North.*10 m.*0 m", "North.*10 m.*1\\.25",
    "East.*5 m.*0 m",   "East.*5 m.*1\\.25",
    "East.*10 m.*0 m",  "East.*10 m.*1\\.25",
    "South.*5 m.*0 m",  "South.*5 m.*1\\.25",
    "South.*10 m.*0 m", "South.*10 m.*1\\.25",
    "West.*5 m.*0 m",   "West.*5 m.*1\\.25",
    "West.*10 m.*0 m",  "West.*10 m.*1\\.25"
  ),
  name = c(
    "N_5m_0m", "N_5m_125m", "N_10m_0m", "N_10m_125m",
    "E_5m_0m", "E_5m_125m", "E_10m_0m", "E_10m_125m",
    "S_5m_0m", "S_5m_125m", "S_10m_0m", "S_10m_125m",
    "W_5m_0m", "W_5m_125m", "W_10m_0m", "W_10m_125m"
  ),
  stringsAsFactors = FALSE
)

# Spatial offsets from NW corner (metres)
offset_defs <- data.frame(
  cb_ground = c("N_5m_0m", "N_10m_0m", "E_5m_0m", "E_10m_0m",
                "S_5m_0m", "S_10m_0m", "W_5m_0m", "W_10m_0m"),
  cb_chest  = c("N_5m_125m", "N_10m_125m", "E_5m_125m", "E_10m_125m",
                "S_5m_125m", "S_10m_125m", "W_5m_125m", "W_10m_125m"),
  dx = c(0, 0, 5, 10, 0, 0, -5, -10),
  dy = c(5, 10, 0, 0, -5, -10, 0, 0),
  direction = c("N", "N", "E", "E", "S", "S", "W", "W"),
  distance_m = c(5, 10, 5, 10, 5, 10, 5, 10),
  stringsAsFactors = FALSE
)


# ==============================================================================
# Load and clean VegPlot_2016
# ==============================================================================

veg2_raw <- read_excel("VegPlot_2016.xlsx", sheet = 1,
                       col_names = FALSE, col_types = "text")

# Find header row (contains "Date" and "Plot")
header_idx <- which(sapply(1:5, function(r) {
  any(grepl("^Date$", veg2_raw[r, ], ignore.case = TRUE)) &
    any(grepl("Plot", veg2_raw[r, ], ignore.case = TRUE))
}))

cnames <- as.character(veg2_raw[header_idx[1], ])
cnames[is.na(cnames)] <- paste0("V", which(is.na(cnames)))
cnames <- make.unique(cnames, sep = "_")

veg2 <- veg2_raw[(header_idx[1] + 1):nrow(veg2_raw), ]
names(veg2) <- cnames

# Identify columns
cb_cols_2 <- names(veg2)[grepl("m away.*m high", names(veg2), ignore.case = TRUE)]
x2 <- names(veg2)[grepl("^x$", names(veg2), ignore.case = TRUE)][1]
y2 <- names(veg2)[grepl("^y$", names(veg2), ignore.case = TRUE)][1]
id2 <- names(veg2)[grepl("Plot.?ID", names(veg2), ignore.case = TRUE)][1]

# Clean coordinates
veg2$x_clean <- clean_coords(veg2[[x2]])
veg2$y_clean <- clean_coords(veg2[[y2]])
ll2 <- assign_lonlat(veg2$x_clean, veg2$y_clean)
veg2$longitude <- ll2$longitude
veg2$latitude  <- ll2$latitude

# Parse checkerboard
for (cc in cb_cols_2) veg2[[cc]] <- parse_cb(veg2[[cc]])

# Deduplicate: one row per plot location
valid2 <- valid_coords(veg2$longitude, veg2$latitude)

veg2_cb <- veg2 %>%
  filter(valid2) %>%
  group_by(across(all_of(id2)), longitude, latitude) %>%
  summarise(across(all_of(cb_cols_2), ~ first(na.omit(.))), .groups = "drop")

message("2016: ", nrow(veg2_cb), " unique plot-locations")


# ==============================================================================
# Load and clean Veg_Plot_Master_2012_2015 (transposed)
# ==============================================================================

raw <- read_excel("Veg_Plot_Master_2012_2015.xlsx",
                  col_names = FALSE, col_types = "text")
row_labels_A <- raw[[1]]  # section headers: "Checkerboard", "Visibility (m)", etc.
row_labels_B <- raw[[2]]  # variable labels: "x", "y", "North - 5 m away, 0 m high", etc.

find_row_B <- function(pattern) {
  idx <- which(grepl(pattern, row_labels_B, ignore.case = TRUE))
  if (length(idx) == 0) return(NA)
  idx[1]
}

row_x <- find_row_B("^\\s*x\\s*$")
row_y <- find_row_B("^\\s*y\\s*$")
cb_start <- which(grepl("^Checkerboard", row_labels_A, ignore.case = TRUE))[1]

message("2012-2015: row_x=", row_x, " row_y=", row_y, " cb_start=", cb_start)

# Find checkerboard rows (labels in column B, after cb_start)
cb_row_idx <- sapply(cb_defs$pattern, function(pat) {
  candidates <- which(grepl(pat, row_labels_B, ignore.case = TRUE))
  candidates <- candidates[candidates >= cb_start]
  if (length(candidates) == 0) return(NA)
  candidates[1]
})

# Data starts from column 3 (A=section headers, B=variable labels, C onward=plot data)
plot_cols <- 3:ncol(raw)

plots1 <- data.frame(
  x_raw = as.character(unlist(raw[row_x, plot_cols])),
  y_raw = as.character(unlist(raw[row_y, plot_cols])),
  stringsAsFactors = FALSE
)

for (i in seq_along(cb_defs$name)) {
  ri <- cb_row_idx[i]
  if (is.na(ri)) {
    plots1[[cb_defs$name[i]]] <- NA_character_
  } else {
    plots1[[cb_defs$name[i]]] <- as.character(unlist(raw[ri, plot_cols]))
  }
}

# Clean coordinates
plots1$x_clean <- clean_coords(plots1$x_raw)
plots1$y_clean <- clean_coords(plots1$y_raw)

message("2012-2015 x non-NA: ", sum(!is.na(plots1$x_clean)),
        " | y non-NA: ", sum(!is.na(plots1$y_clean)))
message("2012-2015 x range: ", round(min(plots1$x_clean, na.rm=TRUE),2),
        " - ", round(max(plots1$x_clean, na.rm=TRUE),2))
message("2012-2015 y range: ", round(min(plots1$y_clean, na.rm=TRUE),2),
        " - ", round(max(plots1$y_clean, na.rm=TRUE),2))

ll1 <- assign_lonlat(plots1$x_clean, plots1$y_clean)
plots1$longitude <- ll1$longitude
plots1$latitude  <- ll1$latitude

message("2012-2015 assigned lon non-NA: ", sum(!is.na(plots1$longitude)),
        " | lat non-NA: ", sum(!is.na(plots1$latitude)))

# Parse checkerboard
for (cn in cb_defs$name) plots1[[cn]] <- parse_cb(as.numeric(plots1[[cn]]))

valid1 <- valid_coords(plots1$longitude, plots1$latitude)
plots1_clean <- plots1 %>% filter(valid1)

message("2012-2015: ", nrow(plots1_clean), " valid plots")


# ==============================================================================
# Standardise column names and combine
# ==============================================================================

# Map 2016 column names -> standard names
cb_map <- data.frame(
  std = cb_defs$name,
  dir = rep(c("North", "East", "South", "West"), each = 4),
  dist = rep(c("5", "5", "10", "10"), 4),
  ht = rep(c("0", "1.25", "0", "1.25"), 4),
  stringsAsFactors = FALSE
)

cb_rename <- setNames(character(0), character(0))
for (i in 1:nrow(cb_map)) {
  pat <- paste0(cb_map$dir[i], ".*", cb_map$dist[i], " m.*", cb_map$ht[i], " m")
  matched <- cb_cols_2[grepl(pat, cb_cols_2, ignore.case = TRUE)]
  if (length(matched) >= 1) cb_rename[matched[1]] <- cb_map$std[i]
}

#  rename() needs new_name = old_name
cb_rename_fwd <- setNames(names(cb_rename), cb_rename)

veg2_std <- veg2_cb %>%
  dplyr::select(longitude, latitude, all_of(names(cb_rename))) %>%
  rename(!!!cb_rename_fwd) %>%
  mutate(dataset = "2016")

veg1_std <- plots1_clean %>%
  dplyr::select(longitude, latitude, all_of(cb_defs$name)) %>%
  mutate(dataset = "2012-2015")

veg_all <- bind_rows(veg1_std, veg2_std)
message("Combined: ", nrow(veg_all), " plots")


# ==============================================================================
# Convert to UTM and offset to measurement locations
# ==============================================================================

veg_sf_cb <- st_as_sf(veg_all, coords = c("longitude", "latitude"), crs = 4326)
veg_utm_cb <- st_transform(veg_sf_cb, crs = 32735)
coords_cb <- st_coordinates(veg_utm_cb)
veg_all$utm_x <- coords_cb[, 1]
veg_all$utm_y <- coords_cb[, 2]

# Build long-format: one row per measurement per height
build_pts <- function(col_name, height_label) {
  lapply(1:nrow(offset_defs), function(i) {
    col <- offset_defs[[col_name]][i]
    veg_all %>%
      dplyr::select(utm_x, utm_y, dataset, vis = !!sym(col)) %>%
      filter(!is.na(vis)) %>%
      mutate(
        meas_x = utm_x + offset_defs$dx[i],
        meas_y = utm_y + offset_defs$dy[i],
        direction = offset_defs$direction[i],
        distance_m = offset_defs$distance_m[i],
        height = height_label
      )
  }) %>% bind_rows()
}

ground_pts <- build_pts("cb_ground", "ground")
chest_pts  <- build_pts("cb_chest", "chest")

all_pts <- bind_rows(ground_pts, chest_pts) %>%
  rename(plot_utm_x = utm_x, plot_utm_y = utm_y) %>%
  mutate(vis_prop = vis / 225) %>%
  dplyr::select(plot_utm_x, plot_utm_y, meas_x, meas_y,
         direction, distance_m, height, vis, vis_prop, dataset)

message("Ground points: ", nrow(ground_pts))
message("Chest points: ", nrow(chest_pts))
message("Total: ", nrow(all_pts))

write.csv(all_pts, "checkerboard_visibility_points.csv", row.names = FALSE)

