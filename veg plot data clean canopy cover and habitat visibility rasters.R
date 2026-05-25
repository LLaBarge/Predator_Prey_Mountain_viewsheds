

library(readxl)
library(dplyr)
library(sf)
library(terra)


veg_plots <- read_excel("vegplot_phenology.xlsx", sheet = 1)
veg_plots[is.na(veg_plots)] <- 0

# found some issues with missing "-" manually - correct this here:
# Find latitude column
lat_col <- names(veg_plots)[grepl("^Lat", names(veg_plots), ignore.case = TRUE)][1]
lon_col <- names(veg_plots)[grepl("^Long", names(veg_plots), ignore.case = TRUE)][1]

if (is.null(lat_col) || is.null(lon_col)) {
  stop("Could not find Lat/Long columns in vegetation data")
}



# Check for positive latitudes (should all be negative)
positive_lats <- veg_plots[[lat_col]] > 0
n_positive <- sum(positive_lats, na.rm = TRUE)

if (n_positive > 0) {
  cat("  ⚠ WARNING: Found", n_positive, "plots with POSITIVE latitude\n")
  cat("  These should be negative (Southern Hemisphere)\n")
  cat("  FIXING: Converting to negative...\n\n")
  
  # Show examples before fix
  cat("  Examples BEFORE fix:\n")
  bad_idx <- which(positive_lats)[1:min(5, n_positive)]
  for (i in bad_idx) {
    cat("    Plot", i, ": Lat =", veg_plots[[lat_col]][i], 
        "Long =", veg_plots[[lon_col]][i], "\n")
  }
  
  # Make all latitudes negative
  veg_plots[[lat_col]][positive_lats] <- -abs(veg_plots[[lat_col]][positive_lats])
  
  # Show examples
  for (i in bad_idx) {
    cat("    Plot", i, ": Lat =", veg_plots[[lat_col]][i], 
        "Long =", veg_plots[[lon_col]][i], "\n")
  }
  
} else {
  
}

# Verify
still_positive <- sum(veg_plots[[lat_col]] > 0, na.rm = TRUE)

still_positive

# vegetation data
species_cols <- names(veg_plots)[grepl("^[A-Z]\\.", names(veg_plots))]
veg_plots$total_trees <- rowSums(veg_plots[, species_cols])
veg_plots$species_richness <- rowSums(veg_plots[, species_cols] > 0)

# Convert to UTM 
veg_sf <- st_as_sf(veg_plots, coords = c(lon_col, lat_col), crs = 4326)
veg_utm <- st_transform(veg_sf, crs = 32735)
coords <- st_coordinates(veg_utm)
veg_plots$UTM_X <- coords[, 1]
veg_plots$UTM_Y <- coords[, 2]
plot(coords)
# Extract variables
veg_plots$HV <- veg_plots$`Horizontal.vis(m)`
veg_plots$CC <- veg_plots$`Canopy.Cover(%)`


if (!exists("vervet_clean")) {
}

# Get vervet extent
if ("x_utm" %in% names(vervet_clean)) {
  vervet_x <- vervet_clean$x_utm
  vervet_y <- vervet_clean$y_utm
} else {
  # Convert to UTM if needed
  ver_lon_col <- names(vervet_clean)[grepl("longitude|lon", names(vervet_clean), ignore.case = TRUE)][1]
  ver_lat_col <- names(vervet_clean)[grepl("latitude|lat", names(vervet_clean), ignore.case = TRUE)][1]
  
  vervet_sf <- st_as_sf(vervet_clean, coords = c(ver_lon_col, ver_lat_col), crs = 4326)
  vervet_utm <- st_transform(vervet_sf, crs = 32735)
  coords <- st_coordinates(vervet_utm)
  vervet_x <- coords[, 1]
  vervet_y <- coords[, 2]
}

# Get bounding box
vervet_bbox <- c(
  xmin = min(vervet_x, na.rm = TRUE),
  xmax = max(vervet_x, na.rm = TRUE),
  ymin = min(vervet_y, na.rm = TRUE),
  ymax = max(vervet_y, na.rm = TRUE)
)


# create template
RESOLUTION <- 50  # metres

raster_template <- rast(
  xmin = vervet_bbox["xmin"],
  xmax = vervet_bbox["xmax"],
  ymin = vervet_bbox["ymin"],
  ymax = vervet_bbox["ymax"],
  resolution = RESOLUTION,
  crs = "EPSG:32736"
)

cat("Raster template (VERVET EXTENT ONLY):\n")
cat("  Resolution:", RESOLUTION, "m\n")
cat("  Size:", nrow(raster_template), "x", ncol(raster_template), "\n")
cat("  Total cells:", ncell(raster_template), "\n")
cat("  Area:", round((vervet_bbox["xmax"] - vervet_bbox["xmin"]) * 
                       (vervet_bbox["ymax"] - vervet_bbox["ymin"]) / 1e6, 2), "km²\n\n")

# filter veg plots to vervet extent
buffer_m <- 500

veg_filtered <- veg_plots %>%
  filter(
    UTM_X >= (vervet_bbox["xmin"] - buffer_m),
    UTM_X <= (vervet_bbox["xmax"] + buffer_m),
    UTM_Y >= (vervet_bbox["ymin"] - buffer_m),
    UTM_Y <= (vervet_bbox["ymax"] + buffer_m)
  )


# Convert to terra points
points_vect <- vect(veg_filtered, geom = c("UTM_X", "UTM_Y"), crs = "EPSG:32736")

# interpolate rast
idw_fast <- function(points, template, field, idp = 2, nmax = 12) {
  
  cat("  Interpolating", field, "...\n")
  
  # Get point data
  coords_mat <- crds(points)
  values_vec <- values(points)[[field]]
  n_points <- nrow(coords_mat)
  
  # Get all cell coordinates
  cell_coords <- xyFromCell(template, 1:ncell(template))
  n_cells <- nrow(cell_coords)
  
  cat("    Points:", n_points, "\n")
  cat("    Cells:", n_cells, "\n")
  
  # Initialize output
  output_values <- rep(NA_real_, n_cells)
  
  # Progress tracking
  progress_interval <- max(1, round(n_cells / 20))
  
  # Process each cell
  for (i in 1:n_cells) {
    
    if (i %% progress_interval == 0) {
      cat("\r    Progress:", round(100 * i / n_cells, 1), "%")
    }
    
    # Calculate distances
    dx <- coords_mat[, 1] - cell_coords[i, 1]
    dy <- coords_mat[, 2] - cell_coords[i, 2]
    distances <- sqrt(dx^2 + dy^2)
    
    # Handle exact match
    if (any(distances == 0)) {
      output_values[i] <- values_vec[which.min(distances)]
      next
    }
    
    # Get nearest neighbors
    if (n_points > nmax) {
      nearest_idx <- order(distances)[1:nmax]
    } else {
      nearest_idx <- 1:n_points
    }
    
    # IDW
    d_nearest <- distances[nearest_idx]
    v_nearest <- values_vec[nearest_idx]
    weights <- 1 / (d_nearest^idp)
    
    output_values[i] <- sum(weights * v_nearest) / sum(weights)
  }
  
  cat("\r    Progress: 100%  \n")
  
  # Create raster
  result <- rast(template)
  values(result) <- output_values
  
  return(result)
}


# Horizontal Visibility
hv_raster <- idw_fast(points_vect, raster_template, "HV", idp = 2, nmax = 12)
names(hv_raster) <- "Horizontal_Visibility"

# Canopy Cover
cc_raster <- idw_fast(points_vect, raster_template, "CC", idp = 2, nmax = 12)
names(cc_raster) <- "Canopy_Cover"

plot(hv_raster)
plot(cc_raster)
writeRaster(hv_raster, "horizontal_visibility.tif", overwrite = TRUE)
writeRaster(cc_raster, "canopy_cover.tif", overwrite = TRUE)

# Load rasters
hv_raster <- rast("horizontal_visibility.tif")
cc_raster <- rast("canopy_cover.tif")

