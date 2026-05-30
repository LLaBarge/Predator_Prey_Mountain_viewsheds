
library(readxl)
library(dplyr)
library(mgcv)
library(terra)
library(sf)
library(elevatr)
library(amt)

# baboon data goes further in each direction than vervets so use this for the bounding box
baboon_sf <- st_as_sf(baboon_clean, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 32735)

bbox_utm <- st_bbox(baboon_sf)
str(bbox_utm)

# Create bounding box as sf object for DEM download
bbox_sf <- st_bbox(bbox_utm, crs = 32735) %>%
  st_as_sfc() %>%
  st_sf()

# Transform to WGS84 for elevatr (it requires lat/long)
bbox_wgs84 <- st_transform(bbox_sf, crs = 4326)

dem_download <- get_elev_raster(
  locations = bbox_wgs84,
  z = 12,  # Adjust this if needed: down=coarse, increase=fine
  src = "aws"  # Uses AWS Terrain Tiles - free
)

# Convert to terra SpatRaster and reproject to UTM
dem_wgs84 <- rast(dem_download)

dem <- project(dem_wgs84, "EPSG:32735", method = "bilinear")

# Crop to exact study area
dem <- crop(dem, ext(bbox_utm))

# Clean up to save RAM
rm(dem_download, dem_wgs84)
gc()

# Save DEM for future use
writeRaster(dem, "dem.tif", overwrite = TRUE)


plot(dem, main = "Elevation (m)")


# Convert to sf objects 
vervet_sf <- st_as_sf(vervet_clean,  coords = c("longitude", "latitude"), crs = 4326)
baboon_sf <- st_as_sf(baboon_clean,
                      coords = c("longitude", "latitude"), 
                      crs = 4326)

# Transform to UTM 35S
vervet_sf <- st_transform(vervet_sf, crs = 32735)
baboon_sf <- st_transform(baboon_sf, crs = 32735)

# Convert to terra vectors for plotting
vervet_pts <- vect(vervet_sf)
baboon_pts <- vect(baboon_sf)

# create tracks
vervet_track <- vervet_clean %>%
  make_track(longitude, latitude, New_Timestamp,
             crs = 4326, all_cols = TRUE) %>%
  transform_coords(32735)

baboon_track <- baboon_clean %>%
  make_track(longitude, latitude, New_Timestamp,
             crs = 4326, all_cols = TRUE) %>%
  transform_coords(32735)  


# Effects of topography on detection
#Calculate elevation difference between two points
calc_elevation_difference <- function(x1, y1, x2, y2, dem) {
  elev1 <- extract(dem, cbind(x1, y1))[,1]
  elev2 <- extract(dem, cbind(x2, y2))[,1]
  return(abs(elev1 - elev2))
}

#find maximum terrain obstacle between two points
calc_max_terrain_obstacle <- function(x1, y1, x2, y2, dem, n_samples = 20) {
  elev1 <- extract(dem, cbind(x1, y1))[,1]
  elev2 <- extract(dem, cbind(x2, y2))[,1]
  
  if (is.na(elev1) || is.na(elev2)) return(NA)
  
  sample_x <- seq(x1, x2, length.out = n_samples)
  sample_y <- seq(y1, y2, length.out = n_samples)
  line_elevs <- extract(dem, cbind(sample_x, sample_y))[,1]
  
  if (any(is.na(line_elevs))) return(NA)
  
  max_elev <- max(line_elevs)
  baseline <- min(elev1, elev2)
  obstacle_height <- max_elev - baseline
  
  return(obstacle_height)
}

#Calculate all 3 topographic metrics between vervet and baboon
calculate_topo_metrics <- function(ver_x, ver_y, bab_x, bab_y, dem) {
  
  # Horizontal distance
  dist_euclidean <- sqrt((ver_x - bab_x)^2 + (ver_y - bab_y)^2)
  
  # Elevation difference
  elev_diff <- calc_elevation_difference(ver_x, ver_y, bab_x, bab_y, dem)
  
  # Maximum terrain obstacle
  max_terrain_obstacle <- calc_max_terrain_obstacle(ver_x, ver_y, bab_x, bab_y, dem)
  
  return(data.frame(
    dist_euclidean = dist_euclidean,
    elev_diff = elev_diff,
    max_terrain_obstacle = max_terrain_obstacle
  ))
}




# Calculate for a single pair of locations
metrics <- calculate_topo_metrics(ver_x = 750000, ver_y = 7449000, 
                                  bab_x = 750500, bab_y = 7449500, 
                                  dem = dem)

# Calculate for multiple pairs
results <- data.frame()
for (i in 1:nrow(vervet_track)) {
  metrics <- calculate_topo_metrics(
    vervet_track$x_[i], vervet_track$y_[i],
    baboon_track$x_[i], baboon_track$y_[i],
    dem
  )
  results <- rbind(results, metrics)
}


vervet_track$dist_euclidean <- results$dist_euclidean
vervet_track$elev_diff <- results$elev_diff
vervet_track$max_terrain_obstacle <- results$max_terrain_obstacle

# check the data is correct  
print(vervet_track)


