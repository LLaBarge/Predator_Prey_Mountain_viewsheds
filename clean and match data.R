
library(ggplot2)
remove.packages("rlang")
install.packages("rlang")

# Ensure New_Timestamp is POSIXct
matching_data$New_Timestamp <- as.POSIXct(matching_data$New_Timestamp, 
                                          format = "%Y-%m-%d %H:%M:%S")
matching_data$Date <- as.Date(matching_data$New_Timestamp)

# Separate by species
baboon_data <- matching_data %>% 
  filter(species == "Baboon") %>%
  arrange(New_Timestamp)

vervet_data <- matching_data %>% 
  filter(species == "Vervet") %>%
  arrange(New_Timestamp)

cat(sprintf("Baboon points: %d\n", nrow(baboon_data)))
cat(sprintf("Vervet points: %d\n\n", nrow(vervet_data)))

# species trajectories plots

p1 <- ggplot(matching_data, aes(x = longitude, y = latitude, color = species)) +
  geom_path(aes(group = interaction(species, Date)), alpha = 0.3) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Baboon" = "#D55E00", "Vervet" = "#009E73")) +
  theme_minimal() +
  labs(title = "Movement Trajectories by Species",
       x = "Longitude", y = "Latitude", color = "Species") +
  coord_fixed(ratio = 1)
p1
ggsave("movement_analysis/trajectories_all.png", p1, width = 10, height = 8, dpi = 300)

# Baboon trajectories only
p2 <- ggplot(baboon_data, aes(x = longitude, y = latitude)) +
  geom_path(aes(group = Date), alpha = 0.5, color = "#D55E00") +
  geom_point(size = 1, alpha = 0.6, color = "#D55E00") +
  theme_minimal() +
  labs(title = "Baboon Movement Trajectories",
       x = "Longitude", y = "Latitude") +
  coord_fixed(ratio = 1)

ggsave("movement_analysis/trajectories_baboon.png", p2, width = 10, height = 8, dpi = 300)

# Vervet trajectories only
p3 <- ggplot(vervet_data, aes(x = longitude, y = latitude)) +
  geom_path(aes(group = Date), alpha = 0.5, color = "#009E73") +
  geom_point(size = 1, alpha = 0.6, color = "#009E73") +
  theme_minimal() +
  labs(title = "Vervet Movement Trajectories",
       x = "Longitude", y = "Latitude") +
  coord_fixed(ratio = 1)
p3
ggsave("movement_analysis/trajectories_vervet.png", p3, width = 10, height = 8, dpi = 300)

cat("Trajectory plots saved\n\n")

# filter unrealistic movement

# Function to calculate speeds and filter (5km/hr)
filter_by_speed <- function(data, species_name, max_speed_kmh = 5) {
  
  data <- data %>%
    mutate(
      next_lon = lead(longitude),
      next_lat = lead(latitude),
      next_time = lead(New_Timestamp),
      
      distance_m = geosphere::distHaversine(
        cbind(longitude, latitude),
        cbind(next_lon, next_lat)
      ),
      
      time_diff_hours = as.numeric(difftime(next_time, New_Timestamp, units = "hours")),
      
      speed_kmh = (distance_m / 1000) / time_diff_hours,
      
      unrealistic = case_when(
        is.na(speed_kmh) ~ FALSE,
        is.na(time_diff_hours) ~ FALSE,
        time_diff_hours >= 2 ~ FALSE,
        speed_kmh > max_speed_kmh ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  n_unrealistic <- sum(data$unrealistic, na.rm = TRUE)
  
  cat(sprintf("  Total points: %d\n", nrow(data)))
  cat(sprintf("  Unrealistic: %d (%.2f%%)\n", 
              n_unrealistic, 100 * n_unrealistic / nrow(data)))
  
  # Plot speed distribution
  if (sum(!is.na(data$speed_kmh)) > 0) {
    p <- ggplot(data %>% filter(!is.na(speed_kmh), speed_kmh < 100), 
                aes(x = speed_kmh)) +
      geom_histogram(bins = 50, fill = ifelse(species_name == "Baboon", "#D55E00", "#009E73"), 
                     alpha = 0.7) +
      geom_vline(xintercept = max_speed_kmh, color = "red", 
                 linetype = "dashed", size = 1) +
      theme_minimal() +
      labs(title = paste(species_name, "Movement Speeds"),
           subtitle = paste("Red line: threshold =", max_speed_kmh, "km/h"),
           x = "Speed (km/h)", y = "Count")
    
    ggsave(paste0("movement_analysis/speed_distribution_", tolower(species_name), ".png"), 
           p, width = 8, height = 6, dpi = 300)
  }
  
  # Remove unrealistic points
  data_clean <- data %>%
    filter(unrealistic == FALSE) %>%
    dplyr::select(-next_lon, -next_lat, -next_time, -distance_m, 
                  -time_diff_hours, -speed_kmh, -unrealistic)
  
  cat(sprintf("  Remaining: %d\n\n", nrow(data_clean)))
  
  return(data_clean)
}

# Filter each species
baboon_clean <- filter_by_speed(baboon_data, "Baboon", max_speed_kmh = 5)
vervet_clean <- filter_by_speed(vervet_data, "Vervet", max_speed_kmh = 5)
str(vervet_clean)
# Combine cleaned data
matching_data_clean <- bind_rows(baboon_clean, vervet_clean)

# save
saveRDS(vervet_clean, "vervet_clean.rds")
saveRDS(baboon_clean, "baboon_clean.rds")
