# explore encounters missed by human observers

# check how many potential "missed" encounters there were

library(dplyr)
library(ggplot2)

# encounter data
if (!file.exists("encounter_behavioral_analysis.csv")) {
  stop("Run encounter_behavioral_analysis.R first to generate encounter data")
}

encounters <- read.csv("encounter_behavioral_analysis.csv", stringsAsFactors = FALSE)

cat("Total encounters:", nrow(encounters), "\n")
cat("  Explicit (observed):", sum(encounters$detection_type == "explicit"), "\n")
cat("  Implicit (not observed):", sum(encounters$detection_type == "implicit"), "\n\n")

#analyse across distance thresholds
distance_thresholds <- c(50, 100, 150, 200, 300, 500)

results_table <- data.frame()

for (threshold in distance_thresholds) {
  
  # Filter to encounters within this distance
  within_threshold <- encounters %>% filter(min_distance_m <= threshold)
  
  n_total <- nrow(within_threshold)
  n_explicit <- sum(within_threshold$detection_type == "explicit")
  n_implicit <- sum(within_threshold$detection_type == "implicit")
  pct_missed <- if (n_total > 0) 100 * n_implicit / n_total else 0
  
  results_table <- rbind(results_table, data.frame(
    distance_threshold_m = threshold,
    total_encounters = n_total,
    observed = n_explicit,
    missed = n_implicit,
    pct_missed = pct_missed
  ))
  
  cat("Within", threshold, "m:\n")
  cat("  Total encounters:", n_total, "\n")
  cat("  Observed (explicit):", n_explicit, "\n")
  cat("  Missed (implicit):", n_implicit, "\n")
  cat("  % Missed:", sprintf("%.1f%%", pct_missed), "\n\n")
}

# chose 200 as almost all encounters vervets respond to (by moving away/running)
within_200m <- encounters %>% filter(min_distance_m <= 200)

n_total_200 <- nrow(within_200m)
n_observed_200 <- sum(within_200m$detection_type == "explicit")
n_missed_200 <- sum(within_200m$detection_type == "implicit")
pct_missed_200 <- if (n_total_200 > 0) 100 * n_missed_200 / n_total_200 else 0

cat("Encounters within 200m:\n")
cat("  Total:", n_total_200, "\n")
cat("  Observed:", n_observed_200, "\n")
cat("  Missed:", n_missed_200, "\n")
cat("  % Missed:", sprintf("%.1f%%", pct_missed_200), "\n\n")

if (n_missed_200 > 0) {
  cat("Details of missed encounters (<200m):\n")
  missed_200 <- within_200m %>% 
    filter(detection_type == "implicit") %>%
    arrange(min_distance_m)
  
  cat("  Closest missed encounter:", round(min(missed_200$min_distance_m)), "m\n")
  cat("  Median distance of missed:", round(median(missed_200$min_distance_m)), "m\n")
  cat("  Mean distance of missed:", round(mean(missed_200$min_distance_m)), "m\n\n")
  
  # Show first few
  cat("First 10 missed encounters within 200m:\n")
  print(missed_200[1:min(10, nrow(missed_200)), 
                   c("encounter_id", "date", "min_distance_m", "duration_min")])
}

# Detection probability by distance
# Create distance bins
encounters$distance_bin <- cut(encounters$min_distance_m, 
                               breaks = c(0, 50, 100, 150, 200, 300, 500),
                               labels = c("<50", "50-100", "100-150", "150-200", 
                                          "200-300", "300-500"))

detection_by_distance <- encounters %>%
  group_by(distance_bin) %>%
  summarize(
    n_encounters = n(),
    n_observed = sum(detection_type == "explicit"),
    n_missed = sum(detection_type == "implicit"),
    pct_detected = 100 * n_observed / n(),
    .groups = "drop"
  )

print(detection_by_distance)

write.csv(results_table, "missed_encounters_by_distance.csv", row.names = FALSE)

if (n_missed_200 > 0) {
  write.csv(missed_200, "missed_encounters_within_150m.csv", row.names = FALSE)
}


# Detection rate by distance
Dect_dist <- ggplot(detection_by_distance, aes(x = distance_bin, y = pct_detected)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.0f%%", pct_detected)), 
            vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "Detection Rate by Distance",
       x = "Distance to Baboons (m)",
       y = "% of Encounters Observed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Dect_dist
ggsave("detection_rate_by_distance.png", Dect_dist, width = 8, height = 6, dpi = 300)

# Plot 2: Observed vs Missed by distance threshold
ob_vs_m <- ggplot(results_table, aes(x = distance_threshold_m)) +
  geom_line(aes(y = observed, color = "Observed"), size = 1) +
  geom_point(aes(y = observed, color = "Observed"), size = 3) +
  geom_line(aes(y = missed, color = "Missed"), size = 1) +
  geom_point(aes(y = missed, color = "Missed"), size = 3) +
  scale_color_manual(values = c("Observed" = "darkgreen", "Missed" = "red")) +
  theme_minimal() +
  labs(
    x = "Distance Threshold (m)",
    y = "Number of Encounters",
    color = "Detection Type") +
  geom_vline(xintercept = 200, linetype = "dashed", alpha = 0.5)

ggsave("observed_vs_missed_by_distance.png", ob_vs_m, width = 8, height = 6, dpi = 300)
ob_vs_m
# Percentage missed
missed <- ggplot(results_table, aes(x = distance_threshold_m, y = pct_missed)) +
  geom_line(color = "red", size = 1) +
  geom_point(size = 3, color = "red") +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Distance Threshold (m)",
    y = "% Missed by human observers") +
  scale_y_continuous(limits = c(0, 100))

ggsave("pct_missed_by_distance.png", missed, width = 8, height = 6, dpi = 300)
missed

# Overall detection rate w/in 500m
overall_pct_detected <- 100 * sum(encounters$detection_type == "explicit") / nrow(encounters)
cat("Overall detection rate:", sprintf("%.1f%%", overall_pct_detected), "\n\n")

# Mean distance of observed vs missed
mean_dist_observed <- mean(encounters$min_distance_m[encounters$detection_type == "explicit"])
max_dist_observed <- max(encounters$min_distance_m[encounters$detection_type == "explicit"])
mean_dist_missed <- mean(encounters$min_distance_m[encounters$detection_type == "implicit"])

mean_dist_observed
mean_dist_missed
max_dist_observed
