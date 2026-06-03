# Observer detection validation
# Does viewshed predict whether human observers detected baboons?
# Follows from the missed encounters analysis.


library(dplyr)
library(sf)
library(raster)
library(brms)
library(ggplot2)

dem <- raster("dem.tif")
canopy_height_raster <- crop(raster("meta_tree_height.tif"), dem)

encounters <- read.csv("encounter_behavioral_analysis.csv", stringsAsFactors = FALSE)

# Load GPS data to reconstruct positions at closest approach
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
  vervet_clean$timestamp <- as.POSIXct(vervet_clean$New_Timestamp, tz = "UTC")
  baboon_clean$timestamp <- as.POSIXct(baboon_clean$New_Timestamp, tz = "UTC")
} else {
  vervet_clean$timestamp <- as.POSIXct(vervet_clean$timestamp, tz = "UTC")
  baboon_clean$timestamp <- as.POSIXct(baboon_clean$timestamp, tz = "UTC")
}

# Find vervet and baboon UTM positions at the moment of closest approach
encounters$vervet_x <- NA
encounters$vervet_y <- NA
encounters$baboon_x <- NA
encounters$baboon_y <- NA

for (i in 1:nrow(encounters)) {
  t_start <- as.POSIXct(encounters$start_time[i], tz = "UTC")
  t_end   <- as.POSIXct(encounters$end_time[i], tz = "UTC")
  
  # For zero-duration encounters, expand window by +-20 min (a scan window)
  if (t_start == t_end) {
    t_start <- t_start - 1200
    t_end   <- t_end + 1200
  }
  
  ver_fixes <- vervet_clean %>%
    filter(timestamp >= t_start & timestamp <= t_end)
  bab_fixes <- baboon_clean %>%
    filter(timestamp >= t_start & timestamp <= t_end)
  
  if (nrow(ver_fixes) == 0 || nrow(bab_fixes) == 0) next
  
  min_dist <- Inf
  for (v in 1:nrow(ver_fixes)) {
    time_diffs <- abs(difftime(bab_fixes$timestamp, ver_fixes$timestamp[v], units = "mins"))
    close_bab <- which(time_diffs <= 20)
    for (b in close_bab) {
      d <- sqrt((ver_fixes$x_utm[v] - bab_fixes$x_utm[b])^2 +
                  (ver_fixes$y_utm[v] - bab_fixes$y_utm[b])^2)
      if (d < min_dist) {
        min_dist <- d
        encounters$vervet_x[i] <- ver_fixes$x_utm[v]
        encounters$vervet_y[i] <- ver_fixes$y_utm[v]
        encounters$baboon_x[i] <- bab_fixes$x_utm[b]
        encounters$baboon_y[i] <- bab_fixes$y_utm[b]
      }
    }
  }
}

print(paste("Encounters with reconstructed positions:",
            sum(!is.na(encounters$vervet_x)), "/", nrow(encounters)))

# Sector viewshed: DEM + tree height only (no canopy cover)
# Observer at 1.5m (human on ground), target at 1.5m (baboon on ground)
# 45 degree sector, 9 rays

viewshed_sector <- function(obs_x, obs_y, obs_h, tgt_x, tgt_y, tgt_h,
                            dem, canopy_height_raster,
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
    ch[is.na(ch)] <- 0
    
    obstacle <- ter + ch
    frac <- seq(0, 1, length.out = n_pts)
    los <- obs_elev + frac * ((rg + tgt_h) - obs_elev)
    cl <- los[2:(n_pts-1)] - obstacle[2:(n_pts-1)]
    
    total <- total + 1
    if (min(cl, na.rm = TRUE) >= 0) clear <- clear + 1
  }
  if (total == 0) return(NA)
  clear / total
}

# Calculate viewshed at each encounter
encounters$observer_viewshed <- NA

for (i in 1:nrow(encounters)) {
  encounters$observer_viewshed[i] <- viewshed_sector(
    obs_x = encounters$vervet_x[i],
    obs_y = encounters$vervet_y[i],
    obs_h = 1.7,
    tgt_x = encounters$baboon_x[i],
    tgt_y = encounters$baboon_y[i],
    tgt_h = 1.7,
    dem, canopy_height_raster
  )
}

encounters$detected <- ifelse(encounters$detection_type == "explicit", 1, 0)

# Identify complete cases
encounters$complete <- !is.na(encounters$viewshed_sc) & !is.na(encounters$dist_sc)
print(paste("Complete:", sum(encounters$complete), "/ Missing:", sum(!encounters$complete)))

# Fit ALL models on the complete subset for fair comparison
enc_complete <- encounters %>% filter(complete)
enc_complete$dist_sc     <- scale(enc_complete$min_distance_m)[, 1]
enc_complete$viewshed_sc <- scale(enc_complete$observer_viewshed)[, 1]


# Priors (standardised scale)
# Intercept: logit(0.29) ~ -0.9, evaluated at mean distance and mean viewshed
#   Normal(-1, 1.5) centres near the observed detection rate
# Distance: detection decreases with distance. A 1-SD increase in distance
#   should moderately reduce log-odds of detection.
#   Normal(-1, 1) — weakly informative, expects negative effect
# Viewshed: clearer view should increase detection. A 1-SD increase in
#   visibility should moderately increase log-odds.
#   Normal(1, 1) — weakly informative, expects positive effect
# Interaction: no strong prior expectation on how viewshed modifies
#   the distance effect for human observers.
#   Normal(0, 1) — weakly informative, centred at zero

priors_additive <- c(
  prior(normal(-1, 1.5), class = "Intercept"),
  prior(normal(-1, 1),   class = "b", coef = "dist_sc"),
  prior(normal(1, 1),    class = "b", coef = "viewshed_sc")
)

priors_interaction <- c(
  prior(normal(-1, 1.5), class = "Intercept"),
  prior(normal(-1, 1),   class = "b", coef = "dist_sc"),
  prior(normal(1, 1),    class = "b", coef = "viewshed_sc"),
  prior(normal(0, 1),    class = "b", coef = "dist_sc:viewshed_sc")
)

# Prior sensitivity analysis

priors_weakly <- c(
  prior(normal(-1, 3),  class = "Intercept"),
  prior(normal(-1, 2),  class = "b", coef = "dist_sc"),
  prior(normal(1, 2),   class = "b", coef = "viewshed_sc")
)

priors_diffuse <- c(
  prior(normal(0, 5),   class = "Intercept"),
  prior(normal(0, 5),   class = "b", coef = "dist_sc"),
  prior(normal(0, 5),   class = "b", coef = "viewshed_sc")
)

# Prior predictive check: sample from priors only (no likelihood)
m_prior_only <- brm(
  detected ~ dist_sc + viewshed_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = priors_additive,
  sample_prior = "only",
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

print("Prior predictive distributions (no data)")
print(summary(m_prior_only))

# Fit additive model under each prior specification
m_informative <- brm(
  detected ~ dist_sc + viewshed_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = priors_additive,
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

m_weakly <- brm(
  detected ~ dist_sc + viewshed_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = priors_weakly,
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

m_diffuse <- brm(
  detected ~ dist_sc + viewshed_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = priors_diffuse,
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

# Compare posteriors across prior specifications
fe_inf <- fixef(m_informative)
fe_wk  <- fixef(m_weakly)
fe_dif <- fixef(m_diffuse)

# Get column names for upper CI 
ci_upper_col <- grep("97", colnames(fe_inf), value = TRUE)[1]

sensitivity_comparison <- data.frame(
  prior_set = rep(c("Informative", "Weakly informative", "Diffuse"), each = 3),
  parameter = rep(c("Intercept", "dist_sc", "viewshed_sc"), 3),
  estimate  = c(fe_inf[, "Estimate"], fe_wk[, "Estimate"], fe_dif[, "Estimate"]),
  lower     = c(fe_inf[, "Q2.5"],     fe_wk[, "Q2.5"],     fe_dif[, "Q2.5"]),
  upper     = c(fe_inf[, ci_upper_col], fe_wk[, ci_upper_col], fe_dif[, ci_upper_col])
)


print(sensitivity_comparison)

# Plot sensitivity
sensitivity_comparison$prior_set <- factor(sensitivity_comparison$prior_set,
                                           levels = c("Diffuse", "Weakly informative", "Informative"))

p_sens <- ggplot(sensitivity_comparison,
                 aes(x = estimate, y = prior_set, colour = prior_set)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ parameter, scales = "free_x") +
  scale_colour_manual(values = c("Diffuse" = "grey50",
                                 "Weakly informative" = "#2E86AB",
                                 "Informative" = "#D32F2F")) +
  labs(x = "Posterior estimate (95% CrI)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave("prior_sensitivity.png", p_sens, width = 10, height = 4, dpi = 300)
p_sens
#look similar so use informative

# Model fitting with chosen priors
# Model 1: distance only
m1 <- brm(
  detected ~ dist_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = c(
    prior(normal(-1, 1.5), class = "Intercept"),
    prior(normal(-1, 1),   class = "b", coef = "dist_sc")
  ),
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

# Model 2: viewshed only
m2 <- brm(
  detected ~ viewshed_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = c(
    prior(normal(-1, 1.5), class = "Intercept"),
    prior(normal(1, 1),    class = "b", coef = "viewshed_sc")
  ),
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

# Model 3: distance + viewshed (additive) — already fitted during sensitivity analysis
m3 <- m_informative

# Model 4: distance * viewshed (interaction)
m4 <- brm(
  detected ~ dist_sc * viewshed_sc,
  family = bernoulli(),
  data = enc_complete,
  prior = priors_interaction,
  chains = 4, iter = 4000, warmup = 1000,
  seed = 42, silent = 2
)

# Results
print("Model 1: distance only")
print(summary(m1))

print("Model 2: viewshed only")
print(summary(m2))

print("Model 3: distance + viewshed")
print(summary(m3))

print("Model 4: distance * viewshed interaction")
print(summary(m4))

# Model comparison via LOO
m1 <- add_criterion(m1, "loo", moment_match=TRUE)
m2 <- add_criterion(m2, "loo", moment_match=TRUE)
m3 <- add_criterion(m3, "loo", moment_match=TRUE)
m4 <- add_criterion(m4, "loo", moment_match=TRUE)


loo_compare(m1, m2, m3, m4)

# Descriptive summary
print("Detection by type — viewshed and distance")
encounters %>%
  group_by(detection_type) %>%
  summarize(
    n = n(),
    mean_viewshed = mean(observer_viewshed, na.rm = TRUE),
    sd_viewshed = sd(observer_viewshed, na.rm = TRUE),
    mean_distance = mean(min_distance_m, na.rm = TRUE),
    sd_distance = sd(min_distance_m, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

# Conditional effects plot
p_ce <- plot(conditional_effects(m4), plot = FALSE)
p_ce

# Conditional effects plot for m4 (interaction model)
# Predicted detection probability across distance at distinct viewshed levels

library(tidybayes)

# Create prediction grid
dist_range <- seq(min(enc_complete$dist_sc), max(enc_complete$dist_sc), length.out = 200)
viewshed_bins <- c(-1, 0, 1)

pred_grid <- expand.grid(
  dist_sc     = dist_range,
  viewshed_sc = viewshed_bins
)

pred_grid$viewshed_label <- factor(
  pred_grid$viewshed_sc,
  levels = viewshed_bins,
  labels = c("Blocked (-1 SD)", "Median viewshed", "Clear (+1 SD)")
)

fits <- fitted(m4, newdata = pred_grid, summary = TRUE)
pred_grid$estimate <- fits[, "Estimate"]
pred_grid$lower    <- fits[, "Q2.5"]
pred_grid$upper    <- fits[, "Q97.5"]

p_cond <- ggplot(pred_grid, aes(x = dist_sc, y = estimate, colour = viewshed_label,
                                fill = viewshed_label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = c("#D32F2F", "grey50", "#1565C0"),
    name = "Sector viewshed"
  ) +
  scale_fill_manual(
    values = c("#D32F2F", "grey50", "#1565C0"),
    name = "Sector viewshed"
  ) +
  labs(
    x = "Distance (scaled)\n\u2190 Closer                     Farther \u2192",
    y = "P(observer detected)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave("observer_detection_conditional.png", p_cond, width = 9, height = 5.5, dpi = 300)
p_cond
# Save results
write.csv(encounters[, c("encounter_id", "date", "detection_type", "detected",
                         "min_distance_m", "observer_viewshed",
                         "dist_sc", "viewshed_sc")],
          "observer_detection_validation.csv", row.names = FALSE)