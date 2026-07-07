library(INLA)
library(dplyr)

# =============================================================================
# Load and prepare data

ver_ssf_data <- read.csv("ver_ssf_data_viewshed_sector.csv", stringsAsFactors = FALSE)

ssf_complete <- ver_ssf_data %>%
  filter(!is.na(dist_to_baboon),
         !is.na(sensoryshed),
         !is.na(hv_chest))

ssf_complete$dist_scaled        <- scale(ssf_complete$dist_to_baboon)[, 1]
ssf_complete$sensoryshed_scaled <- scale(ssf_complete$sensoryshed)[, 1]
ssf_complete$hv_chest_scaled    <- scale(ssf_complete$hv_chest)[, 1]
ssf_complete$dist_x_sensoryshed <- ssf_complete$dist_scaled * ssf_complete$sensoryshed_scaled

ssf_complete <- ssf_complete %>%
  group_by(step_id_) %>%
  mutate(
    n_steps_in_stratum = n(),
    weight = 1e6 / n_steps_in_stratum
  ) %>%
  ungroup()

ssf_complete$case_numeric <- as.integer(ssf_complete$case_)
ssf_complete$stratum_id   <- as.factor(ssf_complete$step_id_)

# Correlations
round(cor(ssf_complete[, c("dist_scaled", "sensoryshed_scaled",
                           "hv_chest_scaled")], use = "complete.obs"), 3)

# =============================================================================
# Prior sensitivity analysis

# INLA default: N(0, prec=0.001) i.e. sd~31.6 (effectively flat)
# Weakly informative: N(0, sd=1)   -- plausible ecological effect sizes
# Moderate:           N(0, sd=0.5) -- shrinks toward no effect
# Sceptical:          N(0, sd=0.3) -- strong regularisation

prior_specs <- list(
  "Default (sd~31.6)" = list(prec = 0.001),
  "N(0, sd=1.0)"      = list(prec = 1),
  "N(0, sd=0.5)"      = list(prec = 4),
  "N(0, sd=0.3)"      = list(prec = 11.1)
)

sensitivity_models <- list()

for (prior_name in names(prior_specs)) {
  sensitivity_models[[prior_name]] <- inla(
    case_numeric ~ -1 +
      dist_scaled +
      hv_chest_scaled +
      sensoryshed_scaled +
      dist_x_sensoryshed +
      f(stratum_id, model = "iid",
        hyper = list(prec = list(initial = -6, fixed = TRUE))),
    family = "poisson",
    data = ssf_complete,
    E = weight,
    control.fixed = list(mean = 0, prec = prior_specs[[prior_name]]$prec),
    control.compute = list(config = TRUE, mlik = TRUE),
    verbose = FALSE
  )
}

# coefficients under each prior
lapply(names(sensitivity_models), function(nm) {
  list(prior = nm, fixed = round(sensitivity_models[[nm]]$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3))
})

# dist_x_sensoryshed across priors
sapply(sensitivity_models, function(m) {
  dx <- m$summary.fixed["dist_x_sensoryshed", ]
  c(mean = dx$mean, lower = dx$`0.025quant`, upper = dx$`0.975quant`,
    excludes_zero = dx$`0.025quant` > 0 | dx$`0.975quant` < 0)
})

# Choose prior and fit final models
#
# N(0, sd=1): constrains coefficients to ecologically plausible ranges

chosen_prior <- list(prec = 1)  # N(0, sd=1)

# M1: Distance + local habitat (no sensoryshed)
m1 <- inla(
  case_numeric ~ -1 +
    dist_scaled +
    hv_chest_scaled +
    f(stratum_id, model = "iid",
      hyper = list(prec = list(initial = -6, fixed = TRUE))),
  family = "poisson",
  data = ssf_complete,
  E = weight,
  control.fixed = list(mean = 0, prec = chosen_prior$prec),
  control.compute = list(config = TRUE, mlik = TRUE),
  verbose = FALSE
)

# M2: Full model -- adds sensoryshed pathway
m2 <- inla(
  case_numeric ~ -1 +
    dist_scaled +
    hv_chest_scaled +
    sensoryshed_scaled +
    dist_x_sensoryshed +
    f(stratum_id, model = "iid",
      hyper = list(prec = list(initial = -6, fixed = TRUE))),
  family = "poisson",
  data = ssf_complete,
  E = weight,
  control.fixed = list(mean = 0, prec = chosen_prior$prec),
  control.compute = list(config = TRUE, mlik = TRUE),
  verbose = FALSE
)


# Fixed effects
round(m1$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3)
round(m2$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3)

# Distance coefficient shift between models
dist_m1 <- m1$summary.fixed["dist_scaled", ]
dist_m2 <- m2$summary.fixed["dist_scaled", ]
dist_shift <- dist_m2$mean - dist_m1$mean

# Bayes factor from marginal likelihoods
lml_m1 <- m1$mlik[1]
lml_m2 <- m2$mlik[1]
log_BF <- lml_m2 - lml_m1
# Kass & Raftery (1995): 2*|log BF| > 10 = very strong evidence

# Prior predictive check
#
# Simulate coefficients from the prior and show what they imply
# about relative selection strength (odds ratios).

set.seed(42)
n_sim <- 1000
prior_sd <- 1 / sqrt(chosen_prior$prec)

prior_draws <- data.frame(
  dist        = rnorm(n_sim, 0, prior_sd),
  sensoryshed = rnorm(n_sim, 0, prior_sd),
  hv_chest    = rnorm(n_sim, 0, prior_sd),
  interaction = rnorm(n_sim, 0, prior_sd)
)

# exp(beta) = odds ratio for 1 SD change
prior_or <- exp(prior_draws)

png("prior_predictive_check.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (param in names(prior_draws)) {
  hist(prior_or[[param]], breaks = 50, col = "grey80", border = "grey60",
       main = param, xlab = "Odds ratio (1 SD change)",
       xlim = quantile(prior_or[[param]], c(0.01, 0.99)))
  abline(v = 1, col = "red", lwd = 2, lty = 2)
}
mtext(paste0("Prior predictive check: N(0, sd=", prior_sd, ")"),
      outer = TRUE, line = -1.5, cex = 1.2)
dev.off()

# Prior vs posterior comparison

n_post <- 1000
post_samples <- inla.posterior.sample(n_post, m2)

param_names <- rownames(m2$summary.fixed)
post_draws <- matrix(NA, nrow = n_post, ncol = length(param_names))
colnames(post_draws) <- param_names

for (j in seq_along(param_names)) {
  # INLA latent field names have :1 suffix (e.g. "dist_scaled:1")
  latent_name <- paste0(param_names[j], ":1")
  post_draws[, j] <- sapply(post_samples, function(s) {
    idx <- which(rownames(s$latent) == latent_name)
    if (length(idx) == 1) s$latent[idx, 1] else NA
  })
}

display_names <- c(
  dist_scaled        = "Distance",
  hv_chest_scaled    = "Local Habitat Visibility",
  sensoryshed_scaled = "Sensoryshed",
  dist_x_sensoryshed = "Interaction"
)

png("prior_vs_posterior.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (param in param_names) {
  post_vals <- post_draws[, param]
  prior_vals <- rnorm(n_post, 0, prior_sd)
  xlims <- range(c(
    quantile(post_vals, c(0.005, 0.995)),
    quantile(prior_vals, c(0.005, 0.995))
  ))
  d_post  <- density(post_vals)
  d_prior <- density(prior_vals, from = xlims[1], to = xlims[2])
  plot_title <- ifelse(param %in% names(display_names), display_names[param], param)
  plot(d_prior, col = "grey60", lwd = 2, lty = 2,
       main = plot_title, xlab = "Coefficient", ylab = "Density",
       xlim = xlims, ylim = c(0, max(d_post$y, d_prior$y) * 1.1))
  lines(d_post, col = "steelblue", lwd = 2)
  abline(v = 0, col = "red", lty = 3)
  legend("topright", c("Prior", "Posterior"), col = c("grey60", "steelblue"),
         lwd = 2, lty = c(2, 1), bty = "n", cex = 0.8)
}
dev.off()


# =============================================================================
# Save
# =============================================================================

saveRDS(m1, "inla_m1_baseline.rds")
saveRDS(m2, "inla_m2_full.rds")

sensitivity_fixed <- lapply(sensitivity_models, function(m) {
  as.data.frame(m$summary.fixed)
})

results <- list(
  m1_fixed    = as.data.frame(m1$summary.fixed),
  m2_fixed    = as.data.frame(m2$summary.fixed),
  log_BF      = log_BF,
  lml_m1      = lml_m1,
  lml_m2      = lml_m2,
  dist_shift  = dist_shift,
  chosen_prior_prec = chosen_prior$prec,
  sensitivity = sensitivity_fixed,
  n_obs       = nrow(ssf_complete)
)

saveRDS(results, "model_comparison_results.rds")
View(results)
