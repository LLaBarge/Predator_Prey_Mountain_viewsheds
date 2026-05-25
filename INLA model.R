
# run INLA model

# load data 
# ver_ssf_data <- read.csv("ver_ssf_data_viewshed.csv", stringsAsFactors = FALSE)

# Prepare for INLA
ssf_complete <- ver_ssf_data %>%
  filter(!is.na(dist_to_baboon),
         !is.na(viewshed_clearance_m))

ssf_complete$dist_scaled      <- scale(ssf_complete$dist_to_baboon)[, 1]
ssf_complete$clearance_scaled <- scale(ssf_complete$viewshed_clearance_m)[, 1]
ssf_complete$vis_scaled       <- scale(ssf_complete$horizontal_visibility_m)[, 1]


ssf_complete <- ssf_complete %>%
  group_by(step_id_) %>%
  mutate(
    n_steps_in_stratum = n(),
    weight = 1e6 / n_steps_in_stratum
  ) %>%
  ungroup()

# create strata for SSF
ssf_complete$case_numeric <- as.integer(ssf_complete$case_)
ssf_complete$stratum_id   <- as.factor(ssf_complete$step_id_)


# interaction effects
ssf_complete$dist_x_vis <- ssf_complete$dist_scaled * ssf_complete$vis_scaled
ssf_complete$dist_x_clear<- ssf_complete$dist_scaled * ssf_complete$clearance_scaled



# =============================================================================
# Fit INLA model
global_model <- inla(
  case_numeric ~ -1 +
    dist_scaled +
    clearance_scaled +
    vis_scaled +
    dist_x_clear + dist_x_vis +
    f(stratum_id, model = "iid",
      hyper = list(prec = list(initial = -6, fixed = TRUE))),
  family = "poisson",
  data = ssf_complete,
  E = weight,
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  verbose = FALSE
)

# =============================================================================
# Results
fixed <- global_model$summary.fixed

for (i in 1:nrow(fixed)) {
  param   <- rownames(fixed)[i]
  beta    <- fixed[i, "mean"]
  ci_low  <- fixed[i, "0.025quant"]
  ci_high <- fixed[i, "0.975quant"]
  sig     <- if (ci_low > 0 | ci_high < 0) " *" else ""
  cat(sprintf("%-25s B = %7.3f  (95%% CI: [%7.3f, %7.3f])%s\n",
              param, beta, ci_low, ci_high, sig))
}

cat(sprintf("\nDIC:  %.1f\n", global_model$dic$dic))
cat(sprintf("WAIC: %.1f\n", global_model$waic$waic))

results_summary <- list(
  fixed_effects = as.data.frame(fixed),
  dic  = global_model$dic$dic,
  waic = global_model$waic$waic,
  n_parameters = global_model$dic$p.eff,
  n_obs = nrow(ssf_complete)
)

saveRDS(results_summary, "model_results.rds")