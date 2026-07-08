
library(ggplot2)
library(dplyr)


# --- Model coefficients (from INLA output) -----------------------------------
betas <- list(
  dist          = 0.053,
  clearance     = 0.536,
  vis           =  -0.017,
  dist_x_clear  = 0.641
)

# 95% CIs
cis <- list(
  dist          = c(-0.381,   0.487),
  clearance     = c( 0.330,   0.742),
  vis           = c(-0.050,  0.017),
  dist_x_clear  = c(0.419,   0.864)
)

# --- Shared theme -------------------------------------------------------------
theme_marginal <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90"),
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 8)),
      plot.subtitle = element_text(colour = "grey40", size = 10, margin = margin(b = 12)),
      axis.title = element_text(size = 12),
      legend.position = "right",
      plot.margin = margin(15, 15, 10, 15)
    )
}


# Interaction — distance effect at different viewshed clearance levels
# Relative selection strength = exp(marginal_effect × dist_scaled)

# Create prediction grid
clearance_levels <- seq(-2, 2, length.out = 200)
dist_levels      <- seq(-2, 3, length.out = 200)

# (a) Marginal β for distance across clearance levels
marginal_dist <- data.frame(
  clearance_sc = clearance_levels,
  beta_dist    = betas$dist + betas$dist_x_clear * clearance_levels,
  beta_lower   = cis$dist[1] + cis$dist_x_clear[1] * clearance_levels,
  beta_upper   = cis$dist[2] + cis$dist_x_clear[2] * clearance_levels
)

p1a <- ggplot(marginal_dist, aes(x = clearance_sc)) +
  geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper), alpha = 0.15, fill = "#2E86AB") +
  geom_line(aes(y = beta_dist), colour = "#2E86AB", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  labs(
    subtitle = "β > 0 = avoidance (selects steps farther from baboon); β < 0 = attraction",
    x = "Viewshed clearance (scaled)\n← Blocked             Clear →",
    y = expression("Marginal " * beta[distance])
  ) +
  annotate("text", x = -1.5, y = max(marginal_dist$beta_dist) * 0.3,
           label = "Cannot see baboon\n→ weak/no avoidance",
           colour = "grey50", size = 3.2, hjust = 0) +
  annotate("text", x = 1.2, y = max(marginal_dist$beta_dist) * 0.85,
           label = "Can see baboon\n→ strong avoidance",
           colour = "#2E86AB", size = 3.2, hjust = 0) +
  theme_marginal()

# (b) Full RSS surface — relative selection strength across distance × clearance
grid <- expand.grid(
  dist_sc = dist_levels,
  clearance_sc = c(-1.5, -0.5, 0, 0.5, 1.5)
)

grid$clearance_label <- factor(
  grid$clearance_sc,
  levels = c(-1.5, -0.5, 0, 0.5, 1.5),
  labels = c("Very blocked (−1.5 SD)", "Blocked (−0.5 SD)",
             "Mean clearance", "Partial view (+0.5 SD)", "Clear view (+1.5 SD)")
)

# log-RSS = β_dist × dist + β_clear × clearance + β_int × dist × clearance
grid$log_rss <- betas$dist * grid$dist_sc +
  betas$clearance * grid$clearance_sc +
  betas$dist_x_clear * grid$dist_sc * grid$clearance_sc

grid$rss <- exp(grid$log_rss)

p1b <- ggplot(grid, aes(x = dist_sc, y = log_rss, colour = clearance_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(
    values = c("#D32F2F", "#FF8F00", "grey50", "#43A047", "#1565C0"),
    name = "Viewshed clearance"
  ) +
  labs(x = "Distance (scaled)\n←Closer                     Farther →",
       y = "log-RSS (relative selection strength)"
  ) +
  theme_marginal() +
  theme(legend.position = "right")
p1a
p1b


# Coefficient summary (forest plot)

coef_df <- data.frame(
  param = c("Distance", "Viewshed clearance", "Horizontal visibility",
            "Distance × Clearance"),
  beta  = c(betas$dist, betas$clearance,
            betas$vis, betas$dist_x_clear),
  lower = c(cis$dist[1], cis$clearance[1],
            cis$vis[1], cis$dist_x_clear[1]),
  upper = c(cis$dist[2], cis$clearance[2],
            cis$vis[2], cis$dist_x_clear[2])
)

coef_df$param <- factor(coef_df$param, levels = rev(coef_df$param))
coef_df$type <- c("Distance", "Viewshed",  "Habitat Visibility", "Interaction")

p4 <- ggplot(coef_df, aes(x = beta, y = param, colour = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_colour_manual(
    values = c(Distance = "#2E86AB", Visual = "#43A047",
               Auditory = "#7B1FA2", Habitat = "#00695C",
               Interaction = "#D32F2F"),
    guide = "none"
  ) +
  labs(
    x = expression(beta * " (log-RSS)"),
    y = NULL
  ) +
  theme_marginal() +
  theme(axis.text.y = element_text(size = 11))
p4
