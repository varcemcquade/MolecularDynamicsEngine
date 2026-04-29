# msd_analysis.R
# Reads msd.csv produced by Box::MeanSquaredDisplacement() and:
#   1. Plots MSD vs. time
#   2. Fits the Einstein relation MSD = 6*D*t to extract diffusion coefficient D
#   3. Compares D to published experimental argon diffusion data
#
# Einstein relation: <|r(t) - r(0)|²> = 6*D*t  (3D, long-time limit)
# So D = slope / 6, where slope = d(MSD)/dt
#
# Experimental argon self-diffusion near 94 K (triple point ~84 K):
#   D ≈ 1.7–2.4 × 10⁻⁹ m²/s  (from Naghizadeh & Rice, J. Chem. Phys. 1962)
#
# Usage: Rscript msd_analysis.R

library(ggplot2)
library(dplyr)

DATA_PATH    <- "msd.csv"
PLOT_DIR     <- "plots"
DT_SECONDS   <- 1e-15   # timestep in seconds (must match main.cpp dt)
NSAVC        <- 100     # save frequency (must match main.cpp nsavc)
# Time per frame in seconds
DT_FRAME     <- DT_SECONDS * NSAVC

dir.create(PLOT_DIR, showWarnings = FALSE)

df <- read.csv(DATA_PATH, header = FALSE, col.names = c("step", "MSD_angstroms2"))
cat("Loaded", nrow(df), "MSD frames from", DATA_PATH, "\n")

df$time_s    <- df$step * DT_SECONDS        # seconds
df$MSD_m2    <- df$MSD_angstroms2 * 1e-20  # Å² → m²

cat(sprintf("Time range: %.2e – %.2e seconds\n\n",
            min(df$time_s), max(df$time_s)))

  geom_line(color = "#0072B2", linewidth = 0.8) +
  labs(title    = "Mean Squared Displacement vs. Time",
       subtitle = "LJ Argon — NVE ensemble",
       x        = "Time (ps)",
       y        = expression("MSD (Å"^2*")")) +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "msd_timeseries.png"),
       p_msd, width = 8, height = 4, dpi = 150)
cat("Saved msd_timeseries.png\n")

cutoff_idx <- max(1, round(nrow(df) * 0.20))
df_linear  <- df[cutoff_idx:nrow(df), ]

model  <- lm(MSD_m2 ~ time_s, data = df_linear)
s      <- summary(model)
slope  <- coef(model)[["time_s"]]
slope_se <- s$coefficients["time_s", "Std. Error"]
D      <- slope / 6
D_se   <- slope_se / 6
ci     <- confint(model)["time_s", ]
D_ci   <- ci / 6

cat("─── Diffusion Coefficient (Einstein Relation: MSD = 6Dt) ────────────────\n")
cat(sprintf("  Fit slope (d(MSD)/dt): %.4e ± %.4e m²/s\n", slope, slope_se))
cat(sprintf("  D = slope/6:           %.4e ± %.4e m²/s\n", D,     D_se))
cat(sprintf("  95%% CI for D:          [%.4e, %.4e] m²/s\n", D_ci[1], D_ci[2]))
cat(sprintf("  R²:                    %.6f\n", s$r.squared))
cat(sprintf("\n  Experimental argon near 94 K: D ≈ 1.7–2.4 × 10⁻⁹ m²/s\n"))
cat(sprintf("  Ratio (sim/exp_midpoint): %.3f\n\n", D / 2.05e-9))

df$MSD_fit <- predict(model, newdata = df)

p_fit <- ggplot(df, aes(x = time_s * 1e12)) +
  geom_line(aes(y = MSD_m2 * 1e18), color = "#0072B2",
            linewidth = 0.8, alpha = 0.7) +
  geom_line(aes(y = MSD_fit * 1e18), color = "#CC3311",
            linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = max(df$time_s) * 1e12 * 0.05,
           y = max(df$MSD_m2 * 1e18, na.rm = TRUE) * 0.85,
           label = sprintf("D = %.3e m²/s\nR² = %.4f", D, s$r.squared),
           hjust = 0, color = "#CC3311", size = 4) +
  labs(title    = "MSD — Einstein Relation Fit",
       subtitle = "Red dashed = linear fit on diffusive regime (t > 20% of run)",
       x = "Time (ps)",
       y = expression("MSD (Å"^2*")")) +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "msd_einstein_fit.png"),
       p_fit, width = 8, height = 4, dpi = 150)
cat("Saved msd_einstein_fit.png\n")

df_pos       <- df[df$MSD_m2 > 0 & df$time_s > 0, ]
df_pos$log_t <- log10(df_pos$time_s)
df_pos$log_M <- log10(df_pos$MSD_m2)

p_loglog <- ggplot(df_pos, aes(x = log_t, y = log_M)) +
  geom_line(color = "#0072B2", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = log10(6 * D) + 0,
              color = "#009E73", linewidth = 0.7, linetype = "dotted") +
  geom_abline(slope = 2, intercept = min(df_pos$log_M),
              color = "#E69F00", linewidth = 0.7, linetype = "dotted") +
  annotate("text", x = min(df_pos$log_t) + 0.5,
           y = min(df_pos$log_M) + 1.5,
           label = "slope=2\n(ballistic)", color = "#E69F00", size = 3.3) +
  annotate("text", x = max(df_pos$log_t) - 0.5,
           y = max(df_pos$log_M) - 0.5,
           label = "slope=1\n(diffusive)", color = "#009E73", size = 3.3) +
  labs(title    = "MSD — Log-Log Scale (Ballistic → Diffusive Crossover)",
       x        = "log₁₀(time / s)",
       y        = "log₁₀(MSD / m²)") +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "msd_loglog.png"),
       p_loglog, width = 8, height = 4, dpi = 150)
cat("Saved msd_loglog.png\n")