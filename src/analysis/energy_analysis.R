# energy_analysis.R
# Reads energies.csv produced by Box::Equilibrate() and performs:
#   1. Time-series plots of KE, PE, and total energy
#   2. Linear regression on total energy to test for systematic drift
#   3. Summary statistics (mean, SD, coefficient of variation)
#
# Usage: Rscript energy_analysis.R
# Place this script in the same directory as energies.csv, or update DATA_PATH below.

library(ggplot2)
library(dplyr)

DATA_PATH   <- "energies.csv"   # path to engine output
PLOT_DIR    <- "plots"          # output directory for figures
TARGET_TEMP <- 94.4             # K — Argon simulation target

dir.create(PLOT_DIR, showWarnings = FALSE)

df <- read.csv(DATA_PATH)
cat("Loaded", nrow(df), "frames from", DATA_PATH, "\n")
cat("Steps:", min(df$step), "to", max(df$step), "\n\n")

df_long <- tidyr::pivot_longer(df,
  cols      = c(KE, PE, E_total),
  names_to  = "component",
  values_to = "energy"
)
df_long$component <- factor(df_long$component,
  levels = c("KE", "PE", "E_total"),
  labels = c("Kinetic", "Potential", "Total")
)

p_energy <- ggplot(df_long, aes(x = step, y = energy, color = component)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  facet_wrap(~component, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Kinetic" = "#E69F00",
                                 "Potential" = "#56B4E9",
                                 "Total"     = "#009E73")) +
  labs(title    = "MD Simulation — Energy vs. Time",
       subtitle = paste("Lennard-Jones Argon |", nrow(df), "frames"),
       x        = "Timestep",
       y        = "Energy (J)",
       color    = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text       = element_text(face = "bold"))

ggsave(file.path(PLOT_DIR, "energy_timeseries.png"),
       p_energy, width = 8, height = 7, dpi = 150)
cat("Saved energy_timeseries.png\n")

model <- lm(E_total ~ step, data = df)
s     <- summary(model)

slope    <- coef(model)[["step"]]
slope_se <- s$coefficients["step", "Std. Error"]
p_val    <- s$coefficients["step", "Pr(>|t|)"]
r_sq     <- s$r.squared
ci       <- confint(model)["step", ]

cat("─── Energy Drift Test (lm: E_total ~ step) ──────────────────────────────\n")
cat(sprintf("  Slope:        %+.4e J/step\n",  slope))
cat(sprintf("  Std. Error:   %.4e J/step\n",   slope_se))
cat(sprintf("  p-value:      %.4f %s\n", p_val,
    ifelse(p_val < 0.05, "(DRIFT DETECTED — check timestep or cutoff)", "(no significant drift)")))
cat(sprintf("  R²:           %.6f\n", r_sq))
cat(sprintf("  95%% CI:       [%+.4e, %+.4e] J/step\n\n", ci[1], ci[2]))

p_drift <- ggplot(df, aes(x = step, y = E_total)) +
  geom_line(color = "#009E73", linewidth = 0.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#CC3311",
              linewidth = 0.8, linetype = "dashed") +
  labs(title    = "Total Energy — Drift Detection",
       subtitle = sprintf("slope = %+.3e J/step  |  p = %.4f  |  R² = %.5f",
                          slope, p_val, r_sq),
       x = "Timestep", y = "Total Energy (J)") +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "energy_drift.png"),
       p_drift, width = 8, height = 4, dpi = 150)
cat("Saved energy_drift.png\n")

cat("─── Thermodynamic Averages ───────────────────────────────────────────────\n")
stats <- df %>%
  summarise(
    KE_mean     = mean(KE),   KE_sd     = sd(KE),
    PE_mean     = mean(PE),   PE_sd     = sd(PE),
    E_tot_mean  = mean(E_total), E_tot_sd = sd(E_total),
    T_mean      = mean(Temperature), T_sd = sd(Temperature)
  )

cat(sprintf("  <KE>        = %.6e ± %.4e J\n", stats$KE_mean,    stats$KE_sd))
cat(sprintf("  <PE>        = %.6e ± %.4e J\n", stats$PE_mean,    stats$PE_sd))
cat(sprintf("  <E_total>   = %.6e ± %.4e J\n", stats$E_tot_mean, stats$E_tot_sd))
cat(sprintf("  <T>         = %.2f ± %.2f K  (target: %.1f K)\n",
            stats$T_mean, stats$T_sd, TARGET_TEMP))
cat(sprintf("  ΔT / T_target = %.2f%%\n\n",
            abs(stats$T_mean - TARGET_TEMP) / TARGET_TEMP * 100))