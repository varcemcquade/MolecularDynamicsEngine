# temperature_check.R
# Validates simulated temperature against target and checks for equilibration.
# Uses the instantaneous temperature column written to energies.csv.
#
# Tests performed:
#   1. Temperature time-series plot with target line
#   2. One-sample t-test: is <T> statistically different from target?
#   3. Equilibration check: does temperature stabilize (rolling average)?
#   4. Distribution plot of instantaneous temperatures
#
# Usage: Rscript temperature_check.R

library(ggplot2)
library(dplyr)

DATA_PATH   <- "energies.csv"
PLOT_DIR    <- "plots"
TARGET_TEMP <- 94.4     # K
WINDOW      <- 20       # rolling average window (frames)

dir.create(PLOT_DIR, showWarnings = FALSE)

df <- read.csv(DATA_PATH)

df$T_rolling <- stats::filter(df$Temperature,
                               rep(1 / WINDOW, WINDOW), sides = 1)

p_temp <- ggplot(df, aes(x = step)) +
  geom_line(aes(y = Temperature), color = "#999999", linewidth = 0.4, alpha = 0.6) +
  geom_line(aes(y = T_rolling),   color = "#E69F00", linewidth = 0.9, na.rm = TRUE) +
  geom_hline(yintercept = TARGET_TEMP, color = "#CC3311",
             linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = max(df$step) * 0.02, y = TARGET_TEMP + 1.5,
           label = paste("Target:", TARGET_TEMP, "K"),
           color = "#CC3311", hjust = 0, size = 3.5) +
  labs(title    = "Instantaneous Temperature vs. Time",
       subtitle = paste("Grey = raw  |  Orange = rolling mean (window =", WINDOW, "frames)"),
       x = "Timestep", y = "Temperature (K)") +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "temperature_timeseries.png"),
       p_temp, width = 8, height = 4, dpi = 150)
cat("Saved temperature_timeseries.png\n")

t_result <- t.test(df$Temperature, mu = TARGET_TEMP)

cat("─── One-Sample t-Test: <T> vs. Target ──────────────────────────────────\n")
cat(sprintf("  Mean T:    %.4f K\n",   mean(df$Temperature)))
cat(sprintf("  SD T:      %.4f K\n",   sd(df$Temperature)))
cat(sprintf("  Target T:  %.4f K\n",   TARGET_TEMP))
cat(sprintf("  t-stat:    %.4f\n",     t_result$statistic))
cat(sprintf("  p-value:   %.4f %s\n",  t_result$p.value,
    ifelse(t_result$p.value < 0.05,
           "(significant deviation from target — check thermostat or temperature formula)",
           "(mean temperature consistent with target)")))
cat(sprintf("  95%% CI:    [%.4f, %.4f] K\n\n",
            t_result$conf.int[1], t_result$conf.int[2]))

p_dist <- ggplot(df, aes(x = Temperature)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "#56B4E9", color = "white", alpha = 0.8) +
  geom_density(color = "#0072B2", linewidth = 1) +
  geom_vline(xintercept = TARGET_TEMP, color = "#CC3311",
             linewidth = 0.9, linetype = "dashed") +
  geom_vline(xintercept = mean(df$Temperature), color = "#009E73",
             linewidth = 0.9, linetype = "solid") +
  annotate("text", x = TARGET_TEMP + 0.5, y = Inf, vjust = 2,
           label = "Target", color = "#CC3311", size = 3.5) +
  annotate("text", x = mean(df$Temperature) + 0.5, y = Inf, vjust = 4,
           label = "Mean", color = "#009E73", size = 3.5) +
  labs(title    = "Distribution of Instantaneous Temperatures",
       subtitle = sprintf("N = %d frames  |  mean = %.2f K  |  SD = %.2f K",
                          nrow(df), mean(df$Temperature), sd(df$Temperature)),
       x = "Temperature (K)", y = "Density") +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "temperature_distribution.png"),
       p_dist, width = 7, height = 4, dpi = 150)
cat("Saved temperature_distribution.png\n")

n_half   <- nrow(df) %/% 2
T_first  <- mean(df$Temperature[1:n_half])
T_second <- mean(df$Temperature[(n_half + 1):nrow(df)])
eq_test  <- t.test(df$Temperature[1:n_half], df$Temperature[(n_half + 1):nrow(df)])

cat("─── Equilibration Check (first half vs. second half) ────────────────────\n")
cat(sprintf("  First  half <T>: %.4f K\n", T_first))
cat(sprintf("  Second half <T>: %.4f K\n", T_second))
cat(sprintf("  p-value:         %.4f %s\n", eq_test$p.value,
    ifelse(eq_test$p.value < 0.05,
           "(halves differ — system may not be fully equilibrated)",
           "(halves statistically equivalent — system appears equilibrated)")))