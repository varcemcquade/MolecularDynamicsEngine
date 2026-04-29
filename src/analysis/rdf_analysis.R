# rdf_analysis.R
# Reads rdf.csv produced by Box::PairDistributionFunction() and:
#   1. Plots g(r) vs. r with peak annotation
#   2. Locates the first coordination shell peak and compares to LJ theory
#   3. Computes the coordination number (integral of first peak)
#
# For Lennard-Jones argon near the triple point (~94 K), expect:
#   - First peak at r ≈ 3.7–3.9 Å  (≈ 1.08–1.15 σ where σ = 3.4 Å)
#   - g(r) → 1 at large r (ideal gas limit)
#   - Coordination number ≈ 12 (close-packed-like for liquid near Tm)
#
# Usage: Rscript rdf_analysis.R

library(ggplot2)
library(dplyr)

DATA_PATH <- "rdf.csv"
PLOT_DIR  <- "plots"
SIGMA     <- 3.4    # Angstroms — LJ sigma for argon

dir.create(PLOT_DIR, showWarnings = FALSE)

df <- read.csv(DATA_PATH)
cat("Loaded", nrow(df), "radial bins from", DATA_PATH, "\n")

df$r_reduced <- df$r / SIGMA

first_peak_idx <- which.max(df$g_r)
r_peak         <- df$r[first_peak_idx]
g_peak         <- df$g_r[first_peak_idx]

p_rdf <- ggplot(df, aes(x = r_reduced, y = g_r)) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  geom_hline(yintercept = 1.0, color = "#999999",
             linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 4.5, y = 1.04,
           label = "g(r) = 1  (ideal gas limit)", color = "#999999", size = 3.3) +
  geom_vline(xintercept = r_peak / SIGMA, color = "#E69F00",
             linewidth = 0.7, linetype = "dotted") +
  annotate("point", x = r_peak / SIGMA, y = g_peak,
           color = "#CC3311", size = 3) +
  annotate("text",  x = r_peak / SIGMA + 0.15, y = g_peak,
           label  = sprintf("r₁ = %.2f σ\ng(r₁) = %.2f", r_peak / SIGMA, g_peak),
           color  = "#CC3311", hjust = 0, size = 3.5) +
  labs(title    = "Pair Distribution Function g(r)",
       subtitle = sprintf("LJ Argon | σ = %.1f Å | First peak at r = %.3f Å (%.3f σ)",
                          SIGMA, r_peak, r_peak / SIGMA),
       x = expression("r / " * sigma),
       y = "g(r)") +
  xlim(0, min(6, max(df$r_reduced))) +
  theme_minimal(base_size = 12)

ggsave(file.path(PLOT_DIR, "rdf.png"), p_rdf, width = 8, height = 4, dpi = 150)
cat("Saved rdf.png\n")

df_after_peak <- df[(first_peak_idx + 1):nrow(df), ]
min_after_peak <- which.min(df_after_peak$g_r[1:min(50, nrow(df_after_peak))])
r_min2_idx    <- first_peak_idx + min_after_peak

r_shell <- df$r[1:r_min2_idx]
g_shell <- df$g_r[1:r_min2_idx]
dr      <- diff(df$r)[1]    # bin width (uniform)

BOX_LENGTH_APPROX <- 34.5   # Angstroms — update if ATOM_COUNT or RHO changes
rho_number        <- 1000 / BOX_LENGTH_APPROX^3   # atoms/Å³

coord_number <- 4 * pi * rho_number * sum(g_shell * r_shell^2 * dr)

cat("─── Coordination Number ─────────────────────────────────────────────────\n")
cat(sprintf("  First peak:         r = %.3f Å  (%.3f σ),  g(r) = %.3f\n",
            r_peak, r_peak / SIGMA, g_peak))
cat(sprintf("  First minimum:      r = %.3f Å\n", df$r[r_min2_idx]))
cat(sprintf("  Coordination number N₁ ≈ %.2f atoms\n", coord_number))
cat(sprintf("  (FCC expected ≈ 12, liquid near Tₘ typically 10–13)\n\n"))

large_r    <- df$g_r[df$r_reduced > 4.0]
mean_large <- mean(large_r)
sd_large   <- sd(large_r)

cat("─── Long-range limit (r > 4σ) ──────────────────────────────────────────\n")
cat(sprintf("  mean g(r) = %.4f  (should be ≈ 1.000)\n", mean_large))
cat(sprintf("  SD  g(r) = %.4f  (should be small)\n", sd_large))
cat(sprintf("  Deviation from 1: %.4f%%\n\n",
            abs(mean_large - 1.0) * 100))