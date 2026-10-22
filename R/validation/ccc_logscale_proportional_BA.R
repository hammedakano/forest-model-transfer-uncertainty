# =============================================================================
# AGB agreement: log-scale CCC + proportional (percentage) Bland-Altman
# Companion to the raw-scale epi.ccc() already run on clean_compare
# Assumes: clean_compare has columns `agb` (field) and `est_agb` (UAV)
#          and that epiR is loaded (library(epiR))
# =============================================================================

library(epiR)

library(epiR)

# Clean out any NA or missing values from your comparison columns first
clean_compare <- all_data[!is.na(all_data$agb) & !is.na(all_data$est_agb), ]

# Calculate Concordance Correlation Coefficient
ccc_results <- epi.ccc(clean_compare$agb, clean_compare$est_agb, ci = "z-transform", conf.level = 0.95)

# Extract the core estimate value
ccc_value <- round(ccc_results$rho.c$est, 3)
# Extract confidence interval bounds
ccc_ci_lower <- round(ccc_results$rho.c$lower, 3)
ccc_ci_upper <- round(ccc_results$rho.c$upper, 3)
cat("Concordance Correlation Coefficient (CCC):", ccc_value,
    "(95% CI", ccc_ci_lower, "-", ccc_ci_upper, ")\n")
cat("  C.b (accuracy/bias factor):", round(ccc_results$C.b, 4), "\n")
cat("  scale shift:", round(ccc_results$s.shift, 3),
    "| location shift:", round(ccc_results$l.shift, 4), "\n\n")

# Work on a clean, complete-case copy (defensive; you already filtered NAs)
cc <- clean_compare[!is.na(clean_compare$agb) & !is.na(clean_compare$est_agb), ]

# -----------------------------------------------------------------------------
# 0. Data-integrity check: are there duplicated (field, UAV) pairs?
#    Duplicates inflate n and artificially narrow the CCC confidence interval.
# -----------------------------------------------------------------------------
n_dup <- sum(duplicated(cc[, c("agb", "est_agb")]))
cat("Duplicated (agb, est_agb) pairs:", n_dup, "of", nrow(cc), "rows\n\n")
# If n_dup is large and unexpected, inspect with:
# cc[duplicated(cc[, c("agb","est_agb")]) |
#    duplicated(cc[, c("agb","est_agb")], fromLast = TRUE), ]

# AGB must be strictly positive to log-transform
if (any(cc$agb <= 0 | cc$est_agb <= 0)) {
  warning("Non-positive AGB values present; log transform will drop or fail on these.")
  cc <- cc[cc$agb > 0 & cc$est_agb > 0, ]
}

# -----------------------------------------------------------------------------
# 1. CCC on the log scale
#    Better reflects agreement when error is multiplicative (scatter grows with size)
# -----------------------------------------------------------------------------
ccc_log <- epi.ccc(log(cc$agb), log(cc$est_agb),
                   ci = "z-transform", conf.level = 0.95)

cat("---- Raw-scale vs log-scale CCC ----\n")
cat("Log-scale CCC :", round(ccc_log$rho.c$est, 3),
    "(95% CI", round(ccc_log$rho.c$lower, 3), "-",
    round(ccc_log$rho.c$upper, 3), ")\n")
cat("  C.b (accuracy/bias factor):", round(ccc_log$C.b, 4), "\n")
cat("  scale shift:", round(ccc_log$s.shift, 3),
    "| location shift:", round(ccc_log$l.shift, 4), "\n\n")

# -----------------------------------------------------------------------------
# 2. Formal test for heteroscedasticity (proportional error)
#    Regress |difference| on the pair mean: a significant positive slope
#    confirms the funnel pattern seen in the raw Bland-Altman.
# -----------------------------------------------------------------------------
m       <- (cc$agb + cc$est_agb) / 2          # mean of each pair
raw_d   <- cc$agb - cc$est_agb                # field - UAV (positive => UAV underestimates)
het_fit <- lm(abs(raw_d) ~ m)
cat("---- Heteroscedasticity test: |difference| ~ mean ----\n")
print(summary(het_fit)$coefficients)
cat("Spearman rho(|difference|, mean):",
    round(cor(abs(raw_d), m, method = "spearman"), 3), "\n\n")

# -----------------------------------------------------------------------------
# 3. Proportional (percentage) Bland-Altman
#    Difference expressed as % of the pair mean -> stable across the size range
# -----------------------------------------------------------------------------
d_pct     <- 100 * raw_d / m
bias_pct  <- mean(d_pct)
sd_pct    <- sd(d_pct)
loa_low   <- bias_pct - 1.96 * sd_pct
loa_high  <- bias_pct + 1.96 * sd_pct

cat("---- Proportional Bland-Altman (% of mean) ----\n")
cat("Mean bias        :", round(bias_pct, 2), "%\n")
cat("SD of % diff     :", round(sd_pct, 2), "%\n")
cat("95% LoA          :", round(loa_low, 2), "% to", round(loa_high, 2), "%\n")
cat("(field - UAV; positive bias => UAV underestimates on average)\n\n")

# -----------------------------------------------------------------------------
# 4. Proportional Bland-Altman plot (base R, no extra dependencies)
# -----------------------------------------------------------------------------
plot(m, d_pct,
     xlab = "Mean of field & UAV AGB (kg)",
     ylab = "Difference (% of mean)",
     main = "Proportional Bland-Altman",
     pch = 19, col = adjustcolor("steelblue", 0.6))
abline(h = bias_pct, lwd = 2)
abline(h = c(loa_low, loa_high), lty = 2)
abline(h = 0, col = "grey60")
legend("topright", bty = "n", lty = c(1, 2),
       legend = c(paste0("bias = ", round(bias_pct, 1), "%"),
                  paste0("95% LoA = +/-", round(1.96 * sd_pct, 1), "%")))
