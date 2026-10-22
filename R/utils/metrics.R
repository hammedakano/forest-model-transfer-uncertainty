# Agreement and out-of-sample error metrics used across the repository.
# All functions take observed and predicted numeric vectors of equal length.

# Lin's concordance correlation coefficient. Penalises bias and scale shifts,
# not only scatter. Population variances (/ n) keep the formula consistent.
ccc <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  mo <- mean(obs); mp <- mean(pred)
  so <- mean((obs - mo)^2); sp <- mean((pred - mp)^2)
  cov_op <- mean((obs - mo) * (pred - mp))
  (2 * cov_op) / (so + sp + (mo - mp)^2)
}

# One-row validation panel.
#   r    Pearson correlation (association only)
#   R2   1 - SSE/SST on held-out data (can be negative out of sample)
#   RMSE root mean squared error; rRMSE as % of the observed mean
#   bias mean(pred - obs); positive means over-prediction
#   CCC  Lin's concordance (headline agreement metric)
metrics_oos <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((pred - obs)^2); sst <- sum((obs - mean(obs))^2)
  rmse <- sqrt(mean((pred - obs)^2))
  data.frame(
    n = length(obs),
    r = stats::cor(obs, pred),
    R2 = 1 - sse / sst,
    RMSE = rmse,
    rRMSE_pct = 100 * rmse / mean(obs),
    bias = mean(pred - obs),
    CCC = ccc(obs, pred)
  )
}

# Proportional Bland-Altman summary: differences as % of the pair mean.
bland_altman_pct <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred) & (obs + pred) > 0
  m <- (obs[ok] + pred[ok]) / 2
  d <- 100 * (pred[ok] - obs[ok]) / m
  data.frame(mean_diff_pct = mean(d),
             loa_lower_pct = mean(d) - 1.96 * stats::sd(d),
             loa_upper_pct = mean(d) + 1.96 * stats::sd(d))
}
