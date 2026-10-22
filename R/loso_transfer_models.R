# Site-transfer modelling of tree or plot attributes (DBH, AGB, height).
#
# Question: does a model trained on some sites predict a site it has never seen?
# Random train/test splits mix trees from the same site into both sets and
# overstate skill. Here every model is scored with leave-one-site-out (LOSO)
# cross-validation, the same structure a forest inventory model faces when it
# is moved to a new stand, species, sensor or region.
#
# Models compared under the same LOSO folds
#   lm      linear baseline
#   rf      random forest (ranger), with quantile prediction intervals
#   xgb     gradient boosting (xgboost)
#   stack   non-negative linear stack of the three, fitted on out-of-fold
#           predictions from the TRAINING sites only (no leakage)
#
# Outputs (in out_dir)
#   loso_predictions.csv      held-out prediction for every row and model
#   loso_metrics_pooled.csv   pooled held-out metrics by model
#   loso_metrics_by_site.csv  held-out metrics by model and site
#   loso_sse_share.csv        share of total error contributed by each site
#   loso_interval_coverage.csv  empirical coverage of RF 90% intervals
#
# Run on the synthetic example:
#   Rscript examples/simulate_example_data.R
#   Rscript R/loso_transfer_models.R
# Run on your data:
#   LOSO_DATA=path.csv LOSO_RESPONSE=agb LOSO_PREDICTORS=Z,crown_area_m2 Rscript ...

suppressPackageStartupMessages({
  library(dplyr); library(ranger); library(xgboost)
})
source(file.path("R", "utils", "metrics.R"))

cfg <- list(
  data       = Sys.getenv("LOSO_DATA", unset = "examples/example_tree_data.csv"),
  response   = Sys.getenv("LOSO_RESPONSE", unset = "agb"),
  predictors = strsplit(Sys.getenv("LOSO_PREDICTORS", unset = "Z,crown_area_m2"), ",")[[1]],
  site_col   = Sys.getenv("LOSO_SITE", unset = "site"),
  out_dir    = Sys.getenv("LOSO_OUT", unset = "outputs/loso"),
  seed       = 9560,
  pi_level   = 0.90
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(cfg$seed)

df <- read.csv(cfg$data)
need <- c(cfg$response, cfg$predictors, cfg$site_col)
miss <- setdiff(need, names(df)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
df <- df[stats::complete.cases(df[, need]), ]
df$.site <- factor(df[[cfg$site_col]])
form <- stats::reformulate(cfg$predictors, response = cfg$response)
cat("Rows per site:\n"); print(table(df$.site))

# ---- Model fitters: each returns a predict(newdata) function ------------------
fit_lm <- function(tr) { m <- lm(form, data = tr); function(nd) as.numeric(predict(m, nd)) }
fit_rf <- function(tr) {
  m <- ranger(form, data = tr, num.trees = 500, min.node.size = 5,
              quantreg = TRUE, seed = cfg$seed)
  function(nd, type = "mean") {
    if (type == "mean") return(predict(m, nd)$predictions)
    a <- (1 - cfg$pi_level) / 2
    predict(m, nd, type = "quantiles", quantiles = c(a, 1 - a))$predictions
  }
}
fit_xgb <- function(tr) {
  X <- as.matrix(tr[, cfg$predictors]); y <- tr[[cfg$response]]
  # xgb.train with a DMatrix works across old and new xgboost releases.
  m <- xgboost::xgb.train(params = list(objective = "reg:squarederror", eta = 0.05,
                                        max_depth = 3, subsample = 0.8, nthread = 1),
                          data = xgboost::xgb.DMatrix(X, label = y), nrounds = 300)
  function(nd) as.numeric(predict(m, xgboost::xgb.DMatrix(as.matrix(nd[, cfg$predictors]))))
}
base_fitters <- list(lm = fit_lm, rf = fit_rf, xgb = fit_xgb)

# Stack weights from inner LOSO within the training sites only.
fit_stack <- function(tr) {
  inner <- matrix(NA_real_, nrow(tr), length(base_fitters), dimnames = list(NULL, names(base_fitters)))
  for (s in unique(tr$.site)) {
    i <- tr$.site == s
    for (b in names(base_fitters)) inner[i, b] <- base_fitters[[b]](tr[!i, ])(tr[i, ])
  }
  # Convex weights (non-negative, sum to 1) that minimise inner-LOSO squared error.
  y <- tr[[cfg$response]]
  sse <- function(lw) { w <- exp(lw) / sum(exp(lw)); sum((y - inner %*% w)^2) }
  lw <- stats::optim(rep(0, ncol(inner)), sse, method = "BFGS")$par
  w <- exp(lw) / sum(exp(lw)); names(w) <- names(base_fitters)
  fitted <- lapply(base_fitters, function(f) f(tr))
  structure(function(nd) as.numeric(sapply(names(fitted), function(b) fitted[[b]](nd)) %*% w),
            weights = w)
}

# ---- Outer LOSO --------------------------------------------------------------
pred <- data.frame(row = seq_len(nrow(df)), site = df$.site, obs = df[[cfg$response]])
for (b in c(names(base_fitters), "stack")) pred[[b]] <- NA_real_
pred$rf_lo <- NA_real_; pred$rf_hi <- NA_real_
stack_weights <- list()

for (s in levels(df$.site)) {
  te <- df$.site == s; tr <- df[!te, ]
  if (length(unique(tr$.site)) < 2) stop("Need at least three sites for LOSO with stacking.")
  for (b in names(base_fitters)) {
    f <- base_fitters[[b]](tr)
    pred[te, b] <- f(df[te, ])
    if (b == "rf") { q <- f(df[te, ], type = "q"); pred$rf_lo[te] <- q[, 1]; pred$rf_hi[te] <- q[, 2] }
  }
  st <- fit_stack(tr); pred$stack[te] <- st(df[te, ])
  stack_weights[[s]] <- attr(st, "weights")
  message("Held out ", s)
}

models <- c(names(base_fitters), "stack")
pooled <- do.call(rbind, lapply(models, function(m) cbind(model = m, metrics_oos(pred$obs, pred[[m]]))))
by_site <- do.call(rbind, lapply(models, function(m)
  pred %>% group_by(site) %>% group_modify(~ metrics_oos(.x$obs, .x[[m]])) %>%
    ungroup() %>% mutate(model = m, .before = 1)))
sse_share <- do.call(rbind, lapply(models, function(m)
  pred %>% mutate(se = (.data[[m]] - obs)^2) %>% group_by(site) %>%
    summarise(SSE = sum(se), n = dplyr::n(), .groups = "drop") %>%
    mutate(model = m, pct_of_total_SSE = round(100 * SSE / sum(SSE), 1))))
coverage <- pred %>% group_by(site) %>%
  summarise(n = dplyr::n(), coverage = mean(obs >= rf_lo & obs <= rf_hi),
            mean_width = mean(rf_hi - rf_lo), .groups = "drop") %>%
  bind_rows(data.frame(site = "ALL", n = nrow(pred),
                       coverage = mean(pred$obs >= pred$rf_lo & pred$obs <= pred$rf_hi),
                       mean_width = mean(pred$rf_hi - pred$rf_lo)))

write.csv(pred, file.path(cfg$out_dir, "loso_predictions.csv"), row.names = FALSE)
write.csv(pooled, file.path(cfg$out_dir, "loso_metrics_pooled.csv"), row.names = FALSE)
write.csv(by_site, file.path(cfg$out_dir, "loso_metrics_by_site.csv"), row.names = FALSE)
write.csv(sse_share, file.path(cfg$out_dir, "loso_sse_share.csv"), row.names = FALSE)
write.csv(coverage, file.path(cfg$out_dir, "loso_interval_coverage.csv"), row.names = FALSE)
write.csv(do.call(rbind, stack_weights), file.path(cfg$out_dir, "stack_weights_by_heldout_site.csv"))

cat("\nPooled held-out performance:\n"); print(pooled, digits = 3)
cat("\nRF ", 100 * cfg$pi_level, "% interval coverage:\n", sep = ""); print(coverage, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "sessionInfo.txt"))
