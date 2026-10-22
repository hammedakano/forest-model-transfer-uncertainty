# Synthetic tree-level data for running the modelling scripts without
# access to the field data. Six "sites" differ in allometry and in UAV
# measurement bias, so transfer between sites is deliberately imperfect.
# These are NOT real measurements.
set.seed(42)
sites <- data.frame(
  site       = c("S1", "S2", "S3", "S4", "S5", "S6"),
  type       = c("Plantation", "Plantation", "Native", "Native", "Native", "Native"),
  n          = c(60, 55, 40, 25, 35, 45),
  h_mean     = c(18, 20, 24, 15, 22, 28),
  dbh_slope  = c(1.25, 1.30, 1.55, 1.40, 1.50, 1.70),
  uav_h_bias = c(-0.5, -0.8, -1.5, -0.3, -1.0, -2.0)
)
rows <- lapply(seq_len(nrow(sites)), function(i) {
  s <- sites[i, ]
  height <- pmax(5, rnorm(s$n, s$h_mean, 4))
  dbh <- pmax(5, s$dbh_slope * height + rnorm(s$n, 0, 3))                 # cm
  agb <- exp(-2.0 + 2.4 * log(dbh) + rnorm(s$n, 0, 0.15))                 # kg, generic power law
  crown <- pmax(2, 0.35 * dbh + rnorm(s$n, 0, 3))                          # m2
  data.frame(site = s$site, type = s$type, tree = seq_len(s$n),
             height = round(height, 2), dbh = round(dbh, 1), agb = round(agb, 1),
             Z = round(height + s$uav_h_bias + rnorm(s$n, 0, 1), 2),
             crown_area_m2 = round(crown, 1))
})
out <- do.call(rbind, rows)
write.csv(out, file.path("examples", "example_tree_data.csv"), row.names = FALSE)
message("Wrote examples/example_tree_data.csv (", nrow(out), " synthetic trees)")
