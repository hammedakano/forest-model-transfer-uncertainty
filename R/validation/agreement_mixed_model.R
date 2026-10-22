# Field vs UAV agreement, accounting for site clustering
# Replaces the pooled Pearson r / R-squared for height, DBH and AGB
#
# Variable names follow 3D_allData.Rmd:
#   field height = height,  UAV height = Z
#   field DBH    = dbh,     UAV DBH    = est_DBH
#   field AGB    = agb,     UAV AGB    = est_agb
#   grouping     = site,    plot id    = plot
#
# Reported validation uses plot-level medians (n = 65), so we aggregate first.

library(dplyr)
library(lme4)
library(lmerTest)   # gives p value for the fixed-effect slope
library(MuMIn)      # marginal and conditional R-squared
library(rmcorr)

# ---- 1. Aggregate tree-level rows to plot-level medians (unit of analysis) ----
# If you already work from a plot-median data frame, skip this and rename it plot_med.

all_data <- read.csv(Sys.getenv("AGREEMENT_DATA", unset = "data/plot_level_field_uav.csv"))  # 
all_data = all_data[all_data$est_agb != -216.77, ] 
all_data = all_data[all_data$est_agb != 3348.35, ] # remove outliers in est_agb

plot_med <- all_data %>%
  group_by(site, plot) %>%
  summarise(
    field_height = median(height,  na.rm = TRUE),
    uav_height   = median(Z,       na.rm = TRUE),
    field_dbh    = median(dbh,     na.rm = TRUE),
    uav_dbh      = median(est_DBH, na.rm = TRUE),
    field_agb    = median(agb,     na.rm = TRUE),
    uav_agb      = median(est_agb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site = factor(site))

# ---- 2. Mixed-effects agreement model: field ~ uav + (1 | site) ----
# The fixed-effect slope tells you whether field and UAV measurements track each
# other after site-level clustering is accounted for. Marginal R-squared (R2m) is
# the variance explained by the UAV term alone; conditional R-squared (R2c) adds site.
fit_agreement <- function(field, uav, data) {
  d <- data.frame(field = data[[field]], uav = data[[uav]], site = data$site)
  d <- d[stats::complete.cases(d), ]
  m <- lmer(field ~ uav + (1 | site), data = d, REML = TRUE)
  r2 <- MuMIn::r.squaredGLMM(m)
  list(
    n        = nrow(d),
    slope    = unname(fixef(m)["uav"]),
    slope_ci = confint(m, parm = "uav", method = "Wald"),
    p_value  = summary(m)$coefficients["uav", "Pr(>|t|)"],
    R2_marginal    = unname(r2[1, "R2m"]),
    R2_conditional = unname(r2[1, "R2c"]),
    site_sd  = as.data.frame(VarCorr(m))$sdcor[1],
    model    = m
  )
}

height_mix <- fit_agreement("field_height", "uav_height", plot_med)
dbh_mix    <- fit_agreement("field_dbh",    "uav_dbh",    plot_med)
agb_mix    <- fit_agreement("field_agb",    "uav_agb",    plot_med)

height_mix$model
dbh_mix$model
agb_mix$model

# Quick view
lapply(list(height = height_mix, dbh = dbh_mix, agb = agb_mix), function(x) {
  data.frame(n = x$n, slope = round(x$slope, 3),
             ci_low = round(x$slope_ci[1], 3), ci_high = round(x$slope_ci[2], 3),
             p = signif(x$p_value, 3),
             R2m = round(x$R2_marginal, 3), R2c = round(x$R2_conditional, 3),
             site_sd = round(x$site_sd, 3))
}) %>% bind_rows(.id = "attribute")

# ---- 3. Repeated-measures correlation (within-site common slope) ----
# Alternative to a pooled Pearson r. Reports the within-site association.
height_rm <- rmcorr(participant = site, measure1 = field_height, measure2 = uav_height, dataset = plot_med)
dbh_rm    <- rmcorr(participant = site, measure1 = field_dbh,    measure2 = uav_dbh,    dataset = plot_med)
agb_rm    <- rmcorr(participant = site, measure1 = field_agb,    measure2 = uav_agb,    dataset = plot_med)

data.frame(
  attribute = c("height", "dbh", "agb"),
  rmcorr_r  = round(c(height_rm$r, dbh_rm$r, agb_rm$r), 3),
  df        = c(height_rm$df, dbh_rm$df, agb_rm$df),
  p         = signif(c(height_rm$p, dbh_rm$p, agb_rm$p), 3)
)

# Note on height: keep the plantation-excluded comparison as a sensitivity check,
# but run it the same way (filter plot_med to the four native sites, refit).


library(ggplot2)
library(patchwork)
library(dplyr)

# S1-S6 labels
plot_med <- plot_med %>%
  mutate(site = recode(as.character(site),
                       "Byaduk" = "S1", "White" = "S2", "Gelantipy" = "S3",
                       "Redbill" = "S4", "Prison" = "S5", "Otways" = "S6"),
         site = factor(site, levels = c("S1","S2","S3","S4","S5","S6")))

panel <- function(data, xvar, yvar, slope, intercept, lab, origin = FALSE) {
  p <- ggplot(data, aes(.data[[xvar]], .data[[yvar]], colour = site)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_abline(slope = slope, intercept = intercept, linewidth = 0.7) +
    geom_point(size = 2, alpha = 0.85) +
    labs(x = paste("UAV", lab), y = paste("Field", lab), colour = "Site") +
    theme_bw()
  if (origin) p <- p + expand_limits(x = 0, y = 0)   # include 0,0; clip the sub-zero sliver
  p
}

ph <- panel(plot_med, "uav_height", "field_height", 0.239,   9.739, "height (m)", origin = TRUE)               # no origin
pd <- panel(plot_med, "uav_dbh",    "field_dbh",    1.017,  -0.417, "DBH (cm)",  origin = TRUE)
pa <- panel(plot_med, "uav_agb",    "field_agb",    1.021, -15.728, "AGB (kg)",  origin = TRUE)

(ph | pd | pa) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")


####### plot B #######
panel <- function(data, xvar, yvar, slope, intercept, lab, origin = TRUE) {
  lims <- range(c(data[[xvar]], data[[yvar]], if (origin) 0), na.rm = TRUE)
  brks <- scales::breaks_pretty(5)(lims)
  
  ggplot(data, aes(.data[[xvar]], .data[[yvar]], colour = site)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_abline(slope = slope, intercept = intercept, linewidth = 0.7) +
    geom_point(size = 2, alpha = 0.85) +
    scale_x_continuous(breaks = brks) +
    scale_y_continuous(breaks = brks) +
    coord_equal(xlim = lims, ylim = lims, expand = TRUE) +
    labs(x = paste("UAV", lab), y = paste("Field", lab), colour = "Site") +
    theme_bw()
}

ph <- panel(plot_med, "uav_height", "field_height", 0.239,   9.739, "height (m)")
pd <- panel(plot_med, "uav_dbh",    "field_dbh",    1.017,  -0.417, "DBH (cm)")
pa <- panel(plot_med, "uav_agb",    "field_agb",    1.021, -15.728, "AGB (kg)")

(ph | pd | pa) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")
