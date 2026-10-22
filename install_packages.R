pkgs <- c("dplyr", "ranger", "xgboost", "lme4", "lmerTest", "MuMIn", "rmcorr",
          "epiR", "ggplot2", "patchwork")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
