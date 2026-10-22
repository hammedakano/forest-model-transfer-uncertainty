# forest-model-transfer-uncertainty

R tools for testing whether models of tree attributes (DBH, height, aboveground biomass) work on sites they were not trained on, and how much to trust their predictions.

**Skills shown:** machine learning and statistical modelling of forest attributes, transferability across sites, uncertainty quantification, agreement analysis between field and remote sensing measurements.

## Why this matters

Random train and test splits put trees from the same site into both sets, so they overstate accuracy. An operational inventory model has to work on a new stand, species, sensor or region. Every model here is therefore scored with leave-one-site-out (LOSO) cross-validation.

In my six-site PhD data, random splits gave high R2. Under LOSO, pooled R2 for DBH was negative, and most of the error came from a small number of structurally different sites.

## Contents

| File | Purpose |
|---|---|
| `R/loso_transfer_models.R` | Compares a linear model, random forest (ranger), XGBoost and a convex stack under LOSO. Stack weights are fitted only on training sites. |
| `R/utils/metrics.R` | Lin's CCC, out-of-sample R2, RMSE, relative RMSE, bias and proportional Bland-Altman limits. |
| `R/validation/agreement_mixed_model.R` | Field vs UAV agreement with a mixed model (`field ~ uav + (1 | site)`), marginal and conditional R2, and repeated-measures correlation. |
| `R/validation/ccc_logscale_proportional_BA.R` | CCC on raw and log scales, a heteroscedasticity test and proportional Bland-Altman analysis. |
| `examples/simulate_example_data.R` | Synthetic six-site data with site effects, so the code runs without private data. |

## Outputs of the LOSO script

- pooled and per-site held-out r, R2, RMSE, relative RMSE, bias and CCC
- each site's share of total squared error
- coverage and width of random forest 90% prediction intervals for each held-out site
- stack weights for each held-out site

## Quick start

```bash
Rscript examples/simulate_example_data.R
Rscript R/loso_transfer_models.R
```
Results are written to `outputs/loso/`. On the synthetic data it runs in about 10 seconds.

To use your own data:
```bash
LOSO_DATA=data/trees.csv LOSO_RESPONSE=agb LOSO_PREDICTORS=Z,crown_area_m2 LOSO_SITE=site Rscript R/loso_transfer_models.R
```
The validation scripts read a plot-level table set with `AGREEMENT_DATA`. It needs columns `height`, `Z`, `dbh`, `est_DBH`, `agb`, `est_agb`, `site` and `plot`.

## Installation

```r
source("install_packages.R")   # dplyr, ranger, xgboost, lme4, lmerTest, MuMIn, rmcorr, epiR, ggplot2, patchwork
```

## Data

Field and UAV data are not included. `.gitignore` excludes `data/` and `outputs/`.

## Related repositories

- [tls-stem-reconstruction](../tls-stem-reconstruction): TLS stem segmentation and diameter profiles
- [uav-sfm-tree-metrics](../uav-sfm-tree-metrics): drone point cloud tree detection and crown metrics
- [forest-model-transfer-uncertainty](../forest-model-transfer-uncertainty): site-transfer validation and prediction intervals

## Author

Hammed Akano, PhD candidate, Deakin University. Contact hammedaakano@gmail.com.

## Licence

MIT. See `LICENSE`.
