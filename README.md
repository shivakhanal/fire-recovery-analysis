# Post-Fire Forest Canopy Recovery - R Code

This repository contains the R code used to quantify post-fire canopy recovery
in south-east Australian forests using MODIS FPAR (fraction of absorbed
photosynthetically active radiation) time series (2001–2022).

## What the code does

The code takes MODIS FPAR time series for burnt and unburnt forest pixels and:

1. Groups pixels into forest phenology types based on their long-term FPAR
   patterns (the clustering itself was done separately in TerrSet - see note below)
2. For each burnt pixel, detects when FPAR dropped (fire impact) and how large
   the drop was relative to pre-fire conditions
3. Estimates how long it took for the burnt pixel's FPAR to recover to levels
   similar to nearby unburnt forest of the same phenology type
4. Produces summary statistics and time series plots for inspection

Recovery is assessed only during the annual peak-growth period of each
phenology type, to avoid false positives from off-season fluctuations. Two
recovery criteria are applied in sequence: the burnt pixel must first reach
95% of the current-year unburnt level (Stage 1); if that is not achieved
within the study period, it is checked against 90% of the long-term average
peak (Stage 2). Pixels meeting neither criterion by 2022 are recorded as
not recovered.

> **Note on clustering:** The S-mode PCA and k-means clustering (k = 50) that
> assign each pixel to a phenology type were carried out in the Earth Trends
> Modeler (ETM) in TerrSet software. These steps are not reproduced in R. The
> cluster assignment file is a required input to this pipeline.

---

## Files

```
fire-recovery-analysis/
├── 01_data_preprocessing.R          # Load and check input data
├── 02_fire_impact_detection.R       # Measure fire impact for each burnt pixel
├── 03_recovery_estimation.R         # Estimate time to canopy recovery
├── 04_visualization.R               # Time series plots and summary figures
├── utils/
│   ├── helper_functions.R           # Shared utilities (date conversion, logging)
│   ├── reference_value_calculation.R # Identify peak-growth periods
│   └── cluster_homogeneity_analysis.R # Check clustering quality

```

---

**Fire year convention:** Australian fire years run from 1 July to 30 June.
The code uses this convention throughout.

### R packages needed

```r
install.packages(c("terra", "sf", "dplyr", "ggplot2", "data.table", "stringr"))
```

R version 4.0 or later is recommended. At least 16 GB RAM is advised when
working with the full raster stack.

---

## Output

For each burnt pixel the code returns:

- **Impact magnitude** - how much FPAR dropped relative to pre-fire levels
- **Recovery time** - days from the FPAR minimum to recovery (-1 if not
  recovered within 2001–2022)
- **Recovery method** - whether recovery was detected by Stage 1, Stage 2,
  or not at all

Results are saved as CSVs and diagnostic time series plots for each
cluster–fire-year combination.

---

## Contact
Shiva Khanal - 1khanalshiva@gmail.com



