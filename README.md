# Post-Fire Vegetation Recovery Analysis Using MODIS FPAR Time Series

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains R code for quantifying fire impact and post-fire vegetation recovery using MODIS FPAR (Fraction of Photosynthetically Active Radiation) time series data. The methodology uses k-means clustering to group pixels with similar temporal patterns and estimates recovery duration by comparing burnt pixel trajectories with unburnt reference pixels within the same cluster.

**Associated Publication:**  
[Your Name et al. (2025). "Quantifying Post-Fire Recovery of Forest Canopies..." Journal Name. DOI: XX.XXXX/XXXXX]

## Key Features

- **Fire Impact Detection**: Identifies maximum FPAR decline during fire periods
- **Reference Peak Identification**: Calculates modal peak timing across years for each cluster
- **Two-Stage Recovery Detection**:
  - **Stage 1 (Dynamic Catchup)**: Burnt vegetation reaches ≥95% of contemporary unburnt vegetation during peak periods
  - **Stage 2 (Static Threshold)**: Burnt vegetation reaches ≥90% of long-term reference mean during peak periods
- **Peak-Period Restriction**: Recovery detection limited to vegetation peak periods to reduce false positives
- **Automated Diagnostics**: Generates time series plots and quality metrics

## Requirements

### System Requirements
- R ≥ 4.0.0
- 16+ GB RAM recommended for large raster stacks
- Multi-core processor (parallel processing supported)

### R Packages
```r
# Core spatial packages
install.packages(c("terra", "sf", "raster", "sp"))

# Data manipulation and visualization
install.packages(c("dplyr", "tidyr", "ggplot2", "viridis", "gridExtra"))

# Time series analysis
install.packages("changepoint")

# Utilities
install.packages("stringr")
```

## Installation
```r
# Clone the repository
git clone https://github.com/yourusername/fire-recovery-analysis.git
cd fire-recovery-analysis

# Source helper functions
source("scripts/utils/helper_functions.R")
source("scripts/utils/reference_value_calculation.R")
```

## Quick Start

### 1. Prepare Input Data

Required inputs:
- **FPAR Time Series**: Multi-temporal MODIS FPAR raster stack (500m, 8-day resolution)
- **Cluster Map**: K-means cluster assignments for each pixel
- **Fire Polygons**: Spatial polygons with fire year attributes
- **Unburnt Reference**: Mean FPAR time series for unburnt pixels by cluster and fire year
```r
library(terra)
library(sf)

# Load data
fpar_stack <- rast("data/modis_fpar_stack.tif")
clusters <- read.csv("data/cluster_assignments.csv")
fire_polygons <- st_read("data/fire_history.gpkg")
unburnt_reference <- read.csv("data/unburnt_mean_timeseries.csv")
```

### 2. Run Fire Impact Detection
```r
source("scripts/02_fire_impact_detection.R")

# Detect fire impact
fire_impact <- detect_fire_impact(
  burnt_ts = burnt_pixel_timeseries,
  unburnt_ts = unburnt_reference_timeseries,
  fire_indices = list(start_idx = 230, end_idx = 276),
  pre_fire_baseline = 5
)
```

### 3. Calculate Recovery Time
```r
source("scripts/03_recovery_estimation.R")

# Calculate recovery with strict peak-period detection
recovery_result <- calculate_recovery_time(
  burnt_ts = burnt_pixel_timeseries,
  unburnt_ts = unburnt_reference_timeseries,
  fire_impact_result = fire_impact,
  recovery_threshold = 0.95,  # Stage 1: 95% of unburnt
  n_per_year = 46,
  strict_peak_only = TRUE
)

# Access results
recovery_result$recovery_days  # Time to recovery in days
recovery_result$recovery_method  # "peak_period_match" or "no_peak_recovery"
recovery_result$reference_mean  # Long-term reference value
```

### 4. Generate Visualizations
```r
source("scripts/04_visualization.R")

# Create diagnostic plots
create_diagnostic_plots(
  results_df = recovery_results,
  cluster_id = 23,
  fire_year = 200607,
  extracted_values = burnt_pixel_data,
  unburnt_ts = unburnt_reference,
  column_names = colnames(burnt_pixel_data),
  sample_size = 3
)
```

## Methodology

### Fire Impact Detection

Fire impact is quantified as the maximum FPAR decline during the fire period:
```
Impact Magnitude = Pre-fire Baseline - Minimum FPAR during Fire Period
```

Where:
- **Pre-fire Baseline**: Mean FPAR of 5 observations before fire start
- **Fire Period**: June 1st (fire year) to July 31st (following year)
- **Impact Date**: Temporal index when minimum FPAR occurred

### Reference Peak Identification

For each cluster, reference peaks are identified using:

1. **Within-Year Peak Detection**: Find 5-observation window with highest mean FPAR each year
2. **Modal Peak Position**: Identify most common peak position (1-46) across all years
3. **Reference Values**: Collect all FPAR values from modal position across years

This approach uses the entire time series (2000-2022) rather than only pre-fire data, improving statistical robustness.

### Two-Stage Recovery Detection

#### Stage 1: Dynamic Catchup (Primary Criterion)
```
Burnt FPAR ≥ 0.95 × Contemporary Unburnt FPAR (at peak periods)
```

**Advantages**:
- Accounts for inter-annual vegetation variability
- Represents true ecological recovery to ambient conditions
- Adaptive to climate and seasonal variations

#### Stage 2: Static Threshold (Secondary Criterion)
```
Burnt FPAR ≥ 0.90 × Long-term Reference Mean (at peak periods)
```

**Purpose**:
- Captures substantial but incomplete recovery
- Provides consistent benchmark across all pixels
- Identifies meaningful vegetation restoration

#### Peak-Period Restriction

Recovery detection is **strictly limited to peak periods** to ensure:
- Maximum vegetation expression comparison
- Reduced false positives from off-season fluctuations
- Temporal consistency across years

### Recovery Time Calculation
```
Recovery Time (days) = (Recovery Index - Impact Index) × 8 days
```

Pixels that don't recover within the observation period are assigned -1.

## Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `recovery_threshold` | 0.95 | Stage 1 threshold (proportion of unburnt FPAR) |
| `static_threshold` | 0.90 | Stage 2 threshold (proportion of reference mean) |
| `n_per_year` | 46 | Number of 8-day observations per year |
| `window_size` | 5 | Size of peak detection window (observations) |
| `pre_fire_baseline` | 5 | Number of pre-fire observations for baseline |
| `strict_peak_only` | TRUE | Limit recovery detection to peak periods only |

**Adjusting Thresholds**: Lower thresholds (0.80-0.85) detect partial recovery; higher thresholds (0.95-1.00) require complete recovery.

## Example Output
```r
# Recovery result structure
$recovery_days
[1] 584  # Days from impact to recovery

$recovery_idx  
[1] 303  # Temporal index of recovery

$recovered
[1] TRUE

$recovery_method
[1] "peak_period_match"  # or "no_peak_recovery"

$reference_mean
[1] 0.723  # Long-term reference FPAR value

$modal_peak_position
[1] 19  # Peak typically occurs at observation 19 (mid-year)
```

## Validation

The methodology has been validated against:
- Independent fire severity maps (Kruskal-Wallis H = 193.85, p < 0.01)
- Field-based recovery assessments
- Comparison with prescribed vs. wildfire recovery patterns

## Troubleshooting

**Issue**: No recovery detected for pixels that visually appear recovered

**Solution**: Check `recovery_threshold` parameter. Consider lowering to 0.90 or 0.85 for ecosystems with incomplete canopy recovery.

---

**Issue**: Many pixels show very rapid recovery (<100 days)

**Solution**: This may indicate:
1. Low-severity fire with minimal canopy damage
2. Patchiness within fire polygons
3. False positives - increase `recovery_threshold`

---

**Issue**: Peak period detection fails

**Solution**: Verify `n_per_year = 46` matches your data temporal resolution. Check for data quality issues.

See [docs/troubleshooting.md](docs/troubleshooting.md) for complete guide.

## Citation
If you use this code in your research, please cite:

```bibtex
@article{yourname2025,
  title={Quantifying Post-Fire Recovery of Forest Canopies in South-East Australia Using MODIS FPAR Time Series and Clustering Analysis},
  author={Your Name and Co-Authors},
  journal={Journal Name},
  year={2025},
  doi={XX.XXXX/XXXXX}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Shiva Khanal
- **Email**: 1khanalshiva@gmail.com
- **Issues**: [GitHub Issues](https://github.com/shivakhanal/fire-recovery-analysis/issues)

