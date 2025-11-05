# Fire Recovery Analysis Methodology

## Overview

This document provides detailed technical documentation of the fire impact and recovery analysis methodology.

## Data Requirements

### Input Data

1. **MODIS FPAR Time Series**
   - Format: Multi-band GeoTIFF raster stack
   - Resolution: 500m spatial, 8-day temporal
   - Coverage: 2000-2022 (or subset)
   - CRS: Any projected or geographic CRS
   - Quality: Pre-processed and gap-filled

2. **Cluster Assignments**
   - Format: CSV with columns: `class`, `x`, `y`
   - `class`: Cluster ID (integer)
   - `x`, `y`: Coordinates matching FPAR raster CRS
   - One row per pixel

3. **Fire Polygons**
   - Format: GeoPackage or Shapefile
   - Required column: `FireYear` (YYYYMM format, e.g., 200607)
   - Optional: `FireType` (planned/unplanned)
   - CRS: Any (will be transformed to match raster)

4. **Unburnt Reference Time Series**
   - Format: CSV with mean FPAR for unburnt pixels
   - Required columns: `cluster`, `fireyear`, time series columns
   - Time series columns named consistently (e.g., `X2000001`, `X2000002`, ...)

## Algorithm Details

### 1. Fire Impact Detection

**Objective**: Quantify the magnitude and timing of fire impact on vegetation FPAR.

**Steps**:

1. Define fire period boundaries:
   - Start: June 1st of fire year → `fire_start_idx`
   - End: July 31st of following year → `fire_end_idx`

2. Calculate pre-fire baseline:
```
   pre_fire_mean = mean(burnt_ts[(fire_start_idx - 5):(fire_start_idx - 1)])
```

3. Find minimum FPAR during fire period:
```
   impact_idx = which.min(burnt_ts[fire_start_idx:fire_end_idx])
   impact_value = burnt_ts[impact_idx]
```

4. Calculate impact magnitude:
```
   impact_magnitude = pre_fire_mean - impact_value
```

**Output**: 
- `impact_magnitude`: FPAR decline
- `impact_date_idx`: Temporal index of impact
- `pre_fire_mean`: Pre-fire baseline FPAR
- `impact_value`: Minimum FPAR value

### 2. Reference Peak Identification

**Objective**: Identify characteristic peak FPAR timing and values for unburnt vegetation.

**Steps**:

1. For each complete year in time series:
```
   for year in 1:n_years:
       year_data = unburnt_ts[(year-1)*46+1 : year*46]
       
       # Find 5-observation window with highest mean
       for position in 1:(46-5+1):
           window = year_data[position:(position+4)]
           if mean(window) > best_mean:
               best_position = position
               best_mean = mean(window)
```

2. Identify modal peak position:
```
   modal_position = mode(yearly_peak_positions)
```

3. Collect reference values from modal position across all years:
```
   for year in 1:n_years:
       year_start = (year-1) * 46 + 1
       indices = year_start + modal_position + (0:4)
       reference_values.append(unburnt_ts[indices])
```

4. Calculate reference statistics:
```
   reference_mean = mean(reference_values)
   reference_sd = sd(reference_values)
```

**Output**:
- `reference_values`: All FPAR values from modal peak periods
- `modal_peak_position`: Most common peak position (1-46)
- `reference_mean`: Mean peak FPAR
- `yearly_peak_indices`: All peak indices across years

### 3. Two-Stage Recovery Detection

**Objective**: Detect when burnt vegetation recovers to unburnt condition.

#### Stage 1: Dynamic Catchup (Primary)

**Criterion**: Burnt FPAR ≥ 95% of contemporary unburnt FPAR during peak periods

**Algorithm**:
```
# Get all future peak periods after fire
future_peaks = yearly_peak_indices[yearly_peak_indices > impact_idx]

for peak_idx in future_peaks:
    burnt_value = burnt_ts[peak_idx]
    unburnt_value = unburnt_ts[peak_idx]
    catchup_threshold = 0.95 * unburnt_value
    
    if burnt_value >= catchup_threshold:
        recovery_idx = peak_idx
        recovery_method = "stage1_dynamic_catchup"
        BREAK
```

**Advantages**:
- Accounts for inter-annual climate variability
- True ecological recovery to ambient conditions
- Adaptive to seasonal variations

#### Stage 2: Static Threshold (Secondary)

**Criterion**: Burnt FPAR ≥ 90% of long-term reference mean during peak periods

**Algorithm**:
```
static_threshold = 0.90 * reference_mean

for peak_idx in future_peaks:
    burnt_value = burnt_ts[peak_idx]
    
    if burnt_value >= static_threshold:
        recovery_idx = peak_idx
        recovery_method = "stage2_static_threshold"
        BREAK
```

**Purpose**:
- Captures substantial but incomplete recovery
- Consistent benchmark across all pixels
- Backup criterion if Stage 1 fails

#### Recovery Time Calculation
```
if recovery_detected:
    recovery_days = (recovery_idx - impact_idx) * 8
else:
    recovery_days = -1  # Not recovered
```

### 4. Peak-Period Restriction

**Critical Feature**: Recovery detection is **strictly limited** to peak periods.

**Rationale**:
- Ensures comparison at maximum vegetation expression
- Eliminates false positives from off-season fluctuations
- Maintains temporal consistency across years

**Implementation**:
- Only check `burnt_ts` values at indices in `yearly_peak_indices`
- Never check recovery between peak periods
- This ensures ecological validity of recovery detection

## Parameter Selection

### Key Parameters

| Parameter | Default | Range | Impact |
|-----------|---------|-------|--------|
| `recovery_threshold_stage1` | 0.95 | 0.80-1.00 | Higher = stricter recovery criterion |
| `recovery_threshold_stage2` | 0.90 | 0.75-0.95 | Higher = requires more complete recovery |
| `window_size` | 5 | 3-7 | Larger = more robust but less precise |
| `pre_fire_baseline` | 5 | 3-7 | More obs = more stable baseline |
| `n_per_year` | 46 | - | Must match data temporal resolution |

### Threshold Selection Guidelines

**For complete recovery detection**: Use 0.95-1.00 for Stage 1

**For partial recovery detection**: Use 0.85-0.90 for Stage 1

**For conservative estimates**: Use Stage 1 only (`strict_peak_only = TRUE`)

**For lenient estimates**: Enable both stages, lower thresholds

## Quality Control

### Diagnostic Checks

1. **Visual inspection** of time series plots:
   - Check fire impact is correctly identified
   - Verify recovery detection makes ecological sense
   - Confirm peak periods align with vegetation phenology

2. **Statistical validation**:
   - Compare with independent fire severity data
   - Check recovery rates by cluster/fire type
   - Assess relationship between impact and recovery time

3. **Data quality filters**:
   - Remove pixels with >20% missing FPAR values
   - Exclude pixels with impact_magnitude < 0 (data errors)
   - Flag pixels with recovery_days < 50 (suspiciously fast)

## Limitations and Considerations

1. **Temporal Resolution**: 8-day resolution may miss rapid changes immediately after fire

2. **Spatial Resolution**: 500m pixels may contain mixed fire severities

3. **Recovery Definition**: Functional recovery (FPAR) may not equal structural recovery

4. **Climate Variability**: Drought can delay recovery independent of fire effects

5. **Other Disturbances**: Insect outbreaks, subsequent fires can confound recovery

## References

See main manuscript for full methodological justification and validation.