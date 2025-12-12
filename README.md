# BC Wildfire Smoke Trends (2013-2023)

## Overview
This project analyzes the increasing impact of wildfire smoke in British Columbia (2013-2023). It integrates ground-level monitoring (NAPS), hazard mapping (HMS), and satellite data (MAIAC) to assess trends in smoke frequency vs. intensity.

## Key Findings
- **Frequency:** Smoke days are increasing ($\tau = 0.42$).
- **Intensity:** Peak intensity trends are moderate but variable ($\tau = 0.27$).
- **Spatial:** Smoke is spatially clustered (Moran's I = 0.36), and satellite detection accuracy degrades in some areas ($r = 0.55$).

## Repository Structure
- `data/`: Processed daily smoke metrics.
- `scripts/`: R code for statistical tests, map, and plotting.
- `outputs/`: Final visualizations used in the study.

## Methods
I used R for all statistical analyses, utilizing the `Kendall` package for trend analysis and `spdep` for spatial autocorrelation.

## Author
Rachmania
