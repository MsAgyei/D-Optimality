# D-Optimality
Code for the Paper on D  Optimality
# Two-Point Visit Design for Estimating Grapevine Flowering Curves

Code accompanying the paper *"[Paper Title Here]"* (submitted to *Plant Methods*).

## Overview

This repository contains the simulation code used to compare **D-optimal** and 
**C-optimal** two-point visit designs for estimating the Day50% flowering point 
(the day at which 50% of flowers are open) in grapevines, using a logistic 
flowering curve model.

Two simulation studies are included:

- **Simulation 1** — a simulation study under a known true model (α = 0, β = 1), 
  used to validate the D-/C-optimal design comparison across a range of 
  inflorescence sample sizes (N_bunches = 10, 20, 50).
- **Simulation 2** — a simulation study using fitted flowering model parameters 
  from three vineyard sites, evaluated across a wider range of sample sizes 
  (N_bunches = 10–100). Includes a reference comparison against the RMSE 
  achieved by the sites' actual visit schedules.

## Contents

| File | Description |
|---|---|
| `flowering_design_simulations.R` | All simulation code: Simulation 1, Simulation 2 (setup and main loop), and the actual-visit RMSE reference comparison. |

## Requirements

- R (≥ 4.0 recommended)
- Package: `dplyr`

No external data files are required. Site-level model parameters (posterior 
means from the fitted Bayesian flowering model reported in the paper) are 
provided directly in the script.

## Usage

Run `flowering_design_simulations.R` in full, from top to bottom, in R. The 
script prints summary tables to the console as it runs and stores results in 
memory as named lists (`sim1_heat_grid_results`, `sim2_heat_grid_results`, 
`actual_rmse_real`) for further inspection.

## Citation

If you use this code, please cite:

> [Author names]. [Year]. [Paper title]. *Plant Methods*. [DOI/link]

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
