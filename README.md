# POTL
Prediction-Oriented Transfer Learning for Survival Analysis

The `simulation` and `application` folders contain the code for the simulation studies and the real-data application, respectively.

All numerical studies were conducted on a CPU server using parallel computing.


## Simulation 

For single-source simulations with a common set of covariates in the target and source studies, run:

    Rscript sim_single_source.R

To consider covariate shift between the target and source data, set the `shift_s` variable to `TRUE` in line 28 of `sim_single_source.R`.

To consider covariate shift between the training and validation data, set the `shift_v` variable to `TRUE` in line 30 of `sim_single_source.R`.

For single-source simulations with different sets of covariates in the target and source studies, run:

    Rscript sim_single_source_diffcovar.R

For multi-source simulations, run:

    Rscript sim_multi_source.R


## Application

To perform the main analysis and estimate the target model, run:

    Rscript analysis.R

To generate predictions for future patients from the target population, run:

    Rscript prediction.R
