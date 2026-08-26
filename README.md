# POTL
Prediction-Oriented Transfer Learning for Survival Analysis

The `simulations` and `application` folders contain the code for the simulation studies and the real-data application, respectively.

All numerical studies were conducted on a CPU server using parallel computing.


## Simulations 

For single-source simulations with a common set of covariates in the target and source studies, run:

    Rscript sim_single_source.R
    Rscript plot_single_source_figures.R single_source

For SC4' and SC5' where POTL uses a misspecified Cox source model, run:
    Rscript sim_single_source_cox_sc45.R

For source-distribution shift, set `shift_s = TRUE` and `shift_v = FALSE` in `sim_single_source.R`, then run:

    Rscript sim_single_source.R
    Rscript plot_single_source_figures.R single_source_source_shift . "simres_sc%d.RData" "individual_survival_bias_single_source_source_shift_sc%d.RData"

For validation-distribution shift, set `shift_s = FALSE` and `shift_v = TRUE` in `sim_single_source.R`, then run:

    Rscript sim_single_source.R
    Rscript plot_single_source_figures.R single_source_validation_shift . "simres_sc%d.RData" "individual_survival_bias_single_source_validation_shift_sc%d.RData"

For single-source simulations with different sets of covariates in the target and source studies, run:

    Rscript sim_single_source_diffcovar.R
    Rscript plot_single_source_figures.R single_source_diffcovar

For multi-source simulations, first set `pool = "cv"` in `sim_multi_source.R` and run:

    Rscript sim_multi_source.R

Then set `pool = "ns"` in `sim_multi_source.R` and run the same command again:

    Rscript sim_multi_source.R

After both runs have finished, create the merged figures with:

    Rscript plot_multi_source_figures.R


## Application

To perform the main analysis and estimate the target model, run:

    Rscript analysis.R

To generate predictions for future patients from the target population, run:

    Rscript prediction.R
