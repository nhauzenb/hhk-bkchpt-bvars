Code package for Hauzenberger, N., F. Huber, & G. Koop (2023). Macroeconomic forecasting using BVARs. Chapter for *Handbook of Macroeconomic Forecasting*, edited by Mike Clements and Ana Galvão.

**Publication ([draft chapter](https://www.dropbox.com/scl/fi/cry8xuxkwwdtc3matz8g1/HHK_bookchp.pdf?rlkey=45ysy3b2hpqykkxormms9bipe&dl=0)).** 

### Data. 
For the forecast exercise, we use the popular FRED-QD dataset provided by the [*Federal Reserve Bank of St. Louis*](https://research.stlouisfed.org/econ/mccracken/fred-databases/). We provide the data as a .rda file [`fred QD`](./fred_data/fred_QD.rda) with 30 columns of the *"Xraw.stat"* object referring to the different variables. This quarterly sample spans from 1965Q1 to 2019Q4. We deliberately exclude the Covid-19 period and focus exclusively on pre-pandemic data. The data is transformed to stationarity, following the suggestions of McCracken and Ng (2020). Table 1 in Sub-section 4.1 of the book chapter shows the set of variables included for different model sizes.

### Bayesian estimation. 
For each model type, we specify an estimation grid over the full evaluation sample and the information sets. In terms of model specifications, we consider four conjugate and four non-conjugate VAR priors. Detailed information on how we estimate these different models is provided below. There are three main estimation/forecasting files: 

* [`!fcst conjVAR`](!fcst_conjVAR.R) for conjugate BVARs, 
* [`!fcst subVAR`](!!fcst_subVAR.R) for conjugate BVARs with subspace shrinkage and
* [`!fcst nconjVAR`](!fcst_nconjVAR.R) for non-conjugate BVARs.

1.) [`Conjugate BVARs`](!fcst_conjVAR.R): The hyperparameters can be treated as unknown and sampled based on an additional inverse transform sampling step, or a-priori chosen with an optimization algorithm (which maximizes the marginal likelihood). Three of the four conjugate variants can be estimated with this function:   

  * A conjugate VAR with a non-informative prior (our benchmark) with the option *model == "conjVAR-FLAT"*.
  * A classic symmetric conjugate Minnesota prior with the options *model == "conjVAR-MINg" (grid for hyperparameters) or *model == “conjVAR-MINo"* (optimization algorithm for hyperparameters).
  * An asymmetric conjugate Minnesota prior with the options *model == "conjVAR-ASYMg" (grid for hyperparameters) or model == "conjVAR-ASYMo" (optimization algorithm for hyperparameters).

2.) [`A conjugate BVAR with a subspace shrinkage prior`](!!fcst_subVAR.R): This is the fourth variant of the conjugate models. The hyperparameters related to the subspace shrinkage prior are treated as unknown and sampled based on an additional inverse transform sampling step.

3.) [`Non-conjugate BVAR`](!fcst_nconjVAR.R): All four non-conjugate variants can be estimated with this function. We use different variants of the global-local shrinkage priors. We consider an adaptive Minnesota prior (*model == "nconjVAR-MINh"*), a Horseshoe prior (*model == "nconjVAR-HS"*), a LASSO prior (*model == "nconjVAR-lasso"*), and a Normal-Gamma prior (*model == "nconjVAR-NG"*). All non-conjugate VARs are estimated using the Gibbs sampling algorithm proposed in [Chan et al. (2022, JoE)](https://doi.org/10.1016/j.jeconom.2021.11.010). Models can be estimated either with homoskedastic errors (*sv == "homo"*) or with stochastic volatility (*sv == "SV").

The folder [`bvar funcs`](./bvar_funcs/) contains the MCMC samplers for each conjugate and non-conjugate variant: 

* [`Direct sampler for conjugate BVAR with classic symmetric Minnesota prior`](./bvar_funcs/conjVARstd_func.R) 
* [`Direct sampler for BVAR with an asymmetric conjugate Minnesota prior`](./bvar_funcs/conjVARasym_func.R)
* [`Direct sampler for conjugate BVAR with a subspace shrinkage prior`](./bvar_funcs/conjVARsub_func.R)
* [`Gibbs sampler for non-conjugate BVAR with global-local shrinkage priors`](./bvar_funcs/nconjVAR_func.R)


Replication codes come without technical support of any kind.
