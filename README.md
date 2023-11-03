Code package for Hauzenberger, N., F. Huber, & G. Koop (2023). Macroeconomic forecasting using BVARs. Chapter for *Handbook of Macroeconomic Forecasting*, edited by Mike Clements and Ana Galvão.

**Publication ([draft chapter](https://www.dropbox.com/scl/fi/cry8xuxkwwdtc3matz8g1/HHK_bookchp.pdf?rlkey=45ysy3b2hpqykkxormms9bipe&dl=0)).** 

### Data. 
For the forecast exercise, we use the popular FRED-QD dataset provided by the [*Federal Reserve Bank of St. Louis*](https://research.stlouisfed.org/econ/mccracken/fred-databases/). We provide the data as a .rda file (*fred_data/fred_QD.rda*) with 30 columns of the *"Xraw.stat"* object referring to the different variables. Our quarterly sample spans from 1965Q1 to 2019Q4. We deliberately exclude the Covid-19 period and focus exclusively on pre-pandemic data. The data is transformed to stationarity, following the suggestions of McCracken and Ng (2020). Table 1 in Sub-section 4.1 of the book chapter shows the set of variables included for different model sizes.

### Data. 
For each model type, we specify an estimation grid over the full evaluation sample and the information sets. In terms of model specifications, we consider four conjugate and four non-conjugate VAR priors. Detailed information on how we estimate these different models is provided below. There are three main estimation/forecasting files: 

* One for conjugate BVARs: [!fcst_conjVAR.R](!fcst_conjVAR.R). 
* One for conjugate BVARs with subspace shrinkage (!fcst_subVAR.R), and
* One for non-conjugate BVARs (!fcst_nconjVAR.R).



