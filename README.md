# Extending Joint Models with Nonlinear Mixed-Effects and Interval-Censored Dropout

## Overview

This repository implements an extended joint modeling framework for longitudinal and survival data, building upon the [HHJMs framework](https://github.com/oliviayu/HHJMs).

The proposed approach introduces two key methodological contributions:

- A biologically motivated nonlinear mixed-effects model for longitudinal antibody dynamics  
- A unified joint modeling framework supporting interval-censored survival outcomes  

These extensions enable more realistic modeling of immune response trajectories and informative dropout mechanisms, while providing a flexible framework for complex longitudinal and event time processes.

An example application to HIV vaccine data is provided to illustrate the implementation.

## Key Features

### Methodological Extensions

- Extension of the original HHJMs framework to support nonlinear mixed-effects longitudinal models through an exponential-delay formulation  
- Extension of the survival component to accommodate interval-censored Weibull proportional hazards models  
- Generalized likelihood construction for integrating standard GLME components with custom nonlinear longitudinal submodels  

### Framework and Utilities

- The framework provides joint modeling of longitudinal and survival processes under an h-likelihood formulation  
- It includes a unified estimation interface supporting both h-likelihood and adaptive Gauss-Hermite (aGH) methods  
- Multiple longitudinal submodels are supported, including LME, GLMM, and custom NLME formulations  
- Standard error estimation is available via both Hessian-based and aGH-based approaches  
- Example workflows are provided for HIV vaccine data analysis, model comparison, and simulation studies

## Model Framework
The framework jointly models longitudinal outcomes and survival data through shared random effects. The overall structure of the model is illustrated below. 
![Model](https://github.com/zhili617/joint-model-exponential-delay-interval-censoring/blob/main/orin/HHJMs-master/plot/framework.png)

## Repository Structure

## 📁 R
``` 
R/
├── JMfit.R
│ - Wrapper function for fitting the joint model.
│ - Dispatches estimation to either h-likelihood or adaptive Gauss–Hermite (aGH) methods.
│
├── JMfit_HL.R
│ - Core h-likelihood estimation routine for joint models combining
│ longitudinal and survival submodels.
│ - Modified: extended to support interval-censored Weibull PH survival models.
│
├── JMfit_aGH.R
│ - Adaptive Gauss–Hermite (aGH) implementation for joint model estimation.
│ - Uses quadrature-based integration over random effects.
│
├── JMsd_aGH.R
│ - Computes standard errors under the adaptive Gauss–Hermite framework.
│
├── JMsummary.R
│ - Generates summary tables of parameter estimates, standard errors,
│ Z-values, and p-values.
│
├── Jglme_loglike.R
│ - Constructs the joint log-likelihood for multiple longitudinal submodels.
│ - Modified: generalized likelihood construction via model_loglike(),
│ enabling integration of nonlinear mixed-effects (NLME) models.
│
├── model_loglike.R
│ - Unified likelihood interface for longitudinal models.
│ - New: dispatches between standard GLME and custom NLME likelihoods.
│
├── nlme_loglike.R
│ - New: defines the nonlinear mixed-effects (NLME) exponential-delay
│ log-likelihood for longitudinal processes.
│
├── weibull_interval_censored_like.R
│ - New: implements the Weibull proportional hazards log-likelihood
│ for interval-censored survival outcomes.
│
├── cox_loglike.R
│ - Constructs survival log-likelihood for standard Cox PH and Weibull models.
│
├── estBaseHazard.R
│ - Estimates the baseline hazard function for the Cox model.
│ - Modified: aligns time discretization with month-based event times.
│
├── estDisp.R
│ - Estimates dispersion parameters and covariance structures.
│ - Modified: improved initialization, parameter bounds, and added
│ Hessian-validity checks to enhance numerical stability during optimization.
│
├── estFixeff.R
│ - Estimates fixed effects via profile h-likelihood optimization.
│
├── estFixeff_adptGH.R
│ - Estimates model parameters under the adaptive Gauss–Hermite framework.
│
├── estRaneff.R
│ - Estimates subject-specific random effects via h-likelihood optimization.
│
├── get_raneff.R
│ - Wrapper for estimating and scaling random effects across subjects.
│
├── get_sd.R
│ - Computes standard errors using Hessian-based covariance estimation.
│
├── get_aGH_sd2.R
│ - Computes standard errors under the aGH framework via numerical derivatives.
│
├── evalMat.R
│ - Evaluates symbolic matrix expressions after substiting parameter values.
│
├── Vassign.R
│ - Assigns values to named parameter vectors.
│
├── Vderiv.R
│ - Computes symbolic derivatives of likelihood expressions.
│
├── getHessian.R
│ - Computes symbolic Hessian matrices.
│
├── getHmat.R
│ - Constructs Hessian components used in h-likelihood calculations.
│
├── strMat.R
│ - Builds structured symbolic matrices for covariance parameterization.
│
├── fmReverse.R
│ - Parses model formulas into response and covariates.
│ - Modified: supports interval-censored survival responses of the form Surv(L, R).
│
└── mgaussHermite.R
- Generates multivariate Gauss–Hermite quadrature points and weights.
```

### Modifications from the Original HHJMs Framework

This repository is based on the original HHJMs framework, with several methodological and implementation extensions.

#### 1. Nonlinear longitudinal modeling (NLME)

- `model_loglike.R`  
  Introduced a unified likelihood interface to support flexible integration of different longitudinal model types.

- `nlme_loglike.R`  
  Implemented a nonlinear mixed-effects exponential-delay model to capture biologically realistic longitudinal dynamics.

- `Jglme_loglike.R`  
  Generalized likelihood construction by replacing the original GLME-specific interface with `model_loglike()`.

---

#### 2. Interval-censored survival modeling

- `weibull_interval_censored_like.R`  
  Implemented a Weibull proportional hazards log-likelihood for interval-censored survival outcomes.

- `JMfit_HL.R`  
  Extended the h-likelihood estimation routine to support interval-censored Weibull survival models.

- `fmReverse.R`  
  Extended formula parsing to support survival responses of the form `Surv(L, R)`.

---

#### 3. Numerical and implementation refinements

- `estBaseHazard.R`  
  Modified baseline hazard estimation to better align with month-based time scales.

- `estDisp.R`  
  Improved parameter initialization and constraints, and added safeguards for invalid Hessian determinants to enhance optimization stability.

---
## 📁 example
```
example/
├── Longdata.csv / Survdata.csv
│ - Simulated datasets from the original HHJMs example.
│
├── test.R
│ - Main analysis script for HIV vaccine data.
│ - Performs data preprocessing and constructs:
│     * longitudinal dataset (long_data.csv)
│     * survival dataset (surv_data.csv)
│ - Implements joint models under:
│     * right-censored and interval-censored survival settings
│     * Cox and Weibull models
│ - Includes both linear mixed-effects models and nonlinear mixed-effects
│   (exponential delay model) within the joint modeling framework.
│
├── nlme_code.R
│ - Fits nonlinear mixed-effects models for the longitudinal process.
│ - Implements:
│     * exponential delay model (primary model)
│     * power-law model (comparison)
│ - Provides initial parameter values for joint modeling.
│ - Includes visualization and two-step estimation procedures.
│
├── graph_compare.R
│ - Compares linear mixed-effects and nonlinear exponential delay models.
│ - Generates subject-specific fitted curves.
│ - Computes evaluation metrics:
│     * MSE / SSE (longitudinal fit)
│     * AUC (survival prediction)
│
├── success_simulation.R
│ - Simulation framework for the joint model with nonlinear longitudinal
│   component (exponential delay model).
│ - Generates data including:
│     * nonlinear longitudinal process
│     * GLMM (binary outcome)
│     * survival data (right- or interval-censored Weibull)
│ - Fits joint models using JMfit().
│ - Repeats simulations until 100 successful runs.
│ - Computes performance metrics:
│     * bias, MSE
│     * coverage probability
│     * empirical vs model-based standard errors
│ - Stores outputs in:
│     example/simulationed_data_group/
│
└── simulationed_data_group/
   - Stores simulation outputs (.rds files).
   - Each file contains:
       * simulated longitudinal data
       * survival data
       * random effects and intermediate results
```

## 📁 man/

```
man/
├── *.Rd files
│ - Documentation files automatically generated by roxygen2.
│ - Not intended for manual editing.
```

---

## 📁 plot/

```
plot/
├── *.png / *.pdf / *.jpg
│ - Figures used in the project report/paper.
```

##  Original HHJMs Package
The following section is preserved from the original HHJMs package:
```
HHJMs-master/
├──.Rbuildignore
├──.gitignore
├──DESCRIPTION
├──HHJMs.Rproj
├──NAMESPACE
├──README.md
├──others.md
```


## Usage and Reproducible Analysis

This repository includes a complete reproducible analysis pipeline based on the proposed joint modeling framework.
A complete analysis pipeline is provided in the `example/` directory. The main workflow is implemented in `example/test.R`.

### Setup

Before running the example scripts, source all functions in the `R/` directory:

```r
srcpath <- "HHJMs-master/R"
file.sources <- list.files(srcpath, pattern = "*.r$", full.names = TRUE)
sapply(file.sources, source)
```

### Recommended workflow

To reproduce the main analysis, run the scripts in the following order:

1. **Fit nonlinear longitudinal models (NLME)**
 ```r
source("example/nlme_code.R")
```
   This step fits the exponential-delay model and provides initial parameter estimates for the joint model.

2.  **Run the joint model analysis**
   ```R
   source("example/test.R")
   ```

This script:
* prepares longitudinal and survival datasets

* constructs glmeObject and survObject

* fits joint models under:
  * right-censored and interval-censored settings
  * Cox and Weibull survival models

* compares model outputs

3. **(Optional) Model comparison and visualization**
 ```R
   source("example/graph_compare.R")
 ```
4. **(Optional) Simulation Study and Result**
```R
   source("example/success_simulation.R")
 ```

## Dependencies

The following R packages are required:

### Core dependencies
- `lme4`, `nlme`: mixed-effects modeling
- `survival`, `flexsurv`: survival modeling
- `glmmTMB`: generalized mixed models
- `MASS`, `Matrix`: numerical computation

### Visualization and evaluation
- `ggplot2`, `pROC`, `survivalROC`

### Simulation and utilities
- `truncnorm`, `expm`, `survsim`, `DescTools`, `tictoc`

Base R functions are used for additional data processing and numerical computation.

## References
- Yu, T., Wu, L., & Gilbert, P. B. (2018). *A joint model for mixed and truncated longitudinal data and survival data*. Biostatistics, 19(3), 374–390.
- Saha et al. (2025). *Quantifying the waning of humoral immunity*. medRxiv. https://doi.org/10.1101/2025.05.13.25327542


## Acknowledgement

This work builds upon the HHJMs framework:  
https://github.com/oliviayu/HHJMs  

The original implementation has been extended to support nonlinear longitudinal modeling and interval-censored survival outcomes.
