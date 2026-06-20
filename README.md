# Simultaneous Inference of Susceptibility Distributions and Intervention Effects from Epidemic Curves

This repository contains the R code and analysis for "On the simultaneous inference of susceptibility distributions and intervention effects from epidemic curves" by Ibrahim Mohammed, Chris Robertson, and M. Gabriela M. Gomes.

## Overview

This research investigates parameter identifiability in SEIR epidemic models, with emphasis on:

1. **Susceptibility heterogeneity**: How individual variation in susceptibility affects epidemic dynamics and parameter estimation
2. **Non-pharmaceutical interventions (NPIs)**: Modeling the effects of public health interventions on transmission
3. **Parameter identifiability**: Assessing correlations between model parameters and strategies to improve parameter estimation
4. **Multiple epidemic analysis**: Using data from concurrent epidemics to enhance parameter identifiability

The study demonstrates that parameter identifiability issues exist in both homogeneous and heterogeneous SEIR models, and shows how using multiple epidemic datasets can significantly improve parameter estimation.

## Repository Structure

```
epidemic-susceptibility-inference/
├── README.md                             # This file
├── LICENSE.txt                           # MIT License
├── .gitignore                            # Git ignore rules
├── install_packages.R                    # Script to install required packages
├── R/                                    # Core functions
│   ├── MaxLik_functions_reduced_model.R  # Maximum likelihood fitting functions
│   └── utility_functions.R                   # Utility functions for modeling
├── scripts/                              # Analysis scripts
│   ├── 1_baseline_cases.R                # Baseline analysis (Case I and II)
│   ├── 2_mle_single_epidemic.R           # MLE profiling for single epidemic
│   ├── 3_mle_two_epidemics.R             # MLE profiling for two epidemics
│   ├── 4_single_epidemic_correlation.R   # Initial condition impact analysis (single)
│   └── 5_two_epidemics_correlation.R     # Initial condition impact analysis (two epidemics)
├── results/                              # Output data files
│   └── .gitkeep                          # Track directory in git
├── figures/                              # Generated plots and visualizations
│   └── .gitkeep                          # Track directory in git
└── paper/                                # Publication files
    ├── ibrahim_1wave_epidemics.pdf        # Main manuscript
    └── ibrahim_1wave_SM.pdf             # Supplementary material
```

## Key Models and Methods

### SEIR Model with Heterogeneity
- **Homogeneous SEIR**: Classical model without individual variation (ν = 0)
- **Heterogeneous SEIR**: Model incorporating gamma-distributed susceptibility variation
- **Reduced form**: Analytically derived system for gamma-distributed heterogeneity

### Parameter Estimation
- **Maximum Likelihood Estimation (MLE)**: Primary method for parameter fitting
- **Profile likelihood**: Analysis of parameter confidence intervals
- **Correlation analysis**: Assessment of parameter identifiability through correlation matrices

### Key Parameters
- **R_0**: Basic reproduction number
- **\nu **: Heterogeneity parameter (coefficient of variation of susceptibility)
- **t_0**: Time onset of behavioural change before official lockdown when people begin to social distance
- **c_1**: Intervention strength parameter

## Requirements

To run this code, you'll need:

- **R version**: 4.0.0 or higher
- **Required packages**:
  - `deSolve` - For solving differential equations

  - `ggplot2` - For plotting
  - `dplyr` - For data manipulation
  - `tidyr` - For data tidying
  - `gridExtra` - For arranging multiple plots
  - `MASS` - For statistical functions

Install all required packages by running:

```r
source("install_packages.R")
```

## Usage

### Quick Start

1. **Clone the repository**:
```bash
git clone https://github.com/yourusername/On-the-simultaneous--inference--of--epidemic-susceptibility-distributions-and-NPI.git
cd epidemic-susceptibility-inference
```

2. **Install dependencies**:
```r
source("install_packages.R")
```

3. **Run analysis scripts in order**:
```r
# Run baseline analysis
source("scripts/1_baseline_cases.R")

# Single epidemic MLE profiling  
source("scripts/2_mle_single_epidemic.R")

# Two epidemics MLE profiling
source("scripts/3_mle_two_epidemics.R")

# Correlation analysis for single epidemic
source("scripts/4_single_epidemic_correlation.R")

# Correlation analysis for two epidemics
source("scripts/5_two_epidemics_correlation.R")
```

### Script Descriptions

- **`1_baseline_cases.R`**: Implements baseline analysis comparing homogeneous and heterogeneous models under different scenarios (Cases I and II)

- **`2_mle_single_epidemic.R`**: Performs maximum likelihood estimation and profile likelihood analysis for single epidemic datasets

- **`3_mle_two_epidemics.R`**: Extends MLE analysis to concurrent two-epidemic scenarios to improve parameter identifiability

- **`4_single_epidemic_correlation.R`**: Investigates the impact of initial conditions on parameter correlations in single epidemic fits

- **`5_two_epidemics_correlation.R`**: Analyzes how initial condition variations affect parameter correlations when fitting to two concurrent epidemics

## Key Findings

1. **Parameter correlations exist in both homogeneous and heterogeneous SEIR models**, not just those with heterogeneity parameters

2. **Multiple epidemic data significantly improves parameter identifiability** by reducing correlations between key parameters

3. **Initial condition diversity** across epidemics is crucial for breaking parameter correlations

4. **Profile likelihood analysis** reveals improved confidence intervals when using multiple epidemic datasets

## Scripts Files 

Original files → Repository organization:
 `scripts/1_baseline_cases.R`
- `scripts/2_mle_single_epidemic.R`
`scripts/3_mle_two_epidemics.R`
 `scripts/4_single_epidemic_correlation.R`
`scripts/5_two_epidemics_correlation.R`
`R/MaxLik_fit_functions_reduced_model.R`
 `R/utility_functions.R`

## Citation

If you use this code or findings in your research, please cite:

```bibtex
@article{MOHAMMED2026100911,
title = {On the simultaneous inference of susceptibility distributions and intervention effects from epidemic curves},
journal = {Epidemics},
volume = {55},
pages = {100911},
year = {2026},
issn = {1755-4365},
doi = {https://doi.org/10.1016/j.epidem.2026.100911},
url = {https://www.sciencedirect.com/science/article/pii/S1755436526000277},
author = {Ibrahim Mohammed and Chris Robertson and M. Gabriela M. Gomes},
keywords = {Individual variation, Heterogeneity, Epidemic model, Parameter estimation, Identifiability},
abstract = {Susceptible–Exposed–Infectious–Recovered (SEIR) models with inter-individual variation in susceptibility or exposure to infection were proposed early in the COVID-19 pandemic as a potential element of the mathematical/statistical toolset available to policy development. In comparison with other models employed at the time, those designed to fully estimate the effects of such heterogeneity tended to predict small epidemic waves and hence require less containment to achieve the same outcomes. However, these models never made it to mainstream COVID-19 policy making due to lack of prior validation of their inference capabilities. Here we report the results of the first systematic investigation of this matter in idealized scenarios created with synthetic data. We simulate datasets using the model with strategically chosen parameter values, and then conduct maximum likelihood estimation to assess how well we can retrieve the assumed parameter values. Parameter uncertainties were found to markedly reduce when concurrently fitting multiple epidemics with shared parameters, suggesting a general methodological approach that can be further developed to tackle real-world questions.}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.

## Authors and Affiliations

- **Ibrahim Mohammed**
  - Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK
  - Department of Mathematical Sciences, Abubakar Tafawa Balewa University, Bauchi, Nigeria
  - Funded by Petroleum Development Fund (PTDF), Nigeria

- **Chris Robertson**
  - Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK  
  - Public Health Scotland, Glasgow, UK

- **M. Gabriela M. Gomes**
  - Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK
  - Centre for Mathematics and Applications (NOVA MATH), NOVA School of Science and Technology, Caparica, Portugal
  - Partially funded by FCT – Fundação para a Ciência e a Tecnologia, I.P.

## Support


For questions about the code or methods, please open an issue on this repository or contact the authors through their institutional affiliations.



