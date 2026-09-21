# `SteppedPower` - Power Calculation for Stepped Wedge Designs

<!-- badges: start -->
[![R-CMD-check](https://github.com/PMildenb/SteppedPower/workflows/R-CMD-check/badge.svg)](https://github.com/PMildenb/SteppedPower/actions)
[![CRAN version](https://img.shields.io/cran/v/SteppedPower)](https://cran.r-project.org/package=SteppedPower)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/SteppedPower)](https://cran.r-project.org/package=SteppedPower)
<!-- badges: end -->

**Tools for power and sample size calculation, as well as design diagnostics.  
For longitudinal mixed model settings, with a focus on stepped wedge designs.**

`SteppedPower` provides power and sample size calculation for parallel, crossover, and stepped wedge designs. It allows for a flexible definition of the covariance structure.
It further offers visualisations and diagnostics tools, to assess cluster importance across time points. 

## Installation

### Stable release (CRAN)
```r
install.packages("SteppedPower")
```

### Development versions
```r
devtools::install_github("PMildenb/SteppedPower", build_vignettes = TRUE)   ## stable development version
devtools::install_github("PMildenb/SteppedPower", ref = "devel", build_vignettes = TRUE)  ## latest 
```

## Quick Start

```r
library(SteppedPower)

# SWD with 4 clusters, ICC = 0.1 (via tau), 10 subjects per cluster, and a treatment effect of 0.5
result <- glsPower(
  Cl = rep(1, 4),       # 4 clusters in 4 sequences
  mu0 = 0,              # Mean under control
  mu1 = 0.5,            # Mean under treatment
  sigma = 1,            # Residual standard deviation
  tau = sqrt(0.111),    # Random intercept SD (ICC = tau^2 / (tau^2 + sigma^2) ≈ 0.1)
  N = 10,               # Subjects per cluster
  verbose = 2           # Save additional info, e.g., complete covariance matrix
)

# View power calculation
print(result)

# check the design matrix
plot(result$DesignMatrix)

# check the covariance matrix
plot(result$CovarianceMatrix)

# check influence diagnostics 
plot(result)
```
## Documentation

For more details, see the package vignettes:
```r
vignette("Getting_Started", package = "SteppedPower")
```
## References

- Hussey S, Hughes JP (2007). "Design and analysis of stepped wedge cluster randomised trials." *Contemporary Clinical Trials*, 28(2), 182-191. <doi:10.1016/j.cct.2006.05.007>
- Li F, et al. (2020). "Mixed model sample size calculations for stepped wedge cluster randomised trials." *Statistical Methods in Medical Research*. <doi:10.1177/0962280220932962>

## Authors

- Philipp Mildenberger ([ORCID: 0000-0002-7367-1708](https://orcid.org/0000-0002-7367-1708)) - `pmildenb@uni-mainz.de`
- Federico Marini ([ORCID: 0000-0003-3252-7758](https://orcid.org/0000-0003-3252-7758))

## License

MIT