
# `SteppedPower` - Sample Size Calculation in Mixed Model Settings with Focus on Stepped Wedge Designs

<!-- badges: start -->
  [![R-CMD-check](https://github.com/PMildenb/SteppedPower/workflows/R-CMD-check/badge.svg)](https://github.com/PMildenb/SteppedPower/actions)
  <!-- badges: end -->

SteppedPower offers tools for power and sample size calculation as well as 
design diagnostics for 
longitudinal mixed model settings, with a focus on stepped wedge designs.
All calculations are oracle estimates i.e. assume random effect variances 
to be known (or guessed) in advance. 


## Installation
To install the development version type    
`devtools::install_github("PMildenb/SteppedPower", build_vignettes=TRUE)`  
`library(SteppedPower)`

## Vignette 
The vignette can be viewed with
`vignette("Getting_Started_with_SteppedPower")`
