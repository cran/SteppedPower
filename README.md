
# `SteppedPower` - Power Calculation for Stepped Wedge Designs

<!-- badges: start -->
  [![R-CMD-check](https://github.com/PMildenb/SteppedPower/workflows/R-CMD-check/badge.svg)](https://github.com/PMildenb/SteppedPower/actions)
  <!-- badges: end -->

SteppedPower offers tools for power and sample size calculation as well as 
design diagnostics for 
longitudinal mixed model settings, with a focus on stepped wedge designs.
All calculations are oracle estimates i.e. assume random effect variances 
to be known (or guessed) in advance. 


## Installation

Install from CRAN with
`install.packages("SteppedPower")`. 
Current version on CRAN is 0.2.0.

To install the latest stable version type    
`devtools::install_github("PMildenb/SteppedPower", build_vignettes=TRUE)`  

The development version with the most recent changes can be installed with
`devtools::install_github("PMildenb/SteppedPower", ref='devel', build_vignettes=TRUE)`


## Vignette 
The vignette can be viewed with
`vignette("Getting_Started", package="SteppedPower")`
