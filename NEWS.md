
# SteppedPower 0.2.0 

* The function `wlsPower()` now has an argument `alpha_012` that offers an alternative
way to specifiy the correlation matrix.
* In function `wlsPower()`, the argument `AR` now accepts a vector of up to three values. 
This allows to specifiy autoregressive structures for only a subset of: random cluster intercept, 
random intervention effect and random subject intercept. 
* Closed formulae were added. 
* The method `plot.wlsPower` now produces up to three plots, the hatmatrix, the intervention design and the covariance matrix.
* The vignette was extended.


# SteppedPower 0.1.0

* Initial submission to CRAN
