

# SteppedPower 0.4.0

* `glsPower()` now supports count outcomes via `family="poisson"`
* Added ICC (intracluster correlation) transformation functions: `icc_to_RandEff()`, `RandEff_to_icc()`, `RandEff_to_alpha012()`, and `alpha012_to_RandEff()` for converting between random effects variances and ICC/CAC/IAC parameters
* Added vignette on binomial and count outcomes, with pre-calculated contour plots
* In `glsPower()`, the argument `N` now overrides the `N` stored in a supplied
`DesMat` object
* `glsPower()` now fails gracefully (with a warning) if the information content
cannot be calculated
* Character input arguments (e.g. `dsntype`, `family`) now throw an error if no
known option is sufficiently similar, instead of silently choosing the closest match
* Diagnostic output now uses `message()` instead of `print()`
* Vignette plots now use the plotly partial bundle to reduce package size
* Changed covariance matrix construction to use `fbdiag` (fast block diagonal matrix)
* `plot_CellWeights()` now treats `NA` entries in `incompMat` as unobserved
cluster periods
* Replaced `\()` with `function()` for backward compatibility with older R versions
* Fixed typo in `RandEff_to_alpha`
* Added tests for `construct_DesMat()`, `construct_CovMat()`, and `glsPower()`
* Updated vignettes and improved documentation with additional links in help files
* Roxygen documentation now uses markdown format; re-roxygenised all documentation


# SteppedPower 0.3.5

* Addressed CRAN comments
* Fixed roxygen package name bug


# SteppedPower 0.3.4

* Addressed CRAN comments


# SteppedPower 0.3.3

* `N` (subjects per cluster-period cell) now belongs to `DesMat` class
* Added vignette for incomplete designs
* Fixed bug for handling `NA` in `incompMat` and `trtMat`


# SteppedPower 0.3.2

* The most noticeable change in this version is that the abbrevation `wls` 
(weighted least squares) in function names is now replaced with `gls`
(generalised least squars) to more properly reflect the scope of the functionality.
For example, the function `wlsPower()` is now called `glsPower()` - although the
former version still works and throws a warning. 
* The closed formula for the computation of information content is now a dedicated formula, 
called `compute_InfoContent()` 
* In `plot.glsPower()` there now is an option to manually set the font size of the
annotation in the influence plots


# SteppedPower 0.3.1

* The function `wlsPower()` now also computes the information content of 
cluster-period cells. Computation is currently done twice, once with a general formula
and once explicitly. Information content of whole periods or clusters is also computed.
* The method `plot.wlsPower()` recieved multiple updates:
  * It now produces up to four plots: the projection matrix, 
  the information content, the intervention design and the covariance matrix.
  * Incomplete designs (SWD where some cluster-period cells are omitted) are now visualised
  * Plots of projection matrix and information content can now be annotated with particular values in each cell;
  This is the default for smaller designs and can be turned on/off via `annotations = <TRUE/FALSE>`
  * An option `show_colorbar` to hide colour bars was added
  * An option `marginal_plots` to hide marginal plots on whole periods or clusters was added.
  * Various aesthetic improvements, e.g.: Improved hover information, dynamic gap size between cells.
* Vignette extended
  


# SteppedPower 0.2.0 

* The function `wlsPower()` now has an argument `alpha_012` that offers an alternative
way to specifiy the correlation matrix.
* In function `wlsPower()`, the argument `AR` now accepts a vector of up to three values. 
This allows to specifiy autoregressive structures for only a subset of: random cluster intercept, 
random intervention effect and random subject intercept. 
* Closed formulae were added. 
* The method `plot.wlsPower` now produces up to three plots, the projection matrix, the intervention design and the covariance matrix.
* The vignette was extended.


# SteppedPower 0.1.0

* Initial submission to CRAN
