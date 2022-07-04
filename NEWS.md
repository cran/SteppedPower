

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
