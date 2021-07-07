#' SteppedPower
#'
#' SteppedPower offers tools for power and sample size
#' calculation as well as design diagnostics for
#' longitudinal mixed model settings, with a focus on stepped wedge designs.
#' All calculations are oracle estimates i.e. assume random effect variances
#' to be known (or guessed) in advance.
#'
#' @importFrom stats coef family gaussian optim pnorm qnorm dnorm rbinom rnorm
#' uniroot integrate pt qt binomial
#' @import Matrix
#' @importFrom grDevices colorRamp
#' @importFrom utils adist
#' @importFrom plotly colorbar config layout plot_ly plotly_empty subplot TeX
#' "%>%"
#'
#' @author Philipp Mildenberger \email{pmildenb@@uni-mainz.de}
#' @name SteppedPower-pkg
#' @docType package
NULL
