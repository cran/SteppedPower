## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  warning = FALSE
)
library(knitr)
library(SteppedPower)
library(Matrix)
library(plotly)
assign("knit_print.plotly", function(x, ...) {
  x <- plotly::partial_bundle(x)
  NextMethod()
}, envir = .GlobalEnv)

# A function for captioning and referencing images
fig <- local({
    i <- 0
    ref <- list()
    list(
        cap=function(refName, text) {
            i <<- i + 1
            ref[[refName]] <<- i
            paste("Figure ", i, ": ", text, sep="")
        },
        ref=function(refName) {
            paste("(Fig.", ref[[refName]],")")
        })
})

## ----design-plot, fig.cap=fig$cap("design","Design matrix for 4 sequences with 1 cluster each. 5 periods."), fig.width=2, fig.height=1, echo=FALSE----
dsn_bin <- construct_DesMat(Cl = rep(1, 4), N = 50)
plot(dsn_bin)

## -----------------------------------------------------------------------------
p <- 0.2; effect <- 0.15; icc <- 0.2
p_mid <- p + effect * 0.5

DM <- construct_DesMat(Cl = rep(1,4), N = 50)
pow_re <- glsPower(DesMat = DM, mu0 = p, mu1 = p + effect, 
                   tau = sqrt(p_mid * (1 - p_mid) * icc / (1 - icc)),
                   family = "binomial", verbose = 2)

pow_al <- glsPower(DesMat = DM, mu0 = p, mu1 = p + effect,
                   alpha_0_1_2 = c(icc, icc),
                   family = "binomial", verbose = 2)

c(power_randomEffects = pow_re$power,
  power_alpha012      = pow_al$power)

## ----fig.cap=fig$cap("covmat","Covariance matrix under random effects (left) and alpha_0_1_2 (right), both targeting ICC = 0.2.")----
plotly::subplot(
  plot(pow_re, which = 4, show_colorbar = FALSE)$CMplot,
  plot(pow_al, which = 4, show_colorbar = FALSE)$CMplot,
  nrows = 1, margin = 0.05)

## ----variance computation, echo=FALSE-----------------------------------------
# data-raw/binom_outcomes_grid.R for the computation.
# vals <- readRDS("inst/vignettes/BinomOutcomes_vals.rds")
vals <- readRDS(system.file("vignettes", 
                            "BinomOutcomes_vals.rds",
                            package = "SteppedPower"))

p_vals <- unique(vals$p0)
reMat  <- matrix(vals[, "re"],  length(p_vals), length(p_vals))
alMat  <- matrix(vals[, "al"],  length(p_vals), length(p_vals))
linMat <- matrix(vals[, "lin"], length(p_vals), length(p_vals))

## ----contour_diff, fig.width=6.5, fig.height=5.5, echo=FALSE, fig.cap=fig$cap("contour_midpoint","Relative difference of variance estimators. ICC=0.2, N=50, 4 sequences. blue: midpoint method has higher variance than cell-by-cell method.")----
diff_mat <- (reMat - linMat)/linMat * 100
maxabs <- max(abs(diff_mat), na.rm = TRUE)
plot_ly(
  x = p_vals, y = p_vals, z = diff_mat,
  type = "contour",
  colorscale = list(list(0,   "steelblue"),
               list(0.5, "white"),
               list(1,   "firebrick") ),
  zmin = -maxabs, zmax = maxabs,
  contours = list(coloring = "fill", showlabels = TRUE),
  line  = list(smoothing = 1.3) ) |>
    add_trace(x = range(p_vals), y = range(p_vals), type = "scatter",
                mode = "lines", line = list(color = "gray", dash = "dash"),
                showlegend = FALSE, hoverinfo = "none") |>
    layout( xaxis = list(title = "p0"),
            yaxis = list(title = "p1") )

## ----contour_diff2, fig.width=6.5, fig.height=5.5, echo=FALSE, fig.cap=fig$cap("contour_alpha","Relative difference of variance estimators. ICC=0.2, N=50, 4 sequences. blue: alpha method has higher variance than random effects method.")----
diff_mat <- (reMat - alMat)/alMat * 100
maxabs <- max(abs(diff_mat), na.rm = TRUE)
plot_ly(
  x = p_vals, y = p_vals, z = diff_mat,
  type = "contour",
  colorscale = list(list(0,   "steelblue"),
               list(0.5, "white"),
               list(1,   "firebrick") ),
  zmin = -maxabs, zmax = maxabs,
  contours = list(coloring = "fill", showlabels = TRUE),
  line  = list(smoothing = 1.3)) |>
    add_trace(x = range(p_vals), y = range(p_vals), type = "scatter",
                   mode = "lines", line = list(color = "gray", dash = "dash"),
                   showlegend = FALSE, hoverinfo = "none") |>
    layout( xaxis = list(title = "p0"),
            yaxis = list(title = "p1") )

