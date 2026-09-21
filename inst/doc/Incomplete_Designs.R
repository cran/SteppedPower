## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>",
                     fig.width = 7, fig.height = 2, warning = FALSE)
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

## ----echo=FALSE, fig.cap=fig$cap("Incomp1","An incomplete stepped wedge design with four clusters with a total study duration of six periods. Clusters are only observed from two periods prior to treatment switch to two periods post switch.")----
plot(construct_DesMat(Cl=c(1,1,1,1,0), incomplete = 2))

## ----echo=FALSE, fig.cap=fig$cap("Incomp2","An incomplete stepped wedge design with four clusters with a total study duration of six periods. The first period after the switch to interventional treatment - the transition period - is not observed in each cluster.")----
plot(construct_DesMat(Cl=c(1,1,1,1,0), trtDelay = c(NA)))

## -----------------------------------------------------------------------------
Dsn1.1 <- construct_DesMat(Cl=rep(2,4), incomplete=2)

## -----------------------------------------------------------------------------
TM  <- toeplitz(c(1,1,0,0))
incompleteMat1 <- cbind(TM[,1:2],rep(1,4),TM[,3:4])
incompleteMat2 <- incompleteMat1[rep(1:4,each=2),]

## ----echo=FALSE---------------------------------------------------------------
suppressWarnings(knitr::kable(incompleteMat1))

## ----echo=FALSE---------------------------------------------------------------
suppressWarnings(knitr::kable(incompleteMat2))

## -----------------------------------------------------------------------------
Dsn1.2 <- construct_DesMat(Cl=rep(2,4), incomplete=incompleteMat1)
Dsn1.3 <- construct_DesMat(Cl=rep(2,4), incomplete=incompleteMat2)

all.equal(Dsn1.1$trtMat,Dsn1.2$trtMat)
all.equal(Dsn1.1$trtMat,Dsn1.3$trtMat)

## -----------------------------------------------------------------------------
Dsn2 <- construct_DesMat(Cl=rep(2,4), trtDelay = c(NA) )
Dsn2

## -----------------------------------------------------------------------------
Dsn3 <- construct_DesMat(Cl=rep(2,4), incomplete=2, trtDelay=c(NA) )
Dsn3

