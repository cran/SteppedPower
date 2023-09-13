## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(SteppedPower)

## -----------------------------------------------------------------------------
Dsn1.1 <- construct_DesMat(Cl=rep(2,4), incomplete=2)

## -----------------------------------------------------------------------------
TM  <- toeplitz(c(1,1,0,0))
incompleteMat1 <- cbind(TM[,1:2],rep(1,4),TM[,3:4])
incompleteMat2 <- incompleteMat1[rep(1:4,each=2),]

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(incompleteMat1))

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(incompleteMat2))

## -----------------------------------------------------------------------------
Dsn1.2 <- construct_DesMat(Cl=rep(2,4), incomplete=incompleteMat1)
Dsn1.3 <- construct_DesMat(Cl=rep(2,4), incomplete=incompleteMat2)

all.equal(Dsn1.1,Dsn1.2)
all.equal(Dsn1.1,Dsn1.3)

## -----------------------------------------------------------------------------
Dsn2 <- construct_DesMat(Cl=rep(2,4), trtDelay = c(NA) )
Dsn2

## -----------------------------------------------------------------------------
Dsn3 <- construct_DesMat(Cl=rep(2,4), incomplete=2, trtDelay=c(NA) )
Dsn3

