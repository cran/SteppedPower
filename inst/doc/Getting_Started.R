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

## -----------------------------------------------------------------------------
wlsPwr <- wlsPower(Cl=c(3,2,3), mu0=0, mu1=1, sigma=1, tau=.5, verbose=2)
plot(wlsPwr)

## -----------------------------------------------------------------------------
wlsPower(Cl=c(3,2,3), mu0=0, mu1=.2, sigma=1, tau=0, Power=.8)

## -----------------------------------------------------------------------------
wlsPower(Cl=c(10,10), mu0=0,mu1=.6,sigma=1, tau=0, N=1, 
              dsntype="parallel", timepoints=1)

## the same:
wlsPower(Cl=c(1,1), mu0=0,mu1=.6, sigma=1, tau=0, N=10,
              dsntype="parallel", timepoints=1)

## -----------------------------------------------------------------------------
pwr::pwr.norm.test(.3,n=20)$power

## -----------------------------------------------------------------------------
wlsPower(Cl=c(10,10),timepoints=5,mu0=0,mu1=.25,
         sigma=.5,dsntype="parallel")

wlsPower(Cl=c(10,10),timepoints=5,mu0=0,mu1=.25,
         sigma=.5,tau=.2,dsntype="parallel")

## ---- warning=FALSE-----------------------------------------------------------
mod1 <- wlsPower(Cl=c(1,1,1,0), mu0=0, mu1=1, 
                 sigma=0.4, tau=0, verbose=2)

knitr::kable(mod1$DesignMatrix$trtMat)

## -----------------------------------------------------------------------------
mod2 <- wlsPower(Cl=c(2,2,2,2), mu0=0, mu1=1, 
              sigma=1, N=100, tau=1, tauAR=.6, verbose=2)

mod3 <- wlsPower(Cl=c(2,2,2,2), mu0=0, mu1=1, 
              sigma=1, N=100, tau=1, tauAR=.95, verbose=2)

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(as.matrix(mod2$CovarianceMatrix[1:5,1:5])))

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(as.matrix(mod3$CovarianceMatrix[1:5,1:5])))
rm(mod1,mod2,mod3)

## -----------------------------------------------------------------------------
mod4 <- wlsPower(Cl=c(1,1,1), mu0=0, mu1=1, N=c(1,3,10), tau=.5, verbose=2)
plot(mod4)

## ---- echo=FALSE--------------------------------------------------------------
rm(mod4)

## -----------------------------------------------------------------------------
incompletePwr <- wlsPower(Cl=rep(2,4), sigma=2, tau=.6, mu0=0,mu1=.5, N=80, 
                             incomplete=2, verbose=2)
incompletePwr

## -----------------------------------------------------------------------------
TM  <- toeplitz(c(1,1,0,0))
incompleteMat1 <- cbind(TM[,1:2],rep(1,4),TM[,3:4])
incompleteMat2 <- incompleteMat1[rep(1:4,each=2),]

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(incompleteMat1))

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(incompleteMat2))

## -----------------------------------------------------------------------------
incompletePwr1 <- wlsPower(Cl=rep(2,4), sigma=2, tau=.6, mu0=0, mu1=.5, N=80, 
                        incomplete=incompleteMat1, verbose=2)
incompletePwr2 <- wlsPower(Cl=rep(2,4), sigma=2, tau=.6, mu0=0, mu1=.5, N=80, 
                        incomplete=incompleteMat2, verbose=2)

all.equal(incompletePwr,incompletePwr1)
all.equal(incompletePwr,incompletePwr2)

## -----------------------------------------------------------------------------
plot(incompletePwr)

## ---- echo=FALSE--------------------------------------------------------------
rm(incompletePwr,incompletePwr1,incompletePwr2,incompleteMat1,incompleteMat2)

## -----------------------------------------------------------------------------
TimeAdj1 <- wlsPower(Cl=rep(2,4), mu0=0, mu1=1, tau=0, 
                     timeAdjust="linear", verbose=2)

TimeAdj2 <- wlsPower(Cl=rep(2,4), mu0=0, mu1=1, tau=0, 
                     timeAdjust="factor", verbose=2)

## -----------------------------------------------------------------------------
knitr::kable(head(TimeAdj1$DesignMatrix$dsnmatrix, 5))

## -----------------------------------------------------------------------------
knitr::kable(head(TimeAdj2$DesignMatrix$dsnmatrix, 5))

## -----------------------------------------------------------------------------
Closed1 <- wlsPower(mu0=0, mu1=5, Cl=rep(3,3), sigma=5, tau=1, psi=2, gamma=1,
                      N=3, verbose=2)
a <- plot(Closed1$DesignMatrix)

## -----------------------------------------------------------------------------
Closed2 <- wlsPower(mu0=0, mu1=5, Cl=rep(3,3), sigma=5, tau=1, psi=2, gamma=1,
                      N=3, verbose=2, INDIV_LVL = TRUE)
plot(Closed2)

## -----------------------------------------------------------------------------
sessionInfo()

