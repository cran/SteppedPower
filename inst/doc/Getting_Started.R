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

## -----------------------------------------------------------------------------
glsPwr <- glsPower(Cl=c(3,2,3), mu0=0, mu1=1, sigma=1, tau=.5, verbose=2)
plot(glsPwr,which=1, show_colorbar=FALSE)$WgtPlot 

## -----------------------------------------------------------------------------
plot(glsPwr,which=2, show_colorbar=FALSE)$ICplot

## -----------------------------------------------------------------------------
glsPower(Cl=c(3,3,3), mu0=0, mu1=.2, sigma=1, tau=0, power=.8)

## -----------------------------------------------------------------------------
glsPower(Cl=c(10,10), mu0=0,mu1=1.2,sigma=1, tau=0, N=1, 
              dsntype="parallel", timepoints=1)$power

## the same:
glsPower(Cl=c(1,1), mu0=0,mu1=1.2, sigma=1, tau=0, N=10,
              dsntype="parallel", timepoints=1)$power

## -----------------------------------------------------------------------------
pwr::pwr.norm.test(.6,n=20)$power

## -----------------------------------------------------------------------------
glsPower(Cl=c(10,10),timepoints=5,mu0=0,mu1=.25,
         sigma=.5,dsntype="parallel")

glsPower(Cl=c(10,10),timepoints=5,mu0=0,mu1=.25,
         sigma=.5,tau=.2,dsntype="parallel")

## ---- warning=FALSE-----------------------------------------------------------
mod1 <- glsPower(Cl=c(1,1,1,0), mu0=0, mu1=1, 
                 sigma=0.4, tau=0, verbose=2)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(mod1$DesignMatrix$trtMat)

## -----------------------------------------------------------------------------
mod2 <- glsPower(Cl=c(2,2,2,2), mu0=0, mu1=1, 
              sigma=1, N=100, tau=1, AR=.6, verbose=2)

mod3 <- glsPower(Cl=c(2,2,2,2), mu0=0, mu1=1, 
              sigma=1, N=100, tau=1, AR=.95, verbose=2)

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(as.matrix(mod2$CovarianceMatrix[1:5,1:5])))

## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(knitr::kable(as.matrix(mod3$CovarianceMatrix[1:5,1:5])))
rm(mod1,mod2,mod3)

## -----------------------------------------------------------------------------
mod4 <- glsPower(Cl=c(1,1,1), mu0=0, mu1=1, N=c(1,3,10), 
                 sigma=1, tau=.5, verbose=2)
plot(mod4, which=2, show_colorbar=FALSE)$ICplot

## ---- echo=FALSE--------------------------------------------------------------
rm(mod4)

## -----------------------------------------------------------------------------
incompletePwr <- glsPower(Cl=rep(2,4), sigma=2, tau=.6, mu0=0,mu1=.5, N=80, 
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
incompletePwr1 <- glsPower(Cl=rep(2,4), sigma=2, tau=.6, mu0=0, mu1=.5, N=80, 
                        incomplete=incompleteMat1, verbose=2)
incompletePwr2 <- glsPower(Cl=rep(2,4), sigma=2, tau=.6, mu0=0, mu1=.5, N=80, 
                        incomplete=incompleteMat2, verbose=2)

all.equal(incompletePwr,incompletePwr1)
all.equal(incompletePwr,incompletePwr2)

## -----------------------------------------------------------------------------
plot(incompletePwr, show_colorbar=FALSE)$WgtPlot

## ---- echo=FALSE--------------------------------------------------------------
rm(incompletePwr,incompletePwr1,incompletePwr2,incompleteMat1,incompleteMat2)

## -----------------------------------------------------------------------------
TimeAdj1 <- glsPower(Cl=rep(2,4), mu0=0, mu1=1, sigma=1, tau=0, 
                     timeAdjust="linear", verbose=2)

TimeAdj2 <- glsPower(Cl=rep(2,4), mu0=0, mu1=1, sigma=1, tau=0, 
                     timeAdjust="factor", verbose=2)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(TimeAdj1$DesignMatrix$dsnmatrix, 5))

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(TimeAdj2$DesignMatrix$dsnmatrix, 5))

## -----------------------------------------------------------------------------
Closed1 <- glsPower(mu0=0, mu1=5, Cl=rep(3,3), 
                    sigma=5, tau=1, psi=3,
                    N=3, verbose=2)
Closed1

## -----------------------------------------------------------------------------
Closed2 <- glsPower(mu0=0, mu1=5, Cl=rep(3,3), 
                    sigma=5, tau=1, psi=3,
                    N=3, verbose=2, INDIV_LVL = TRUE)
Closed2
plot(Closed2, annotations=FALSE, show_colorbar=FALSE)$WgtPlot

## -----------------------------------------------------------------------------
Closed1$power - Closed2$power

## -----------------------------------------------------------------------------
Open1 <- glsPower(mu0=0, mu1=5, Cl=rep(3,3), 
                  sigma=5, tau=1, psi=3, AR=c(1,1,.75), N=3)

Closed1$power
Open1$power

## ---- fig.height=6------------------------------------------------------------
Open2Indiv <- glsPower(mu0=0, mu1=10, Cl=c(1,1,1,0), 
                       sigma=1, tau=5, psi=10, AR=c(1,1,.60),
                       N=3, verbose=2, INDIV_LVL=TRUE)

plot(Open2Indiv, which=4, show_colorbar=FALSE)$CMplot

## -----------------------------------------------------------------------------
trtMat <- construct_DesMat(c(6,6,6,6))$trtMat
mu0 <- 0.05 ; mu1 <- 0.032 ; N <- 100 

tau <- .025 ; sigma <- sqrt(.041*.959) 
gamma <- 0.01 ; psi <- .1 ; chi <- .5 ; AR <- .5

## -----------------------------------------------------------------------------
tmp <- VarClosed_Li(trtMat, tau=tau, psi=psi, N=N, AR=AR)
tTestPwr(mu0-mu1, se=sqrt(tmp), df=Inf)

## -----------------------------------------------------------------------------
a <- SteppedPower::glsPower(Cl=rep(6,4), mu0=mu0, mu1=mu1, AR=AR,
                       sigma=0, tau=tau, N=N, psi=psi, verbose=1, INDIV_LVL = TRUE)
a

## -----------------------------------------------------------------------------
tmp <- VarClosed_Kasza(trtMat, sigma=sigma, tau=tau, gamma=gamma, psi=psi, N=N, chi=0)
tTestPwr(mu0-mu1, se=sqrt(tmp), df = Inf)
glsPower(Cl = rep(6,4), N=N, mu0=mu0, mu1=mu1, verbose=0,
         sigma=sigma, tau=tau, gamma=gamma, psi=psi)

tmp <- VarClosed_Kasza(trtMat, sigma=sigma, tau=tau, gamma=gamma, psi=psi, N=N, chi=1)
tTestPwr(mu0-mu1, sqrt(tmp), df = Inf)
glsPower(Cl = rep(6,4), N=N, mu0=mu0, mu1=mu1, verbose=0,
         sigma=sigma, tau=tau, 
         gamma=sqrt(gamma^2+psi^2/N), psi=0)

tmp <- VarClosed_Kasza(trtMat, sigma=sigma, tau=tau, gamma=gamma, psi=psi, N=N, chi=chi)
tTestPwr(mu0-mu1, sqrt(tmp), df = Inf)
glsPower(Cl = rep(6,4), N=N, mu0=mu0, mu1=mu1, verbose=0,
         sigma=sigma, tau=tau, 
         gamma=sqrt(gamma^2+chi*psi^2/N), psi=sqrt(1-chi)*psi)

## ---- echo=FALSE--------------------------------------------------------------
print(sessionInfo(),locale=FALSE)

