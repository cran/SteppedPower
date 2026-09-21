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

## -----------------------------------------------------------------------------
glsPwr <- glsPower(Cl=c(3,2,3), mu0=0, mu1=1, sigma=1, tau=.5, verbose=2)


## ----echo=FALSE, fig.height=8, fig.cap=fig$cap("Influence_Plot_1","Treatment allocation plot (top). Influence of cluster-period cells (center) and information content (bottom) of a stepped wedge design with 8 clusters in 3 sequences")----

tmpplt <- plot(glsPwr,which=1:3, 
               marginal_plots = FALSE,
               show_colorbar  = FALSE)
subplot(tmpplt[[3]],tmpplt[[1]],tmpplt[[2]], 
        titleY=TRUE,
        nrows=3,
        margin=c(0,0,.03,.03)  )

## ----fig.cap=fig$cap("Info_Plot_1","Information content of cluster-period cells of a stepped wedge design with 8 clusters in 3 sequences.")----
plot(glsPwr,which=2, show_colorbar=FALSE)$IMplot

## -----------------------------------------------------------------------------
glsPower(Cl=c(3,3,3), mu0=0, mu1=.2, sigma=1, tau=0, power=.8)

## -----------------------------------------------------------------------------
mod4 <- glsPower(Cl=c(1,1,1), mu0=0, mu1=1, N=c(1,3,10), 
                 sigma=1, tau=.5, verbose=2)
plot(mod4, which=1, show_colorbar=FALSE)[[1]]

## ----echo=FALSE---------------------------------------------------------------
rm(mod4)

## -----------------------------------------------------------------------------
Incomp1 <- plot(construct_DesMat(Cl=c(1,1,1,1,0), incomplete = 2))
Incomp2 <- plot(construct_DesMat(Cl=c(1,1,1,1,0), trtDelay = c(NA)))

plotly::subplot(Incomp1, Incomp2, nrows=1, titleX=TRUE, titleY=TRUE, margin=0.05) 

## -----------------------------------------------------------------------------
plt1 <- glsPower(Cl=c(1,1,1), mu0=0, mu1=1, 
                 sigma=1, N=1, tau=1, AR=1, verbose=2) |>
  plot( which=4, show_colorbar = FALSE)

plt2 <- glsPower(Cl=c(1,1,1), mu0=0, mu1=1, 
                 sigma=1, N=1, tau=1, AR=.6, verbose=2) |>
  plot( which=4, show_colorbar = FALSE)

plotly::subplot(plt1[[1]], plt2[[1]], 
                nrows=1 , titleX=TRUE, titleY= TRUE)

## -----------------------------------------------------------------------------
TimeAdj1 <- glsPower(Cl=rep(2,4), mu0=0, mu1=1, sigma=1, tau=0, 
                     timeAdjust="linear", verbose=2)

TimeAdj2 <- glsPower(Cl=rep(2,4), mu0=0, mu1=1, sigma=1, tau=0, 
                     timeAdjust="factor", verbose=2)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(head(TimeAdj1$DesignMatrix$dsnmatrix, 5))

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(head(TimeAdj2$DesignMatrix$dsnmatrix, 5))

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
         sigma=.5,tau=.2,dsntype="parallel")

## ----warning=FALSE------------------------------------------------------------
mod1 <- glsPower(Cl=c(1,1,1,0), mu0=0, mu1=1, 
                 sigma=0.4, tau=0, verbose=2)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(mod1$DesignMatrix$trtMat)

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

## ----fig.height=6-------------------------------------------------------------
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

## ----echo=FALSE---------------------------------------------------------------
print(sessionInfo(),locale=FALSE)

