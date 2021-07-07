
#' @title Compute Power of a Wald Test
#'
#' @description
#' Computes the power of a scaled Wald test given a standard error,
#' an effect size, the degrees of freedom of the t-distribution and
#' a significance level.
#' Computes the exact power, see second example
#'
#' @param d numeric, raw effect
#' @param se numeric, standard error
#' @param df numeric, degrees of freedom of the t-distribution
#' @param sig.level numeric, significance level, defaults to 0.05
#'
#' @return a scalar
#' @export
#'
#' @examples tTestPwr(4,1,10) ; tTestPwr(4,1,30) ; tTestPwr(4,1,Inf)

tTestPwr <- function(d,se,df,sig.level=0.05){
  dsz <- abs(d/se)
  q   <- qt(sig.level/2,df=df)
  Pwr <- pt(dsz + q, df=df) + pt(-dsz + q, df=df)
  return(Pwr)
}

# tTestPwr uses a central t-distribution.
# As applied sometimes in the literature,
# e.g. Li, Turner, Preisser 2018 | Li, Redden 2015
# This is NOT the same as a t-Test, which uses a non-central t-distribution !!
# Hence the name of the above function should be changed to 'scaledWaldPwr'
# tTestPwr2 tries to mimic a t-test when only df are given.

tTestPwr2 <- function(d,se,df,sig.level=.05){
  d   <- abs(d/se)
  q   <- qt(sig.level/2, df=df)
  ncp <- -sqrt((df+2)/4)*d ## TODO: ?? how to calculate n using df ??
  pwr <- pt(q,  df, ncp=ncp) + pt(-q, df,ncp=ncp,lower.tail=FALSE)
  return(pwr)
}

choose_character_Input <- function(Options, Input){
  Options[which.min(adist(Input, Options,
                          costs=c(insertions    = 1,
                                  deletions     = 100,
                                  substitutions = 100),
                          ignore.case=TRUE))]
}



BYROW <- function(input, SumCl, tp){
  len <- length(input)
  out <- if(len %in% c(1,SumCl, SumCl*tp)) FALSE else if(len==tp) TRUE else
           stop(paste('length of input is ', len,
                      '. This does not fit to given number of timepoints, ',
                      'which is ',tp,
                      ', or to the given number of clusters, which is ', SumCl))
  # if(tp==SumCl & len==tp)
  #   warning("Input is assumed to change between clusters. If you wanted it
  #             to change over time, please provide as matrix of dimension
  #             #Cluster x timepoints")
  return(out)
}

input_to_Mat <- function(input, SumCl, tp){
  return(matrix(input, nrow=SumCl, ncol=tp, byrow=BYROW(input,SumCl,tp)))
}

input_to_List <- function(input, SumCl, tp){
  if(is.null(input)){
    out <- vector("list", SumCl)
  } else {
    Mat <- input_to_Mat(input=input, SumCl=SumCl, tp=tp)
    out <- split(Mat, 1:nrow(Mat))
  }
  return(out)
}

################################################################################
## auxiliary functions for binomial outcome ####
logit.deriv      <- function(x) 1/(x-x^2)
sdLin_invlogit   <- function(sdLin,mu){sdLin/logit.deriv(mu)}
sd_to_logit      <- function(sd,mu){sd * logit.deriv(mu)}
logit    <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))
invlogit.deriv <- function(x) exp(-x)/(1+exp(-x))^2

func <- function(tLin, tauLin, muCondLin){
  stats::dnorm(tLin,0,tauLin) * invlogit(muCondLin + tLin)
}
muCond_to_muMarg <- function(muCond, tauLin){
  muCondLin <- logit(muCond)
  out <- if (tauLin<1e-5) { invlogit(muCondLin) } else {
    integrate(func, -Inf, Inf, tauLin=tauLin, muCondLin=muCondLin )
  }
  return(out$value)
}
muMarg_to_muCond <- function(muMarg,tauLin){
  uniroot(function(x,muMarg,tauLin) {muCond_to_muMarg(x,tauLin=tauLin)-muMarg},
          c(1e-6,1-1e-6),
          muMarg=muMarg, tauLin=tauLin,
          tol=5*.Machine$double.eps)$root
}


#
# x <- seq(0,1,length.out=1001)
# plot(x, tau_to_tauLin(x,mu=.13), type="l")
# lines(x,tau_to_tauLin(x,mu=.5), col=3)
# lines(x,tau_to_tauLin(x,mu=.05), col=4)
# lines(x,x,col=2)
#
# binomial()$linkfun(.65) - binomial()$linkfun(.65+.05)
# binomial()$linkfun(.65) - binomial()$linkfun(.65-.05)
#
# tau_to_tauLin(.22,.65)
#
# plot(x,binomial()$linkfun(x), type="l")
# deriv(expression(binomial()$linkfun),"x")
#
# logit <- function(x) log(x/(1-x))
# a <- Deriv::Deriv(logit)

################################################################################
## Alternative input options for covariance structure ####

### Transform icc and cac to random effects

icc_to_RandEff <- function(icc, cac=1, sigResid){
  if(any(c(icc,cac)<0,c(icc,cac)>1))
    stop("ICC and CAC must be between 0 and 1.")

  if (cac == 1) {
    gamma <- 0
    tau   <- sqrt(sigResid^2 * icc/(1-icc))
  }
  else {
    gamma <- sqrt(icc * sigResid^2 * (1-cac)/(1-icc))
    tau   <- sqrt(gamma^2 * cac / (1-cac))
  }
  return(list(gamma =gamma,
              tau   =tau))
}

RandEff_to_icc <- function(sigResid, tau, gamma=0){
  sigMarg <- sqrt(sigResid^2+tau^2+gamma^2)
  cac     <- tau^2 / (tau^2 + gamma^2)
  icc     <- (tau^2+gamma^2) / sigMarg^2
  return(list(icc=icc,
              cac=cac,
              sigMarg=sigMarg))
}


#' Correlation structure: transform random effects to alpha
#'
#' @param sigResid Residual standard deviation on individual leve
#' @param tau standard deviation of random cluster intercept
#' @param gamma standard deviation of random time effect
#' @param psi standard deviation of random subject specific intercept
#'
#' @return a list containing four named elements (possibly matrices):
#' `alpha0`, `alpha1`, `alpha2` specify a correlation structure and SigMarg
#' denotes the marginal standard deviation
#' @export
#'
#' @examples
#' RandEff_to_alpha012(sigResid=sqrt(11), tau=4, gamma=3, psi=2)
#'
#' ## The function is vectorised:
#' RandEff_to_alpha012(sigResid = matrix(c(0,1,2,3,4,5), 2, 3),
#'                     tau      = matrix(c(1,1,1,0,0,0), 2, 3),
#'                     gamma    = matrix(c(0,0,1,0,0,1), 2, 3),
#'                     psi      = matrix(c(0,1,1,0,0,1), 2, 3))

RandEff_to_alpha012 <- function(sigResid, tau, gamma, psi){
  SigMargSq <- (tau^2 + gamma^2 +psi^2 + sigResid^2)

  alpha0 <- (tau^2 + gamma^2) / SigMargSq
  alpha1 <- (tau^2) / SigMargSq
  alpha2 <- (tau^2 + psi^2)   / SigMargSq

  return(list(alpha0 = alpha0,
              alpha1 = alpha1,
              alpha2 = alpha2,
              SigMarg= sqrt(SigMargSq)))
}


#' Correlation structure: transform alpha to random effects
#'
#' @param alpha012 A vector or a list of length 3. Each list element must have
#' the same dimension.
#' @param sigResid Residual standard deviation on individual level. Either
#' residual sd or marginal sd needs to be specified.
#' @param sigMarg Marginal standard deviation on individual level. Either
#' residual sd or marginal sd needs to be specified.
#'
#' @return a list containing four named elements (possibly matrices):
#' random cluster intercept `tau`, random time effect `gamma`, random subject
#' intercept and residual standard deviation
#'
#' @export
#'
#' @examples
#' alpha012_to_RandEff(alpha012=c(.1,.1,.1), sigMarg=1)
#' alpha012_to_RandEff(alpha012=c(.1,.1,.1), sigResid=.9486833)
#'
#'## The function is vectorised:
#' alpha012_to_RandEff(alpha012=list(matrix(c(0,.1,.1,.2), 2, 2),
#'                                   matrix(c(0,0,.1,.2) , 2, 2),
#'                                   matrix(c(0,0,.2,.2) , 2, 2)),
#'                     sigMarg=1)
#'
alpha012_to_RandEff <- function(alpha012, sigResid=NULL, sigMarg=NULL){

  if(is.null(sigResid)==is.null(sigMarg))
    stop("Either `sigResid` or `sigMarg` must be declared (But not both).")

  a0 <- alpha012[[1]]
  a1 <- alpha012[[2]]
  a2 <- alpha012[[3]]

  sigMargSq <- if(is.null(sigResid))  sigMarg^2 else ( sigResid^2/(1-a0-a2+a1) )

  tau      <- sqrt( sigMargSq * a1 )
  gamma    <- sqrt( sigMargSq * (a0-a1) )
  psi      <- sqrt( sigMargSq * (a2-a1) )
  sigResid <- sqrt( sigMargSq * (1-a0-a2+a1) )

  return(list(tau      = tau,
              gamma    = gamma,
              psi      = psi,
              sigresid = sigResid))
}
