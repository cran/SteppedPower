
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
  dists <- adist(Input, Options,
                 costs = c(insertions = 1, deletions = 1, substitutions = 1),
                 ignore.case = TRUE)
  best <- which.min(dists)
  ## Reject if the closest option is still too far away
  rel_dist <- dists[best] / nchar(Options[best])
  if (rel_dist > 0.5)
    stop("'", Input, "' does not match any of ",
         paste(Options, collapse = ", "),
         ". Did you misspell it?")
  Options[best]
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

fbdiag <- function(lmat) {
  ## fast block diagonal sparse matrix construction.
  ## no checks at all. Isanely unsave, for internal use only.
  ## assumes symmetrical blocks of equal length.

  l <- length(lmat)         ## number of blocks
  N <- nrow(lmat[[1]])      ## block size N x N
  s  <- rep(seq(0,(l-1)*N,by=N), each=N*(N+1)/2) ## create block shift vector

  i0 <- rep(1:N, times=N:1) ## row indices of upper triangle of first block
  i  <- rep(i0,l) + s

  j0 <- unlist(lapply(seq_len(N), function(k) seq(k, N))) ## column indices
  j  <- rep(j0,l) + s

  idx <- lower.tri(lmat[[1]], diag = TRUE) ## lower, bc unlist works columnwise
  x   <- unlist(lapply(lmat, "[", idx)) ## values in tiangles of blocks

  sparseMatrix(i, j, x=x, symmetric=TRUE, check=FALSE)
}


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


## Alternative input options for covariance structure ####


#' Transform ICC, CAC, IAC into random effects
#'
#' This function transforms a covariance structure specified by
#' intracluster correlation (ICC), cluster autocorrelation (CAC)
#' and individual autocorrelation (IAC) as in Hooper et al. (2016)
#' into the random effects notation. The latter specification type
#' is used by the workhorse function of this package, [`glsPower()`].
#'
#'
#' @param icc intracluster correlation
#' @param cac cluster autocorrelation
#' @param iac individual autocorrelation
#' @param sigMarg Marginal standard deviation
#' @param sigResid Residual standard deviation on individual level
#'
#' @return a list containing five named elements (possibly vectors or matrices):
#'  - random cluster intercept `tau`,
#'  - random time effect `gamma`,
#'  - random subject intercept `psi`,
#'  - residual standard deviation `sigResid` and
#'  - marginal standard deviation `sigMarg`.
#'
#' @details
#' The formulae used are \cr
#'  \deqn{ \sigma^2_{resid} = \sigma^2_{marg} (1-ICC) (1-IAC) \\[1ex]
#'         \tau^2 = \sigma^2_{marg} \cdot ICC \cdot CAC \\[1ex]
#'         \gamma^2 = \sigma^2_{marg} \cdot ICC \cdot (1-CAC) \\[1ex]
#'         \psi^2 =  \sigma^2_{marg} \cdot IAC \cdot (1-ICC) }
#'
#' Note that this does not allow to define the standard deviation of
#' a random treatment effect `eta`, but you can manually define it in the call to
#' `glsPower()`, see examples. \cr
#'
#' The **CAC** is sometimes interpreted in two different ways. *Two-period decay*
#' (the more common interpretation) allows correlation within the same period to
#' differ form those across periods. *Discrete time decay* (much less common interpretation),
#' on the other hand, assumes a consistent rate of correlation decay over time.
#' The former introduces an additional random term (i.e. an cluster-period specific
#' random intercept, with standard deviation \eqn{\gamma}), whereas the latter simply defines a fixed autocorrelation structure
#' that dictates how correlations diminish over time. Therefore, if **CAC** is specified,
#' is is interpreted as *two-period decay* and a (non-zero) value for gamma is returned. To define *discrete time decay*, please use
#' the `AR` option in the main function `glsPower()` instead.
#'
#'
#' @references Hooper, R., Teerenstra, S., de Hoop, E., & Eldridge, S. (2016). \cr
#' *Sample size calculation for stepped wedge and other longitudinal cluster randomised trials.* \cr
#' Statistics in medicine, 35(26), 4718-4728. DOI: 10.1002/sim.7028
#' @seealso [RandEff_to_icc()], [RandEff_to_alpha012()] or [alpha012_to_RandEff()]
#' @export
#' @examples
#' ## The function can be applied to vectors
#' tmp <- icc_to_RandEff(icc=c(0.01,0.005,0.001),
#'                       cac=1,
#'                       iac=.001,
#'                       sigMarg=sqrt(.5*(1-.5)) )
#' tmp
#'
#' ## This can then be inserted into `glsPower()` to calculate power for a
#' ## specific setting (albeit not as vectors)
#' ## Possibly with an additional random treatment effect `eta`.
#' glsPower(Cl=rep(2,4),
#'          N = 15,
#'          mu0=.5, mu1=.3,
#'          sigma = tmp$sigResid[1],
#'          tau   = tmp$tau[1],
#'          gamma = tmp$gamma[1],
#'          psi   = tmp$psi[1],
#'          eta   = 0.03 )
#'
#'
#'
icc_to_RandEff <- function(icc,
                           cac=1,
                           iac=0,
                           sigMarg=NULL,
                           sigResid=NULL){

  if(any(c(icc,cac,iac)<0, c(icc,cac,iac)>1))
    stop("ICC, CAC and IAC must be between 0 and 1.")

  if(is.null(sigResid)==is.null(sigMarg))
    stop("Either `sigResid` or `sigMarg` must be declared (But not both).")

  if (!is.null(sigMarg)) {
    sigResid <- sqrt(sigMarg^2 * (1-icc)*(1-iac) )
  } else if (!is.null(sigResid)) {
    sigMarg  <- sqrt(sigResid^2 / ((1-icc)*(1-iac)) )
  }
  ssq <- sigMarg^2

  tau      <- sqrt(ssq * icc*cac )
  gamma    <- sqrt(ssq * icc*(1-cac) )
  psi      <- sqrt(ssq * iac*(1-icc) )

  return(list(sigMarg  = sigMarg,
              sigResid = sigResid,
              gamma    = gamma,
              tau      = tau,
              psi      = psi))
}

#' Transform random effects into ICC, CAC, IAC
#'
#' This function transforms standard deviations of random effects into
#' intracluster correlation (ICC), cluster autocorrelation (CAC) and
#' individual autocorrelation (IAC).
#' It reproduces the formulae in Hooper et al. (2016).
#'
#' @param sigResid Residual standard deviation on individual level
#' @param tau standard deviation of random cluster intercept
#' @param gamma standard deviation of random time effect
#' @param psi standard deviation of random subject specific intercept
#'
#' @return a list containing four named elements (possibly vectors or matrices):
#'  - `icc` intracluster correlation
#'  - `cac` cluster autocorrelation
#'  - `iac` individual autocorrelation
#'  - `sigResid` Residual standard deviation on individual level
#'
#' @details
#' The formulae used are \cr
#'       \deqn{
#'        \sigma^2_{marg} = \sigma^2_{resid} + \tau^2 + \gamma^2 + \psi^2 \\[1ex]
#'        ICC = \frac{\tau^2 + \gamma^2}{\sigma^2_{marg}} \\[1ex]
#'        CAC = \frac{\tau^2}{\tau^2 + \gamma^2} \\[1ex]
#'        IAC = \frac{\psi^2}{\sigma^2_{resid}+\psi^2}
#'       }
#'
#'
#' @references Hooper, R., Teerenstra, S., de Hoop, E., & Eldridge, S. (2016). \cr
#' *Sample size calculation for stepped wedge and other longitudinal cluster randomised trials.* \cr
#' Statistics in medicine, 35(26), 4718-4728. DOI: 10.1002/sim.7028
#' @seealso [icc_to_RandEff()], [RandEff_to_alpha012()] or [alpha012_to_RandEff()]
#' @export
#'
RandEff_to_icc <- function(sigResid,
                           tau,
                           gamma=0,
                           psi=0){

  sigMarg <- sqrt(sigResid^2+tau^2+gamma^2+psi^2)
  icc     <- (tau^2+gamma^2) / sigMarg^2
  cac     <- tau^2 / (tau^2 + gamma^2)
  iac     <- psi^2 / (sigResid^2 + psi^2)

  return(list(icc=icc,
              cac=cac,
              iac=iac,
              sigMarg=sigMarg))
}


#' Transform random effects to alpha
#'
#' Transforms standard deviations of random effects into a
#' correlation structure specified by \eqn{\alpha_0, \alpha_1, \alpha_2}
#' as defined in Li et al. (2018).
#'
#' @param sigResid Residual standard deviation on individual level
#' @param tau standard deviation of random cluster intercept
#' @param gamma standard deviation of random time effect
#' @param psi standard deviation of random subject specific intercept
#'
#' @return a list containing four named elements (possibly vectors or matrices):
#' `alpha0`, `alpha1`, `alpha2` specify a correlation structure. `SigMarg`
#' denotes the marginal standard deviation
#' @seealso [RandEff_to_icc()], [icc_to_RandEff()] or [alpha012_to_RandEff()]
#'
#' @references Li, F., Turner, E. L., & Preisser, J. S. (2018). \cr
#' *Sample size determination for GEE analyses of stepped wedge cluster randomized trials.* \cr
#' Biometrics, 74(4), 1450-1458. DOI: 10.1111/biom.12918
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
  sigMargSq <- (tau^2 + gamma^2 +psi^2 + sigResid^2)

  alpha0 <- (tau^2 + gamma^2) / sigMargSq
  alpha1 <- (tau^2) / sigMargSq
  alpha2 <- (tau^2 + psi^2)   / sigMargSq

  return(list(alpha0 = alpha0,
              alpha1 = alpha1,
              alpha2 = alpha2,
              sigMarg= sqrt(sigMargSq)))
}


#' Transform alpha to random effects
#'
#' Transforms a correlation structure specified by \eqn{\alpha_0, \alpha_1, \alpha_2}
#' as defined in Li et al. (2018) into standard deviations of random effects.
#'
#' @param alpha012 A vector or a list of length 3. Each list element must have
#' the same dimension(s).
#' @param sigResid Residual standard deviation on individual level. Either
#' residual sd or marginal sd needs to be specified.
#' @param sigMarg Marginal standard deviation on individual level. Either
#' residual sd or marginal sd needs to be specified.
#'
#' @return a list containing four named elements (possibly vectors or matrices):
#'  - random cluster intercept `tau`,
#'  - random time effect `gamma`,
#'  - random subject intercept `psi`,
#'  - residual standard deviation `sigResid`.
#'
#' @seealso [RandEff_to_icc()], [icc_to_RandEff()] or [RandEff_to_alpha012()]
#'
#' @references Li, F., Turner, E. L., & Preisser, J. S. (2018). \cr
#' *Sample size determination for GEE analyses of stepped wedge cluster randomized trials.* \cr
#' Biometrics, 74(4), 1450-1458. DOI: 10.1111/biom.12918
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
              sigResid = sigResid))
}
