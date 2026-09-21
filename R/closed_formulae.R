

## formula from Kasza (open cohort)
## includes formulae of HuHu, of Hooper&Girling&Hemming and of Li (2018) as
## special cases

#' Closed formula for treatment variance in open cohort settings
#'
#' @description
#' From Kasza et al "Sample size and power calculations for open cohort
#' longitudinal cluster rondomized trials" 2020
#'
#' @param trtMat a matrix trtMat to define treatment allocation,
#' where rows and columns correspond to cluster and timepoints, respectively
#' @param tau numeric, standard deviation of random intercepts
#' @param gamma numeric, random time effect
#' @param psi numeric, random subject specific intercept.
#' @param sigma numeric, residual error on subject level.
#' @param N numeric, number of individuals per cluster.
#' @param chi Attrition factor
#'
#' @return numeric, variance of the estimator for treatment effect
#' @export
#'
#' @examples
#' ##  test setting, from Hussey&Hughes 2007  ####
#' trtMat <- construct_DesMat(c(6,6,6,6))$trtMat
#' tau <- .025 ; sigma <- sqrt(.041*.959) ; N <- 100 ;
#' gamma <- 0.01 ; psi <- .1 ; chi <- .7
#'
#' tmp <- VarClosed_Kasza(trtMat, tau=tau, sigma=sigma, gamma=0, psi=0, N=N, chi=0)
#' tTestPwr((.05-.032), sqrt(tmp), df = Inf)
#' glsPower(Cl = rep(6,4), N=N, mu0=.05, mu1=.032, verbose=0,
#'         sigma=sigma, gamma=0, tau=tau, psi=0)
#'
#' tmp <- VarClosed_Kasza(trtMat, tau=tau, sigma=sigma, gamma=gamma, psi=psi, N=N, chi=0)
#' tTestPwr((.05-.032), sqrt(tmp), df = Inf)
#' glsPower(Cl = rep(6,4), N=N, mu0=.05, mu1=.032, verbose=0,
#'         sigma=sigma, gamma=gamma, tau=tau, psi=psi)
#'
#' tmp <- VarClosed_Kasza(trtMat, tau=tau, sigma=sigma, gamma=gamma, psi=psi, N=N, chi=1)
#' tTestPwr((.05-.032), sqrt(tmp), df = Inf)
#' glsPower(Cl = rep(6,4), N=N, mu0=.05, mu1=.032, verbose=0,
#'          sigma=sigma, gamma=sqrt(gamma^2+psi^2/N), tau=tau, psi=0)
#'
#' tmp <- VarClosed_Kasza(trtMat, tau=tau, sigma=sigma, gamma=gamma, psi=psi, N=N, chi=chi)
#' tTestPwr((.05-.032), sqrt(tmp), df = Inf)
#' glsPower(Cl = rep(6,4), N=N, mu0=.05, mu1=.032, verbose=0,
#'          sigma=sigma, gamma=sqrt(gamma^2+chi*psi^2/N), tau=tau, psi=sqrt(1-chi)*psi)


VarClosed_Kasza <- function(trtMat, tau, gamma=0, psi=0, sigma, N, chi){

  XX <- trtMat
  II <- dim(XX)[1]   ## number of clusters
  JJ <- dim(XX)[2]   ## number of time periods

  UU <- sum(XX) ## number of cluster-periods under intervention
  VV <- sum((rowSums(XX))^2)
  WW <- sum((colSums(XX))^2)

  sigMargSq <- tau^2 + gamma^2 + psi^2 + sigma^2
  iccA <- (tau^2+psi^2) / sigMargSq
  iccB <-  tau^2 / sigMargSq
  iccW <- (tau^2+gamma^2) / sigMargSq
  tmp  <- chi * (iccA - iccB) ## chi is attrition rate  0:=cohort, 1:=cross-sectional

  lamda1 <- 1 + (N-1)*(iccW-iccB) - iccA + tmp
  lamda2 <- 1 + (N-1)*iccW + (JJ-1)*(N-1)*iccB + (JJ-1)*iccA - tmp*(JJ-1)

  a <- (sigMargSq/N) * II*JJ * lamda1 * lamda2
  b <- (UU^2 + II*JJ*UU - JJ*WW - II*VV) * lamda2 - (UU^2-II*VV)*lamda1

  return(varTheta=a/b)
}


## function from F. Li for proportional decay ####
## from Li (2020) "Design and analysis considerations for cohort swcrt with decay"

#' Closed formula for treatment variance, with proportional decay
#'
#' @description
#' From Li et al "Design and analysis considerations for cohort stepped wedge
#' cluster randomized trials with a decay correlation structure"
#'
#' @inheritParams VarClosed_Kasza
#' @param AR numeric (scalar), It defines the AR(1)-correlation of random effects.
#'
#' @return numeric, variance of the estimator for treatment effect
#' @export
#'
#' @examples
#' ##  test setting, from Hussey&Hughes 2007  ####
#' trtMat <- construct_DesMat(c(6,6,6,6))$trtMat
#' tau <- .025 ; N <- 100 ; psi <- .1 ; AR <- .6

#' tmp <- VarClosed_Li(trtMat, tau=tau, psi=psi, N=N, AR=AR)
#' tTestPwr((.05-.032), se=sqrt(tmp), Inf)
#' glsPower(Cl=rep(6,4), mu0=.05, mu1=.032, AR=AR,
#'          tau=tau, N=N, sigma=0, psi=psi, verbose=0)


VarClosed_Li <- function(trtMat, tau, psi, N, AR){

  XX <- trtMat
  II <- dim(XX)[1]   ## number of clusters
  JJ <- dim(XX)[2]   ## number of time periods

  sigMargSq <- tau^2 + psi^2
  rho_w     <- tau^2/sigMargSq

  UU <- sum(XX) ## number of cluster-periods under intervention
  PP <- sum(XX[,1:(JJ-1)] * XX[,2:(JJ)])
  WW <- sum((colSums(XX))^2)
  QQ <- sum(colSums(XX[,1:(JJ-1)]) * colSums(XX[,2:(JJ)]))

  a <- sigMargSq * II /N * (1-AR^2) *(1 + (N-1)*rho_w)
  b <- ((II*UU - WW)*(1+AR^2) - 2*(II*PP - QQ)*AR)
  return(varTheta= a/b)
}




