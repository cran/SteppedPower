
#' @title Construct a Block of the Covariance Matrix
#'
#' @description Constructs the covariance matrix
#' for multiple measurements of the same cluster.
#' This function is not designed to be used directly.
#'
#' @inheritParams construct_CovSubMat
#'
#' @return a block of a covariance matrix,
#' corresponding to intra-cluster covariance over time for one cluster
#'
#' @export
#'
#' @examples
#' construct_CovBlk(sigma=rep(2,5), tau=rep(1,5))


construct_CovBlk <- function(sigma,
                             tau,
                             eta   = NULL,
                             tauAR = NULL,
                             etaAR = NULL,
                             rho   = NULL){
  if(!(length(sigma)==length(tau)))
    stop ("In construct_CovBlk: sigma and tau must be of same length.")
  if(!is.null(eta) & !length(tau)==length(eta))
    stop ("In construct_CovBlk: sigma, tau and eta must be of same length.")

  timepoints <- length(sigma)

  tauMat <- if(is.null(tauAR)) tau %o% tau
            else toeplitz(tau^2 * tauAR ** c(0:(timepoints-1)))
  out    <- diag(sigma^2, timepoints) + tauMat

  if(!is.null(eta)) {
    etaMat <- if(is.null(etaAR)) eta %o% eta
              else toeplitz(eta^2 * etaAR ** c(0:(timepoints-1)))
    out <- out + etaMat
  }
  if(!is.null(rho)) {
    rhoVec <- (rho * eta * tau)
    out    <- out + outer(rhoVec, rhoVec, "+")
  }
  return(out)
}



#' @title Construct a Block of the Covariance Matrix
#'
#' @description Constructs the covariance matrix
#' for multiple measurements of the same cluster
#' if the same individuals are observed at all time periods.
#' This function is not designed to be used directly.
#'
#' @inheritParams construct_CovMat
#' @param N Number of individuals per cluster
#' @param sigma numeric (vector of length `timepoints`),
#' residual error
#' @param tau numeric (vector of length `timepoints`),
#'  standard deviation of random intercepts
#' @param eta numeric (vector of length `timepoints`),
#' standard deviation of random slope
#' @param gamma numeric (vector of length `timepoints`),
#' standard deviation of a random time effect.
#'
#' @return a block of a covariance matrix with two levels of clustering,
#' corresponding to intra-cluster covariance over time for one cluster
#'

construct_CovSubMat <- function(N,
                                timepoints,
                                sigma,
                                tau,
                                eta       = NULL,
                                tauAR     = NULL,
                                etaAR     = NULL,
                                rho       = NULL,
                                gamma     = 0,
                                trtMat    = NULL,
                                psi       = NULL,
                                INDIV_LVL = FALSE){

  ClBlk  <- construct_CovBlk(sigma = gamma,
                             tau   = tau,
                             eta   = eta,
                             tauAR = tauAR,
                             etaAR = etaAR,
                             rho   = rho)

  if(INDIV_LVL){
    SubMat <- construct_CovMat(SumCl      = N,
                               timepoints = timepoints,
                               sigma      = sigma,
                               tau        = psi)
    out <- matrix(1,N,N) %x% ClBlk + SubMat
  }else {
    SubBlk <- construct_CovBlk(sigma=sigma,
                               tau  =psi)
    out <- ClBlk + (1/N)*SubBlk
  }

  return(out)
}



#' @title Construct a Covariance Matrix
#'
#' constructs a (block diagonal) covariance matrix.
#' This function calls `construct_CovBlk`
#' (or `construct_CovSubMat` in case of repeated
#' observations of the same individuals) for each block.
#'
#' @inheritParams compute_wlsPower
#' @param timepoints numeric (scalar or vector), number of timepoints (periods).
#' If design is swd, timepoints defaults to length(Cl)+1.
#' Defaults to 1 for parallel designs.
#' @param SumCl total number of clusters
#' @param trtMat a matrix of dimension *#Cluster* x *timepoints* as produced by
#' the function `construct_trtMat`, indicating the cluster-periods that receive
#' interventional treatment. Defaults to NULL. If trtMat is given, the arguments
#' `SumCl` and `timepoints` are ignored (!).
#' @param CovBlk a matrix of dimension *timepoints* x *timepoints*.
#'
#' @return a covariance matrix
#' @export
#'
#' @examples
#'
#' ## Two clusters, three timepoints,
#' ## residual standard error sd=3, random slope sd=1.
#' construct_CovMat(SumCl=2, timepoints=3, sigma=3, tau=1)
#' ##
#' ##
#' ## ... with random slope as AR-1 process
#' construct_CovMat(SumCl=2, timepoints=3, sigma=3, tau=1, tauAR=.8)
#' ##
#' ##
#'
#'
#' ## ... with sigma and tau variing over time and between clusters:
#' construct_CovMat(SumCl=2,timepoints=3,
#'                  sigma=matrix(c(1,2,2,1,1,2),nrow=2, byrow=TRUE),
#'                  tau=matrix(c(.2,.1,.1,.2,.2,.1),nrow=2, byrow=TRUE),
#'                  N=c(3,4))




construct_CovMat <- function(SumCl      = NULL,
                             timepoints = NULL,
                             sigma,
                             tau,
                             eta        = NULL,
                             tauAR      = NULL,
                             etaAR      = NULL,
                             rho        = NULL,
                             gamma      = NULL,
                             trtMat     = NULL,
                             N          = NULL,
                             CovBlk     = NULL,
                             psi        = NULL,
                             INDIV_LVL  = FALSE){
  if(!is.null(CovBlk)){
    CovBlks <- rep(list(CovBlk),SumCl)
  } else {
    ## Checks ##
    if(is.null(SumCl)      & !is.null(trtMat)) SumCl      <- nrow(trtMat)
    if(is.null(timepoints) & !is.null(trtMat)) timepoints <- ncol(trtMat)
    timepoints  <- sum(timepoints)

    if(!is.null(rho) & is.null(eta))
      stop("In construct_CovMat: eta is needed if rho is not NULL")


    ## sigma ##
    lenS <- length(sigma)
    sigmaMat <- if(lenS %in% c(1,SumCl,SumCl*timepoints)){
      matrix(sigma, nrow=SumCl, ncol=timepoints)
    }else if(lenS==timepoints){
      matrix(sigma, nrow=SumCl, ncol=timepoints, byrow=TRUE)
    }else stop(paste('length of sigma is ', N,
                     '. This does not fit to given number of timepoints, ',
                     'which is ',timepoints,
                     ' or to the given number of clusters, which is ', SumCl))
    if(timepoints==SumCl & lenS==SumCl)
      warning("sigma is assumed to change between clusters. If you wanted sigma
              to change over time, please provide as matrix of dimension
              #Cluster x timepoints")

    ## N (if psi==NULL on cluster means, if psi!=NULL on individual level ##
    if(is.null(psi) & !INDIV_LVL){
      if(is.null(N)) N <- 1
      if(length(N) %in% c(1,SumCl,SumCl*timepoints)) {
        NMat <- matrix(N, nrow=SumCl, ncol=timepoints)
      }else stop(paste('length of cluster size vector N is ', N,
                       '. This does not fit to given number of clusters, ',
                       'which is ', SumCl,"\n"))

      ## N into sigma (aggregate on cluster means) ##
      sigmaMat <- sigmaMat / sqrt(NMat)
    }

    ## sigma transformed to list ##
    sigmaLst <- split(sigmaMat,row(sigmaMat))

    ## tau (input can be scalar, vector or matrix) ##
    tauMat <- matrix(tau, nrow=SumCl, ncol=timepoints)
    tauLst <- split(tauMat, row(tauMat))

    ## eta (input can be scalar or matrix, is passed as list of vectors) ##
    if(!is.null(eta)) {
      if(is.matrix(eta)){
        if(nrow(eta)==SumCl & ncol(eta)==timepoints)
          etaMat <- eta
        else stop("matrix dimensions of argument `eta` are ",
                  paste(dim(eta),collapse="x"), " but must be ",
                  SumCl,"x",timepoints)
      }else if(!is.null(trtMat) & length(eta)==1){
        etaMat <- trtMat * eta
      }else stop("If argument eta is a scalar, ",
                 "argument trtMat needs to be provided")
      etaLst <- split(etaMat, row(etaMat))
    }else
      etaLst <- vector("list", length=SumCl)

    ## rho (input must be scalar, is passed as scalar) ##
    if(!is.null(rho)) {
      rhoLst <- as.list(rep(rho,SumCl))
    }else
      rhoLst <- vector("list", length=SumCl)

    if(is.null(tauAR)) tauAR <- vector("list", length=SumCl)
    if(is.null(etaAR)) etaAR <- vector("list", length=SumCl)

    if(is.null(psi) & !INDIV_LVL){
      CovBlks <- mapply(construct_CovBlk,
                        sigma = sigmaLst,
                        tau   = tauLst,
                        eta   = etaLst,
                        tauAR = tauAR,
                        etaAR = etaAR,
                        rho   = rhoLst,
                        SIMPLIFY = FALSE)
    }else{
      NMat <- matrix(N, nrow=SumCl, ncol=1)
      NLst <- split(NMat, row(NMat))

      ## gamma ##
      gammaMat <- matrix(ifelse(is.null(gamma),0,gamma),
                         nrow=SumCl, ncol=timepoints)
      gammaLst <- split(gammaMat, row(gammaMat))

      psiMat <- matrix(ifelse(is.null(psi),0,psi), nrow=SumCl, ncol=timepoints)
      psiLst <- split(psiMat, row(psiMat))

      CovBlks <- mapply(construct_CovSubMat,
                        sigma = sigmaLst,
                        tau   = tauLst,
                        eta   = etaLst,
                        rho   = rhoLst,
                        N     = NLst,
                        gamma = gammaLst,
                        psi   = psiLst,
                        tauAR =tauAR,
                        etaAR =etaAR,
                        MoreArgs = list(timepoints = timepoints,
                                        INDIV_LVL  = INDIV_LVL),
                        SIMPLIFY = FALSE)
    }
  }
  CovMat <- Matrix::bdiag(CovBlks)
  if(!is.null(gamma) & is.null(psi)) {
    diag(CovMat) <- Matrix::diag(CovMat) + gamma^2
  }


  return(CovMat)
}
