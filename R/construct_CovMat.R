
#' @title Construct a Single Block of the Covariance Matrix
#'
#' @description Constructs the covariance matrix
#' for multiple measurements of the same cluster.
#' This function is usually called by `construct_CovMat` and is
#'  not designed to be used directly.
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
#'
#' construct_CovBlk(sigma=rep(2,5),
#'                 tau=rep(.5,5), eta=c(0,0,1,1,1),
#'                 AR=c(.5, 1))

construct_CovBlk <- function(sigma,
                             tau   = NULL,
                             eta   = NULL,
                             AR    = NULL,
                             rho   = NULL){
  if(!(length(tau) %in% c(0,length(tau)) ) )
    stop ("In construct_CovBlk: If tau is provided, ",
          "sigma and tau must be of same length.")
  if(!is.null(eta) & !length(tau)==length(eta))
    stop ("In construct_CovBlk: sigma, tau and eta must be of same length.")

  # AR         <- rep(AR, length.out=3) ## NEEDED ???
  timepoints <- length(sigma)

  out  <- diag(sigma^2, timepoints)

  if(!is.null(tau)){
    tauMat <- if(is.null(AR[[1]])) tau %o% tau
              else tau %o% tau * toeplitz(AR[[1]] ^ (0:(timepoints-1)))
    out <- out + tauMat
  }
  if(!is.null(eta)) {
    etaMat <- if(is.null(AR[[2]])) eta %o% eta
              else eta %o% eta * toeplitz(AR[[2]] ^ (0:(timepoints-1)))
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
                                AR        = NULL,
                                rho       = NULL,
                                gamma     = NULL,
                                psi       = NULL,
                                INDIV_LVL = FALSE){

  ClBlk  <- construct_CovBlk(sigma = vector("double", timepoints),
                             tau   = tau,
                             eta   = eta,
                             AR    = c(AR[[1]],AR[[2]]),
                             rho   = rho)
  if(!is.null(gamma)) diag(ClBlk) <- diag(ClBlk) + gamma^2


  if(INDIV_LVL){
    SubMat <- construct_CovMat(sumCl      = N,
                               timepoints = timepoints,
                               sigma      = sigma,
                               tau        = psi,
                               gamma      = NULL,
                               AR         = c(AR[[3]],AR[[2]]))
    out <- matrix(1,N,N) %x% ClBlk + SubMat
  }else {
    SubBlk <- construct_CovBlk(sigma = sigma,
                               tau   = psi,
                               AR    = c(AR[[3]],AR[[2]]))
    out <- ClBlk + (1/N)*SubBlk
  }

  return(out)
}



#' @title Construct a Covariance Matrix
#'
#' @description
#' constructs a (block diagonal) covariance matrix.
#' This function calls \code{\link{construct_CovBlk}}
#' (or \code{\link{construct_CovSubMat}} in case of repeated
#' observations of the same individuals) for each block.
#'
#' @inheritParams glsPower
#' @param timepoints numeric (scalar or vector), number of timepoints (periods).
#' If design is swd, timepoints defaults to length(Cl)+1.
#' Defaults to 1 for parallel designs.
#' @param sumCl total number of clusters
#' @param trtMat a matrix of dimension *#Cluster* x *timepoints* as produced by
#' the function \code{\link{construct_trtMat}}, indicating the cluster-periods that receive
#' interventional treatment. Defaults to NULL. If trtMat is given, the arguments
#' `sumCl` and `timepoints` are ignored (!).
#' @param CovBlk a matrix of dimension *timepoints* x *timepoints*.
#'
#' @return a covariance matrix
#' @export
#'
#' @examples
#'
#' ## Two clusters, three timepoints,
#' ## residual standard error sd=3, random slope sd=1.
#' construct_CovMat(sumCl=2, timepoints=3, sigma=3, tau=1)
#' ##
#' ##
#' ## ... with random slope as AR-1 process
#' construct_CovMat(sumCl=2, timepoints=3, sigma=3, tau=1, AR=.8)
#' ##
#' ##
#' ## ... with sigma and tau variing over time and between clusters:
#' construct_CovMat(sumCl=2,timepoints=3,
#'                  sigma=matrix(c(1,2,2,1,1,2),nrow=2, byrow=TRUE),
#'                  tau=matrix(c(.2,.1,.1,.2,.2,.1),nrow=2, byrow=TRUE),
#'                  N=c(3,4))




construct_CovMat <- function(sumCl      = NULL,
                             timepoints = NULL,
                             sigma,
                             tau,
                             eta        = NULL,
                             AR         = NULL,
                             rho        = NULL,
                             gamma      = NULL,
                             trtMat     = NULL,
                             N          = NULL,
                             CovBlk     = NULL,
                             psi        = NULL,
                             INDIV_LVL  = FALSE){
  cross_sectional <- ifelse(is.null(psi), TRUE, ifelse(psi==0, TRUE, FALSE))

  if(!is.null(CovBlk)){
    CovBlks <- rep(list(CovBlk),sumCl)
  } else {
    ## Checks ##
    if(is.null(sumCl)      & !is.null(trtMat)) sumCl      <- nrow(trtMat)
    if(is.null(timepoints) & !is.null(trtMat)) timepoints <- ncol(trtMat)
    timepoints  <- sum(timepoints)

    if(!is.null(rho) & is.null(eta))
      stop("In construct_CovMat: eta is needed if rho is not NULL")

    ## sigma ##
    sigmaMat <- input_to_Mat(sigma, sumCl, timepoints)

    ## N into sigma (aggregate on cluster means) ##
    if( cross_sectional & !INDIV_LVL){
      if(is.null(N)) N <- 1
      NMat     <- input_to_Mat(N, sumCl, timepoints)
      sigmaMat <- sigmaMat / sqrt(NMat)
    }
    ## sigma transformed to list ##
    sigmaLst <- split(sigmaMat,1:nrow(sigmaMat))

    ## tau (input can be scalar, vector or matrix) ##
    tauLst <- input_to_List(tau, sumCl, timepoints)

    ## eta (input can be scalar or matrix, is passed as list of vectors) ##
    if(!is.null(eta)) {
      if(is.matrix(eta)){
        if(nrow(eta)==sumCl & ncol(eta)==timepoints)
          etaMat <- eta
        else stop("matrix dimensions of argument `eta` are ",
                  paste(dim(eta),collapse="x"), " but must be ",
                  sumCl,"x",timepoints)
      }else if(!is.null(trtMat) & length(eta)==1){
        etaMat <- trtMat * eta
      }else stop("If argument eta is a scalar, ",
                 "argument trtMat needs to be provided")
      etaLst <- split(etaMat, 1:nrow(etaMat))
    }else
      etaLst <- vector("list", length=sumCl)

    ## rho (input must be scalar, is passed as scalar) ##
    ## AR  (input must be scalar, is passed as scalar) ##
    ## gamma (input can be scalar or matrix, is passed as list of vectors) ##

    if(cross_sectional & !INDIV_LVL){

      CovBlks <- mapply(construct_CovBlk,
                        sigma = sigmaLst,
                        tau   = tauLst,
                        eta   = etaLst,
                        MoreArgs = list(AR  = AR,
                                        rho = rho),
                        SIMPLIFY = FALSE)
      # CovMat       <- Matrix::bdiag(CovBlks)
      CovMat       <- fbdiag(CovBlks)

      if(!is.null(gamma)) {
        gammaMat     <- input_to_Mat(gamma, sumCl, timepoints)
        diag(CovMat) <- Matrix::diag(CovMat) + as.numeric(t(gammaMat))^2
      }
    }else{
      gammaLst <- input_to_List(gamma, sumCl, timepoints)
      psiLst   <- input_to_List(psi, sumCl, timepoints)
      NLst     <- input_to_List(N, sumCl, ifelse(INDIV_LVL, 1, timepoints))

      CovBlks <- mapply(construct_CovSubMat,
                        sigma = sigmaLst,
                        tau   = tauLst,
                        eta   = etaLst,
                        N     = NLst,
                        gamma = gammaLst,
                        psi   = psiLst,
                        MoreArgs = list(timepoints = timepoints,
                                        INDIV_LVL  = INDIV_LVL,
                                        AR         = AR,
                                        rho        = rho),
                        SIMPLIFY = FALSE)
      # CovMat       <- Matrix::bdiag(CovBlks)
      CovMat       <- fbdiag(CovBlks)
    }
  }
  return(CovMat)
}


#' @title Visualise a Covariance Matrix
#'
#' @description Currently not exported.
#'
#' @param CovMat A covariance matrix (possibly in sparse matrix notation)
#' @param show_colorbar logical, should the colorbar be shown?
#'
#' @return a plotly object
#'

plot_CovMat <- function(CovMat, show_colorbar=FALSE){

  CovMat    <- as.matrix(CovMat)
  CMcols    <- diag(CovMat!=Inf)

  seqLength <- 1:(dim(CovMat)[1])
  seqLength <- seqLength[CMcols]
  gaps <- 20/dim(CovMat)[1]

  ## Work-around for incomplete designs, .. is there a nicer way?
  # tmpaux <- colSums(CovMat)==Inf
  # CovMat[tmpaux,] <- 0
  # CovMat[,tmpaux] <- 0

  plot_ly(type="heatmap", colors=c("white","steelblue"),
          x=~seqLength, y=~seqLength, z=~CovMat[CMcols,CMcols],
          xgap=gaps, ygap=gaps,
          showscale=show_colorbar) %>%
    layout(xaxis=list(title="", visible=FALSE),
           yaxis=list(title="", visible=FALSE, autorange="reversed") ) %>%
    colorbar(len=1, title="")
}

