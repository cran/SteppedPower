#'@title
#'Compute power
#'
#' @description
#' This is the main function of the SteppedPower package.
#' It calls the constructor functions for the design matrix and
#' covariance matrix, and then calculates the variance of the
#' intervention effect estimator. The latter is then used
#' to compute the power of a Wald test of a (given) intervention effect.
#'
#'
#'
#' @param Cl integer (vector), number of clusters per sequence group (in SWD),
#' or number in control and intervention (in parallel designs)
#' @param timepoints numeric (scalar or vector), number of timepoints (periods).
#' If design is swd, timepoints defaults to length(Cl)+1.
#' Defaults to 1 for parallel designs.
#' @param DesMat matrix of dimension ... , if supplied,
#' `timepoints`,`Cl`,`trtDelay` are ignored.
#' @param trtDelay numeric (possibly vector), value(s)
#' between 0 and 1 specifying
#' the intervention effect in the first (second ... ) intervention phase
#' @param incomplete integer, either a vector (only for SWD) or a matrix.
#' A vector defines the number of periods before and after the switch from
#' control to intervention that are observed. A matrix consists of 1's for
#' observed clusterperiods and 0's for unobserved clusterperiods.
#' @param timeAdjust character, specifies adjustment for time periods.
#' One of the following: "factor", "linear", "none", "periodic".
#' Defaults to "factor".
#' @param dsntype character, defines the type of design. Options are "SWD",
#' "parallel" and "parallel_baseline", defaults to "SWD".
#' @param mu0 numeric (scalar), mean under control
#' @param mu1 numeric (scalar), mean under treatment
#' @param marginal_mu logical. Only relevant for non-gaussian outcome.
#' Indicates whether mu0 and mu1 are to be interpreted as marginal prevalence
#' under control  and under treatment, respectively, or whether they denote
#' the prevalence conditional on random effects being 0
#' (It defaults to the latter).
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param eta numeric (scalar or matrix), standard deviation of random slopes.
#' If `eta` is given as scalar, `trtMat` is needed as well.
#' @param tauAR numeric (scalar), value between 0 and 1. Defaults to NULL.
#' If `tauAR` is not NULL, the random intercept `tau` is AR1-correlated.
#' *Currently not compatible with `rho`!=0 !*
#' @param rho numeric (scalar), correlation of `tau` and `eta`
#' @param gamma numeric (scalar), random time effect
#' @param psi numeric (scalar), random subject specific intercept.
#' Leads to a closed cohort setting
#' @param alpha_0_1_2 numeric vector of length 3, that consists of
#' alpha_0, alpha_1 and alpha_2. This is an alternative way to define
#' the correlation structure, following Li et al. (2018).
#' @param N numeric, number of individuals per cluster. Either a scalar, vector
#' of length #Clusters or a matrix of dimension #Clusters x timepoints.
#' Defaults to 'rep(1,sum(Cl))' if not passed.
#' @param family character, distribution family. One of "gaussian", "binomial".
#' Defaults to "gaussian"
#' @param Power numeric, a specified target power.
#' If supplied, the minimal `N` is returned.
#' @param N_range numeric, vector specifying the lower and upper bound for `N`,
#' ignored if `Power` is NULL.
#' @param sig.level numeric (scalar), significance level, defaults to 0.05
#' @param dfAdjust character, one of the following: "none","between-within",
#' "containment", "residual".
#' @param verbose integer, how much information should the function return?
#' @param period numeric (scalar)
#' @param CovMat numeric, a positive-semidefinite matrix with
#' (#Clusters \eqn{\cdot} timepoints) rows and columns. If `CovMat` is given,
#' `sigma`, `tau`, `eta`, `rho` and `psi` are ignored.
#' @param INDIV_LVL logical, should the computation be conducted on an
#' individual level? This leads to longer run time and is
#' mainly for diagnostic purposes.
#'
#' @details
#' Let \eqn{\theta:= \mu_1-\mu_0} the treatment effect under investigation.
#' The variance of the treatment effect estimator \eqn{\hat\theta} can then be
#' estimated via weighted least squares (see also vignette 'Getting Started').
#'
#'
#'
#'
#'
#' @return
#' The return depends on the `verbose` parameter.
#' If `verbose`=0, only the power is returned
#' If `verbose`=1 (the default), a list containing power and the
#' parameters of the specific setting is returned.
#' If requested (by `verbose`=2) this list also contains relevant matrices.
#'
#' @export
#'
#' @examples
#' ## See also vignette for more examples
#' ##
#' ##
#' ## stepped wedge design with 5 Clusters in 5 sequences,
#' ## residual standard deviation 2,
#' ## cluster effect sd = 0.33, and 10 individuals per cluster.
#' ## Further, let the mean under the null and alternative hypothesis 0 and 1,
#' ## respectively.
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10)
#' ##
#' ##
#' ## ... with auto-regressive cluster effect `tauAR=0.7`.
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, tauAR=0.7, N=10)
#' ##
#' ##
#' ## ... with varying cluster size
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=c(12,8,10,9,14))
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33,
#'               N=matrix(c(12,8,10,9,14,
#'                          11,8,10,9,13,
#'                          11,7,11,8,12,
#'                          10,7,10,8,11,
#'                           9,7, 9,7,11,
#'                           9,6, 8,7,11),5,6))
#' ##
#' ##
#' ## ... with random treatment effect (with standard deviation 0.2),
#' ## which is correlated with the cluster effect with `rho`=0.25.
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, eta=.2, rho=.25, N=10)
#' ##
#' ##
#' ## ... with missing observations (a.k.a. incomplete stepped wedge design)
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10, incomplete=3)
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,
#'              incomplete=matrix(c(1,1,1,0,0,
#'                                  1,1,1,1,0,
#'                                  1,1,1,1,1,
#'                                  1,1,1,1,1,
#'                                  0,1,1,1,1,
#'                                  0,0,1,1,1),5,6))
#'## -> the same.
#'##
#'## ... with two levels of clustering. This arises if the patients are
#'## observed over the whole  study period
#'## (often referred to as closed cohort design) or if subclusters exist
#'## (such as wards within clinics). For
#'mod_aggr  <- wlsPower(mu0=0, mu1=1, Cl=rep(1,5),
#'                           sigma=2, tau=0.33, psi=.25,
#'                           N=10, incomplete=3, verbose=2)
#'mod_indiv <- wlsPower(mu0=0, mu1=1, Cl=rep(1,5),
#'                           sigma=2, tau=0.33, psi=.25,
#'                           N=10, incomplete=3, verbose=2, INDIV_LVL=TRUE)
#'mod_aggr
#'mod_indiv
#'## Compare covariance matrices of first cluster
#'mod_aggr$CovarianceMatrix[1:6,1:6] ; mod_indiv$CovarianceMatrix[1:60,1:60]
#'##
#'##
#'## stepped wedge design with 5 Clusters in 5 sequences, residual sd = 2,
#'## cluster effect sd = 0.33. How many Individuals are needed to achieve a
#'## power of 80% ?
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, Power=.8)
#'##
#'## ... How many are needed if we have a closed cohort design with a random
#'## individuum effect of .7?
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, psi=.7, Power=.8)
#'##
#'##
#'## longitudinal parallel design, with 5 time periods, 3 clusters in treatment
#'## and control arm each.
#' wlsPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
#'               dsntype="parallel", timepoints=5)
#'##
#'##
#'##
#'## ... with one baseline period and four parallel periods
#' wlsPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
#'               dsntype="parallel_baseline", timepoints=c(1,4))
#'##
#'##
#'##
#'## cross-over design with two timepoints before and two after the switch
#' wlsPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
#'               dsntype="crossover", timepoints=c(2,2))
#'##
#'##
#'##
#'## stepped wedge design with 32 Individuals in 8 sequences, binomial outcome,
#'## 50% incidence under control, 25% incidence under interventional treatment.
#'## cluster effect sd = 0.5 (ICC of 1/3 under control),
#'## every individual is its own cluster.
#'## ... with incidences defined conditional on cluster effect=0
#'wlsPower(mu0=0.5, mu1=0.25, Cl=rep(4,8), tau=0.5, N=1,
#'              family="binomial")
#'##
#'##
#'## ... with  marginally defined incidences
#' wlsPower(mu0=0.5, mu1=0.25, Cl=rep(4,8), tau=0.5, N=1,
#'               family="binomial", marginal_mu=TRUE)
#'



wlsPower <- function( Cl            = NULL,
                      timepoints    = NULL,
                      DesMat        = NULL,
                      trtDelay      = NULL,
                      incomplete    = NULL,
                      timeAdjust    = "factor",
                      period        = NULL,
                      dsntype       = "SWD",
                      mu0,
                      mu1,
                      marginal_mu   = FALSE,
                      sigma         = 1,
                      tau           = NULL,
                      eta           = NULL,
                      tauAR         = NULL,
                      rho           = NULL,
                      gamma         = NULL,
                      psi           = NULL,
                      alpha_0_1_2   = NULL,
                      CovMat        = NULL,
                      N             = NULL,
                      Power         = NULL,
                      family        = "gaussian",
                      N_range       = c(1,1000),
                      sig.level     = 0.05,
                      dfAdjust      = "none",
                      INDIV_LVL     = FALSE,
                      verbose       = 1){
  ## CHECKS #####
  if(!is.null(N) & !is.null(Power))
    stop("Both target power and individuals per cluster not NULL.")

  if(!all(sapply(c(tau,eta,rho,gamma,tauAR),is.null))){
    if(any(c(tau,eta,gamma)<0))
      stop("tau, eta and gamma must be >=0")
    if(!is.null(tauAR)){
      if(is.null(tau)) stop("If tauAR is supplied, tau is needed as well.")
      if(tauAR<0 | tauAR>1) stop("tauAR must be between 0 and 1.")
    }
    if(!is.null(rho)){
      if(is.null(eta) | is.null(tau))
        stop("If the correlation rho between random intercept and",
             " slope is not 0, a random slope must be provided.")
      if( (-1)>rho | rho>1 )
        stop("Correlation rho must be between -1 and 1")
    }
    if(is.null(tau)) tau <- 0 ## needed for computational reasons. Add warning?
  }else if (FALSE) {
    ## clause for alternative input options (icc & cac or alpha0-alpha2)
  }else{
    tau <- 0
    if(is.null(CovMat))
      warning("Random cluster effect tau and random treatment effect eta",
              " are assumed to be 0, i.e. the observations across clusters are",
              " assumed to be i.i.d. Declare tau=0 to supress this warning.")
  }

  if(!is.null(psi) & is.null(Power)){
    if(is.null(N))
      stop("If the standard deviation `psi` is not null, N is needed.")
    if(is.matrix(N)){
      N <- N[,1]
      warning("If psi is not NULL, the number of individuals per cluster must",
              "not change over time. Only the first column of N is considered.")
    }
  }

  ## Match string inputs ####
  ### dsntype
  dsntypeOptions <- c("SWD","parallel","parallel_baseline","crossover")
  tmpdsntype     <- dsntypeOptions[which.min(adist(dsntype,
                                                   dsntypeOptions,
                                                   costs=c(insertions    = 1,
                                                           deletions     = 100,
                                                           substitutions = 100),
                                                   ignore.case=TRUE))]
  if(dsntype != tmpdsntype) {
    message("Assumes ", tmpdsntype, " design")
    dsntype <- tmpdsntype
  }

  ###

  ## DesMat #####
  if(is.null(DesMat)){
    if(!is.null(timepoints) & !is.null(trtDelay)) {
      if(length(trtDelay)>max(timepoints)) {
        stop("The length of vector trtDelay must be",
             "less or equal to timepoints.")
    }}
    DesMat    <- construct_DesMat(Cl         = Cl,
                                  trtDelay   = trtDelay,
                                  dsntype    = dsntype,
                                  timepoints = timepoints,
                                  timeAdjust = timeAdjust,
                                  period     = period,
                                  N          = if(INDIV_LVL) N,
                                  INDIV_LVL  = INDIV_LVL )
  }else{
    if(!all(sapply(list(Cl, timepoints, trtDelay, period),is.null)))
      warning("If argument DesMat is provided, Cl, timepoints, trtDelay,",
              "timeAdjust, period and dsntype are ignored.")

    if(inherits(DesMat,"matrix") & !inherits(DesMat,"DesMat")){
      DesMat <- construct_DesMat(trtmatrix=DesMat)
    }else if(!inherits(DesMat,"DesMat"))
      stop("In wlsPower: Cannot interpret input for DesMat. ",
           "It must be either an object of class DesMat or a matrix")
    dsntype <- DesMat$dsntype
  }

  ## temporary variables #####
  timepoints <- DesMat$timepoints
  lenCl      <- length(DesMat$Cl)
  SumCl      <- sum(DesMat$Cl)


  ## incomplete designs #####
  if(!is.null(incomplete) & is.null(CovMat)){

    if(is.vector(incomplete) & dsntype=="SWD"){
      if(incomplete>timepoints) {
        incomplete <- timepoints
        warning("Argument `incomplete` must be less or equal to the number of",
                "timepoints. `incomplete` is set to ", timepoints )
      }
      Toep <- toeplitz(c(rep(1,incomplete),rep(Inf,lenCl-incomplete)))
      lastCols <- (timepoints-lenCl+1):timepoints

      IM <- matrix(1,lenCl,timepoints)
      IM[lower.tri(IM)]                       <- Toep[lower.tri(Toep)]
      IM[,lastCols][upper.tri(IM[,lastCols])] <- Toep[upper.tri(Toep)]

      IM <- IM[rep(seq_len(lenCl),DesMat$Cl),]

    }else if(is.matrix(incomplete)){
      if(!nrow(incomplete) %in% c(lenCl,SumCl) | ncol(incomplete)!=timepoints)
        stop("matrix dimensions of argument `incomplete` are ",
             paste(dim(incomplete),collapse="x"), " but must be ",
             paste(dim(DesMat$trtMat),collapse="x"), " or ",
             paste(dim(unique(DesMat$trtMat)),collapse="x"))
      IM <- incomplete
      IM[which(IM==0)] <- Inf
      if(nrow(incomplete)==lenCl) IM <- IM[rep(seq_len(lenCl),DesMat$Cl),]
    }
    sigma <- matrix(sigma, nrow=SumCl, ncol=timepoints,
                    byrow=ifelse(length(sigma)!=timepoints,TRUE,FALSE)) * IM
  }

  if(family =="binomial"){

    if(marginal_mu){

      mu0 <-muCond_to_muMarg(muCond=mu0, tauLin=tau)
      mu1 <-muCond_to_muMarg(muCond=mu1, tauLin=tau)
      print(paste("mu0=",round(mu0,5),", mu1=",round(mu1,5),"."))

    }

    sig0  <- sqrt(mu0*(1-mu0))
    sig1  <- sqrt(mu1*(1-mu1))
    ## for delayed trt effect only approximate sigma
    sigma <- matrix(sig0,
                    nrow=SumCl,
                    ncol=timepoints) + DesMat$trtMat * (sig1-sig0)

    OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
    print(paste("The assumed odds ratio is",round(OR,4))) ## user information
  }

  EffSize <- mu1-mu0
  if(marginal_mu) print(paste("The (raw) effect is",round(EffSize,5)))

  ## calculate samplesize (if needed, i.e. if power is not NULL ) #####
  if(!is.null(Power)){
    if(Power<0 | Power>1) stop("Power needs to be between 0 and 1.")
    N_opt <- tryCatch(ceiling(
              uniroot(function(N){Power - compute_wlsPower(DesMat    = DesMat,
                                                           EffSize   = EffSize,
                                                           sigma     = sigma,
                                                           tau       = tau,
                                                           eta       = eta,
                                                           tauAR     = tauAR,
                                                           rho       = rho,
                                                           gamma     = gamma,
                                                           psi       = psi,
                                                           N         = N,
                                                           dfAdjust  = dfAdjust,
                                                           sig.level = sig.level,
                                                           CovMat    = CovMat,
                                                           INDIV_LVL = INDIV_LVL,
                                                           verbose   = 0)},
                interval=N_range)$root),
              error=function(cond){
                message(paste0("Maximal N yields power below ",Power,
                               ". Increase argument N_range."))
                return(N_range[2])
              })
    N <- N_opt
  }
  ## calculate Power #####
  out <- compute_wlsPower(DesMat    = DesMat,
                          EffSize   = EffSize,
                          sigma     = sigma,
                          tau       = tau,
                          eta       = eta,
                          tauAR     = tauAR,
                          rho       = rho,
                          gamma     = gamma,
                          psi       = psi,
                          N         = N,
                          dfAdjust  = dfAdjust,
                          sig.level = sig.level,
                          CovMat    = CovMat,
                          INDIV_LVL = INDIV_LVL,
                          verbose   = verbose)
  if(!is.null(Power)) out$N_opt <- N_opt

  if(verbose>0) {
    out$Params <- append(out$Params,
                         list(mu0=mu0,
                              mu1=mu1))
    class(out) <- "wlsPower"
  }

  return(out)
}

#' @title Compute Power via weighted least squares
#'
#' @description This function calls `construct_DesMat` and `construct_CovMat` to
#' construct the design and covariance matrix, repectively. These matrices are
#' used to calculate the variance of the treatment effect estimator which is
#' then used to calculate the power to detect the assumed treatment effect.
#'
#' @inheritParams wlsPower
#' @param DesMat  list, containing a matrix, the design matrix,
#' numeric timepoints, numeric total number of Clusters
#' @param EffSize raw effect, i.e. difference between mean under control and
#' mean under intervention
#' @param etaAR numeric (scalar), value between 0 and 1. Defaults to NULL.
#' If `etaAR` is not NULL, the random slope `eta` is AR1-correlated.
#' *Currently not compatible with `rho`!=0 !*
#'
#' @return
#' The return depends on the `verbose` parameter.
#' If `verbose`=0, only the power is returned
#' If `verbose`=1 (the default), a list containing power and the
#' parameters of the specific setting is returned.
#' If requested (by `verbose`=2) this list also contains relevant matrices.
#'
#' @export

compute_wlsPower <- function(DesMat,
                             EffSize,
                             sigma,
                             tau        = 0,
                             eta        = NULL,
                             tauAR      = NULL,
                             etaAR      = NULL,
                             rho        = NULL,
                             gamma      = NULL,
                             psi        = NULL,
                             N          = NULL,
                             CovMat     = NULL,
                             dfAdjust   = "none",
                             sig.level  = .05,
                             INDIV_LVL  = FALSE,
                             verbose    = 1){
  dsnmatrix  <- DesMat$dsnmatrix
  timepoints <- DesMat$timepoints
  SumCl      <- sum(DesMat$Cl)
  SumSubCl   <- sum(DesMat$N)
  trtMat     <- DesMat$trtMat

  ## Checks ####
  if(!is.null(CovMat) & sum(sapply(c(sigma, tau, eta, rho, gamma, N),
                                   is.null)) > 0 )
    warning("If argument CovMat is provided, sigma, tau, eta, rho, gamma and N",
            "are ignored.")

  ## get covariance matrix #####
  if(is.null(CovMat))
    CovMat   <- construct_CovMat(SumCl      = SumCl,
                                 timepoints = timepoints,
                                 sigma      = sigma,
                                 tau        = tau,
                                 eta        = eta,
                                 tauAR      = tauAR,
                                 etaAR      = etaAR,
                                 rho        = rho,
                                 gamma      = gamma,
                                 psi        = psi,
                                 trtMat     = trtMat,
                                 N          = N,
                                 INDIV_LVL  = INDIV_LVL)

  ## matrices for power calculation #####
  tmpmat <- t(dsnmatrix) %*% Matrix::chol2inv(Matrix::chol(CovMat))
  VarMat <- Matrix::solve(tmpmat %*% dsnmatrix)
  if(verbose==2) ProjMat <- matrix((VarMat %*% tmpmat)[1,],
                                nrow=ifelse(INDIV_LVL,SumSubCl,SumCl),
                                byrow=TRUE)

  ## ddf for power calculation #####
  df <- switch(dfAdjust,
               "none"           = Inf,
               "between-within" = SumCl - rankMatrix(dsnmatrix),
               "containment"    = dim(dsnmatrix)[1] - SumCl,
               "residual"       = dim(dsnmatrix)[1] - rankMatrix(dsnmatrix))
  if(df<3){
    warning(dfAdjust,"-method not applicable. No DDF adjustment used.")
    df <- Inf }

  Pwr <- tTestPwr(d=EffSize, se=sqrt(VarMat[1,1]), df=df, sig.level=sig.level)
  if(verbose==0){
    out <- c(Pwr)
  } else {
    out <- list(Power  =Pwr,
                Params =list(N         = N,
                             sigma     = sigma,
                             tau       = tau,
                             eta       = eta,
                             tauAR     = tauAR,
                             etaAR     = etaAR,
                             rho       = rho,
                             gamma     = gamma,
                             psi       = psi,
                             denomDF   = df,
                             dfAdjust  = dfAdjust,
                             sig.level = sig.level))
  }
  if(verbose==2)
    out <- append(out,
                  list(ProjMatrix       = ProjMat,
                       DesignMatrix     = DesMat,
                       CovarianceMatrix = CovMat))
  return(out)
}

#' @title print.wlsPower
#'
#' @param x object of class wlsPower
#' @param ... Arguments to be passed to methods
#'
#' @method print wlsPower
#'
#' @return Messages, containing information about (at least) power and
#' significance level
#'
#' @export
#'
#'
print.wlsPower <- function(x, ...){
  message("Power                                = ", round(x$Power,4))
  if(x$Params$dfAdjust!="none"){
    message("ddf adjustment                       = ", x$Params$dfAdjust,"\n",
            "Denominator degrees of freedom       = ", x$Params$denomDF)
  }
  message("Significance level (two sided)       = ", x$Params$sig.level)
  if("N_opt" %in% names(x))
  message("Needed N per cluster per period      = ", x$N_opt,"\n")
}



#' @title plot an object of class `wlsPower`
#'
#' @description Plot a matrix that visualises the influence of each cluster for
#' each timepoint
#'
#' @param x object of class wlsPower
#' @param ... Arguments to be passed to methods
#'
#' @method plot wlsPower
#'
#' @return a plotly html widget
#'
#' @export
#'
plot.wlsPower <- function(x,...){
  if(!"ProjMatrix" %in% names(x))
    stop("Please rerun wlsPower() with `verbose=2` ")
  wgt <- x$ProjMatrix
  mx <- max(abs(wgt))
  sumCl <- dim(wgt)[1]
  timep <- dim(wgt)[2]

  subp <- suppressWarnings(subplot(
    plot_ly(data=data.frame(time   = seq_len(dim(wgt)[2]),
                            weight = colSums(abs(wgt))),
            type="bar", x=~time, y=~weight, color=I("grey")) %>%
      layout(yaxis=list(title="Sum|weights|"),
             xaxis=list(title="", showticklabels=FALSE))
    ,
    plotly_empty(type="scatter",mode="marker")
    ,
    plot_ly(x=seq_len(timep), y=seq_len(sumCl), z=wgt, type="heatmap",
            colors=grDevices::colorRamp(c("steelblue","white","firebrick")),
            xgap=.3, ygap=.3, name=" ",
            hovertemplate="Time: %{x}\nCluster: %{y}\nWeight: %{z}") %>%
      colorbar(len=1,limits=c(-mx,mx)) %>%
      layout(xaxis=list(title="time"),
             yaxis=list(title="cluster", autorange="reversed"))
    ,
    plot_ly(data=data.frame(cluster=seq_len(dim(wgt)[1]),
                            weight=rowSums(abs(wgt))),
            type="bar", orientation="h",
            y=~cluster, x=~weight, color=I("grey")) %>%
      layout(xaxis=list(title="Sum|weights|"),
             yaxis=list(title="", showticklabels=FALSE, autorange="reversed"))
    ,
    nrows=2, heights=c(.2,.8), widths=c(.8,.2), titleX=TRUE, titleY=TRUE
  ) %>% layout(showlegend=FALSE)
  )
  return(subp)
}

