#'@title
#' Compute power via weighted least squares
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
#' @param DesMat Either an object of class `DesMat` or a matrix indicating the
#' treatment status for each cluster at each timepoint. If supplied,
#' `timepoints`,`Cl`,`trtDelay` are ignored.
#' @param trtDelay numeric (possibly vector), value(s)
#' between 0 and 1 specifying the proportion of intervention effect
#' in the first (second ... ) intervention phase.
#' @param incomplete integer, either a scalar (only for SWD) or a matrix.
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
#' (It defaults to the latter). *(experimental!)*
#' @param sigma numeric, residual error of cluster means if no N given.
#' @param tau numeric, standard deviation of random intercepts
#' @param eta numeric (scalar or matrix), standard deviation of random slopes.
#' If `eta` is given as scalar, `trtMat` is needed as well.
#' @param AR numeric, vector containing up to three values, each between 0 and 1.
#' Defaults to NULL. It defines the AR(1)-correlation of random effects.
#' The first element corresponds to the cluster intercept, the second to the
#' treatment effect and the third to subject specific intercept.
#' If only one element is provided, autocorrelation of all random effects is
#' assumed to be the same.
#' *Currently not compatible with `rho`!=0 !*
#' @param rho numeric (scalar), correlation of `tau` and `eta`. The default is no correlation.
#' @param gamma numeric (scalar), random time effect
#' @param psi numeric (scalar), random subject specific intercept.
#' Leads to a closed cohort setting
#' @param alpha_0_1_2 numeric vector or list of length 2 or 3, that consists of
#' alpha_0, alpha_1 and alpha_2. Can be used instead of random effects to define
#' the correlation structure, following Li et al. (2018). When omitting alpha_2,
#' this describes a cross-sectional design, where alpha_0 and alpha_1 define
#' the intracluster correlation and cluster autocorrelation, respectively - as
#' defined by Hooper et al. (2016).
#' @param N numeric, number of individuals per cluster. Either a scalar, vector
#' of length #Clusters or a matrix of dimension #Clusters x timepoints.
#' Defaults to 1 if not passed.
#' @param family character, distribution family. One of "gaussian", "binomial".
#' Defaults to "gaussian"
#' @param power numeric, a specified target power.
#' If supplied, the minimal `N` is returned.
#' @param N_range numeric, vector specifying the lower and upper bound for `N`,
#' ignored if `power` is NULL.
#' @param sig.level numeric (scalar), significance level, defaults to 0.05
#' @param dfAdjust character, one of the following: "none","between-within",
#' "containment", "residual".
#' @param verbose integer, how much information should the function return?
#' See also under `Value`.
#' @param period numeric (scalar)
#' @param CovMat numeric, a positive-semidefinite matrix with
#' (#Clusters \eqn{\cdot} timepoints) rows and columns. If `CovMat` is given,
#' `sigma`, `tau`, `eta`, `rho`, `gamma` and `psi` as well as `alpha_0_1_2`
#' must be NULL.
#' @param INDIV_LVL logical, should the computation be conducted on an
#' individual level? This leads to longer run time and is
#' mainly for diagnostic purposes.
#' @param INFO_CONTENT logical, should the information content of cluster cells be
#' computed? The default is `TRUE` for designs with less or equal than 2500
#' cluster cells, otherwise `FALSE`. Ignored if `verbose=0`.
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
#' If `verbose`=1 (the default), a list containing power, projection matrix and the
#' parameters of the specific setting is returned.
#' If explicitly requested (by `verbose`=2) this list also contains
#' the `DesMat`-object and the covariance matrix.
#'
#' If INFO_CONTENT= TRUE, the returned list contains a named list with four elements:
#' `Cells` is explicit computation of the information content in each cell;
#' `Cluster` is the information content of entire clusters;
#' `time` is thie information content of entire time periods and
#' `Closed` is a formula-based computation the information content in each cell,
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
#' ## ... with auto-regressive cluster effect `AR=0.7`.
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, AR=0.7, N=10)
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
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, power=.8)
#'##
#'## ... How many are needed if we have a closed cohort design with a random
#'## individuum effect of .7?
#' wlsPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, psi=.7, power=.8)
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
#'## ... with  marginally defined proportions
#' wlsPower(mu0=0.5, mu1=0.25, Cl=rep(4,8), tau=0.5, N=1,
#'               family="binomial", marginal_mu=TRUE)
#'
#'##
#'##
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
                      sigma         = NULL,
                      tau           = NULL,
                      eta           = NULL,
                      AR            = NULL,
                      rho           = NULL,
                      gamma         = NULL,
                      psi           = NULL,
                      alpha_0_1_2   = NULL,
                      CovMat        = NULL,
                      N             = NULL,
                      power         = NULL,
                      family        = "gaussian",
                      N_range       = c(1,1000),
                      sig.level     = 0.05,
                      dfAdjust      = "none",
                      INDIV_LVL     = FALSE,
                      INFO_CONTENT  = NULL,
                      verbose       = 1){
  ## Match string inputs ####
  ### dsntype
  dsntypeOptions <- c("SWD","parallel","parallel_baseline","crossover")
  tmpdsntype     <- choose_character_Input(dsntypeOptions, dsntype)
  if(dsntype != tmpdsntype) {
    message("Assumes ", tmpdsntype, " design")
    dsntype <- tmpdsntype
  }
  ### family
  familyOptions <- c("gaussian", "binomial")
  tmpfamily     <- choose_character_Input(familyOptions, family)
  if(family != tmpfamily) {
    message("Assumes ", tmpfamily, "distribution")
    family <- tmpfamily
  }

  ## CHECKS #####
  if(!is.null(N) & !is.null(power))
    stop("Both target power and individuals per cluster not NULL. ",
         "Either N or power must be NULL.")

  if(is.null(sigma) & family=="gaussian" & is.null(CovMat))
    stop("For gaussian distribution, sigma must be provided.")

  if(!is.null(sigma) & family=="binomial")
    warning("Argument sigma is not used for binomial distribution.")

  ## Check covariance information #####
  UseRandEff <- !all(sapply(c(tau,eta,rho,gamma,AR), is.null))
  Usealpha   <- !is.null(alpha_0_1_2)
  UseCovMat  <- !is.null(CovMat)
  UsedOptions <- sum(UseRandEff, Usealpha, UseCovMat)

  if (UsedOptions==0) UseRandEff <- TRUE
  if (UsedOptions>=2)
    stop("There are three different alternatives to specify the covaricance, ",
         "structure, \nyou must use exactly one.\nPlease specify EITHER \n",
         "  - random effects: tau, eta, rho, gamma, psi   OR \n",
         "  - alpha_0_1_2                                 OR \n",
         "  - CovMat")

  if (UseRandEff) {
    if(any(c(tau,eta,gamma, psi)<0))
      stop("tau, eta, gamma and psi must be >=0")
    if(!is.null(AR)){
      if(is.null(tau)) stop("If AR is supplied, tau is needed as well.")
      if(is.null(eta) & length(AR)==2)
        stop("If AR has length 2, eta must be supplied.")
      if(is.null(psi) & length(AR)==3)
        stop("If AR has length 3, psi must be supplied.")
      if(min(AR)<0 | max(AR)>1) stop("AR must be between 0 and 1.")
      AR <- rep(AR, length.out=3)
    }
    if(!is.null(rho)){
      if(is.null(eta) | is.null(tau))
        stop("If the correlation rho between random intercept and",
             " slope is not 0, a random slope must be provided.")
      if( (-1)>rho | rho>1 )
        stop("Correlation rho must be between -1 and 1")
    }
    if(!is.null(psi) & is.null(power)){
      if(is.null(N))
        stop("If the standard deviation `psi` is not null, N is needed.")
      if(is.matrix(N)){
        N <- N[,1]
        warning("If psi is not NULL, the number of individuals per cluster must",
                "not change over time. Only the first column of N is considered.")
      }
    }
  }else if (Usealpha) {
    if(length(alpha_0_1_2)==2){
      alpha_0_1_2 <- append(alpha_0_1_2, alpha_0_1_2[[2]])
      message("Since length of alpha_0_1_2 is 2, a cross-sectional design is",
              "assumed. Hence, alpha2 is set to alpha1.")
    }
    if(alpha_0_1_2[[2]] > alpha_0_1_2[[1]] + alpha_0_1_2[[3]])
      stop("Correlation matrix defined by alpha_0_1_2 is not positve definite.",
           "\nThe following must hold:   alpha1 < alpha0 + alpha2")
  }else if (UseCovMat){
    # if(min(eigen(CovMat)$values) > 0)
      # stop("Covariance matrix is not positive definite")
  }


  ## construct Design Matrix #####
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
                                  incomplete = incomplete,
                                  N          = if(INDIV_LVL) N,
                                  INDIV_LVL  = INDIV_LVL )
  }else{
    if(inherits(DesMat, "DesMat")) {
      if(!all(sapply(list(Cl, timepoints ,trtDelay,period,incomplete),is.null)))      ## timeAdjust defaults to factor (!)
        warning("If input to argument DesMat inherits class `DesMat`, \n",
                "Cl, timepoints, trtDelay, incomplete,",
                "timeAdjust, period and dsntype are ignored.")
    } else if(inherits(DesMat,"matrix") & !inherits(DesMat,"DesMat")){
      DesMat <- construct_DesMat(trtmatrix  = DesMat,
                                 timeAdjust = timeAdjust,
                                 period     = period,
                                 incomplete = incomplete,
                                 N          = if(INDIV_LVL) N,
                                 INDIV_LVL  = INDIV_LVL)
      if(!all(sapply(list(Cl, timepoints, trtDelay, dsntype), is.null)))
        warning("If input to argument DesMat is of class `matrix`, \n",
                "Cl, timepoints, trtDelay, dsntype are ignored.")
    }else
      stop("In wlsPower: Cannot interpret input for DesMat. ",
           "It must be either an object of class DesMat or a matrix")
    dsntype <- DesMat$dsntype
  }

  ## declare temporary variables #####
  timepoints <- DesMat$timepoints
  lenCl      <- length(DesMat$Cl)
  sumCl      <- sum(DesMat$Cl)

  ## default for INFO_CONTENT ####
  if(is.null(INFO_CONTENT)){
    INFO_CONTENT <- ifelse(length(DesMat$trtMat)<=2500 & sumCl>2 &
                           verbose>0 & !INDIV_LVL,                   TRUE,FALSE)
  }else if(INFO_CONTENT & sumCl<=2) {
    INFO_CONTENT <- FALSE
    warning("Information Content not (yet) implemented for only two clusters.",
            "INFO_CONTENT set to FALSE.")
  }

  ## distribution family ####
  if(family =="gaussian"){
    if(Usealpha){
      tmp   <- alpha012_to_RandEff(alpha012=alpha_0_1_2, sigResid=sigma)
      tau   <- tmp$tau
      gamma <- tmp$gamma
      psi   <- tmp$psi
    }
  } else if(family =="binomial"){

    if(marginal_mu){
      if(!UseRandEff)
        stop("marginal_mu currently only implemented for random effects")
      mu0 <-muCond_to_muMarg(muCond=mu0, tauLin=tau)
      mu1 <-muCond_to_muMarg(muCond=mu1, tauLin=tau)
      print(paste("mu0=",round(mu0,5),", mu1=",round(mu1,5),"."))
    }

    muMat   <- matrix(mu0, sumCl, timepoints) + DesMat$trtMat*(mu1-mu0)
    sigma   <- sqrt(muMat * (1-muMat))

    if (verbose>0) {
      OR <- (mu1*(1-mu0))/(mu0*(1-mu1))
      print(paste("The assumed odds ratio is",round(OR,4))) ## user information
    }

    if(Usealpha){
      tmp   <- alpha012_to_RandEff(alpha012=alpha_0_1_2, sigResid=sigma)
      tau   <- tmp$tau
      gamma <- tmp$gamma
      psi   <- tmp$psi
    }
  }

  EffSize <- mu1-mu0

  if(marginal_mu & verbose>0)
    print(paste("The (raw) effect is",round(EffSize,5)))


  ## incomplete designs #####
  if(!is.null(DesMat$incompMat) & is.null(CovMat)){
    IM <- DesMat$incompMat
    IM[IM==0] <- Inf

    sigma <- matrix(sigma, nrow=sumCl, ncol=timepoints,
                    byrow=ifelse(length(sigma)!=timepoints,TRUE,FALSE)) * IM
  }

  ## calculate samplesize (if needed, i.e. if power is not NULL ) #####
  if(!is.null(power)){
    if(power<0 | power>1) stop("power needs to be between 0 and 1.")
    N_opt <- tryCatch(ceiling(
              uniroot(function(N){power - compute_wlsPower(DesMat    = DesMat,
                                                           EffSize   = EffSize,
                                                           sigma     = sigma,
                                                           tau       = tau,
                                                           eta       = eta,
                                                           AR        = AR,
                                                           rho       = rho,
                                                           gamma     = gamma,
                                                           psi       = psi,
                                                           N         = N,
                                                           dfAdjust  = dfAdjust,
                                                           sig.level = sig.level,
                                                           CovMat    = CovMat,
                                                           INDIV_LVL = INDIV_LVL,
                                                           INFO_CONTENT = FALSE,
                                                           verbose   = 0)},
                interval=N_range)$root),
              error=function(cond){
                message(paste0("Maximal N yields power below ",power,
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
                          AR        = AR,
                          rho       = rho,
                          gamma     = gamma,
                          psi       = psi,
                          N         = N,
                          dfAdjust  = dfAdjust,
                          sig.level = sig.level,
                          CovMat    = CovMat,
                          INDIV_LVL = INDIV_LVL,
                          INFO_CONTENT = INFO_CONTENT,
                          verbose   = verbose)
  if(!is.null(power)) out$N_opt <- N_opt

  if(verbose>0) {
    out$Params <- append(out$Params,
                         list(mu0         = mu0,
                              mu1         = mu1,
                              family      = family,
                              alpha_0_1_2 = alpha_0_1_2))
    class(out) <- "wlsPower"
  }

  return(out)
}

#' @title Compute power via weighted least squares
#'
#' @description
#' This function is not intended to be used directly, but rather to be called
#' by `wlsPower` - the main function of this package.
#' It expects the design matrix as an input argument `DesMat` and
#' construct the covariance matrix (if not given as well). These matrices are
#' used to calculate the variance of the treatment effect estimator which is
#' then used to calculate the power to detect the assumed treatment effect.
#'
#' @inheritParams wlsPower
#' @param DesMat  object of class `DesMat`.
#' @param EffSize raw effect, i.e. difference between mean under control and
#' mean under intervention
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
                             AR         = NULL,
                             rho        = NULL,
                             gamma      = NULL,
                             psi        = NULL,
                             N          = NULL,
                             CovMat     = NULL,
                             dfAdjust   = "none",
                             sig.level  = .05,
                             INDIV_LVL  = FALSE,
                             INFO_CONTENT = FALSE,
                             verbose    = 1){
  dsn        <- DesMat$dsnmatrix
  tp         <- DesMat$timepoints
  sumCl      <- sum(DesMat$Cl)
  SumSubCl   <- sum(DesMat$N)
  trtMat     <- DesMat$trtMat

  ## get covariance matrix #####
  if(is.null(CovMat))
    CovMat   <- construct_CovMat(sumCl      = sumCl,
                                 timepoints = tp,
                                 sigma      = sigma,
                                 tau        = tau,
                                 eta        = eta,
                                 AR         = AR,
                                 rho        = rho,
                                 gamma      = gamma,
                                 psi        = psi,
                                 trtMat     = trtMat,
                                 N          = N,
                                 INDIV_LVL  = INDIV_LVL)

  ## matrices for power calculation #####
  dsncols <- dim(dsn)[2]

  XW  <- matrix( ((tdsn <- t(dsn)) %*% chol2inv(drop0(chol(CovMat))))@x, dsncols ) ## drop0 for incomplete designs
  Var <- spdinv( VarInv <- XW %*% dsn )
  if(verbose>0) ProjMat <- matrix((Var[1,] %*% XW),
                                  nrow=ifelse(INDIV_LVL,SumSubCl,sumCl),
                                  byrow=TRUE)

  ## Information content, if requested ####
  # currently computed twice, once explicitly, once using the specific formula
  # explicit computation will be removed in near future
  if(INFO_CONTENT){
    I <- 1:sumCl
    J <- 1:tp
    InfoContent <- list(Cells   = matrix(NA,sumCl,tp),
                        Cluster = numeric(sumCl),
                        time    = numeric(tp))
    tp_drop <- array(0,dim=c(dsncols,dsncols,tp))

    if(length(CovMat@x)<tp*tp*sumCl) {  ## ugly. Is needed to produce consistent sparse matrix indexing for tau=0 (and eta >=0)
      i <- rep(J,tp)      + (i_add <- rep(tp*(I-1),each=tp*tp) ) ## wouldnt be necessary with bdiag_m ...
      j <- rep(J,each=tp) +  i_add
      CovMat <- CovMat + sparseMatrix(i,j,x=0)
    }

    for(i in I){
      J_start  <- tp*(i-1)
      J_incomp <- if(is.null(DesMat$incompMat)) J else J[DesMat$incompMat[i,]==1]  ## no computation of empty cells

      Var_drop <-spdinv( VarInv - XW[,(J_start+J)] %*% submatrix(dsn,J_start+1,J_start+tp,1,dsncols) )
      InfoContent$Cluster[i] <- Var_drop[1,1]/Var[1,1]

      if(tp>1){
        for(j in J_incomp){
          J_drop <- (J_ <- J[-j]) + J_start
          XW_drop <- matrix(0, dsncols, tp)  ## add *real* zeros to remain consistent with sparse matrix indexing
          XW_drop[,J_] <- tdsn[,J_drop] %*%
            spdinv( matrix(CovMat@x[tp^2*(i-1) + 1:(tp)^2],tp) [J_,J_] )
          Var_drop <- spdinv( VarInv +
              (tp_drop_updt <- (XW_drop[,J]-XW[,J_start+J]) %*% dsn[J_start+J,]) )
          tp_drop[,,j]           <- tp_drop[,,j] + tp_drop_updt
          InfoContent$Cells[i,j] <- Var_drop[1,1]/Var[1,1]
        }
      }
    }
    if( tp>1 & sum(colSums(DesMat$trtMat)>0)>1 ){
      if(DesMat$timeAdjust=="factor"){ ## TODO: ADD WARNINGS !!
        for(j in J){
          Var_drop <- spdinv( (VarInv + tp_drop[,,j])[-(j+1),-(j+1)] )
          InfoContent$time[j] <- Var_drop[1,1]/Var[1,1]
        }
      }else {
        for(j in J){
          Var_drop <- spdinv( VarInv + tp_drop[,,j] )
          InfoContent$time[j] <- Var_drop[1,1]/Var[1,1]
        }
      }
    }

  ## Formula-based calculation of information content
    W  <- spdinv(as.matrix(CovMat))
    X2 <- dsn[,-1]
    x1 <- dsn[, 1]
    Q  <- W %*% X2 %*% spdinv(t(X2)%*%W%*%X2) %*%t(X2)%*%W
    WQ <- W - Q
    hh <- as.numeric( t(x1)%*%WQ / c(t(x1)%*%WQ%*%x1) )
    InfoContent$Closed <- matrix(1/(1 - hh^2*as.numeric(t(x1)%*%WQ%*%x1) / diag(WQ) ),
                                 sumCl,tp, byrow=TRUE)

    # ## Check consistency of methods
    # maxDiff <- 0
    # maxDiff <- max(abs(InfoContent$Cells - InfoContent$Closed))
    # if(maxDiff>1e-12)  warning("formula-based information content and explicit information content differ")
  }


  ## ddf for power calculation #####
  df <- switch(dfAdjust,
               "none"           = Inf,
               "between-within" = sumCl - rankMatrix(dsn),
               "containment"    = dim(dsn)[1] - sumCl,
               "residual"       = dim(dsn)[1] - rankMatrix(dsn))
  if(df<3){
    warning(dfAdjust,"-method not applicable. No DDF adjustment used.")
    df <- Inf }

  Pwr <- tTestPwr(d=EffSize, se=sqrt(Var[1,1]), df=df, sig.level=sig.level)
  if(verbose==0){
    out <- Pwr
  } else {
    out <- list(power  =Pwr,
                Params =list(Cl         = DesMat$Cl,
                             timeAdjust = DesMat$timeAdjust,
                             timepoints = DesMat$timepoints,
                             designtype = DesMat$dsntype,
                             trtDelay   = DesMat$trtDelay,
                             N          = N,
                             sigma      = sigma,
                             tau        = tau,
                             eta        = eta,
                             AR         = AR,
                             rho        = rho,
                             gamma      = gamma,
                             psi        = psi,
                             denomDF    = df,
                             dfAdjust   = dfAdjust,
                             sig.level  = sig.level),
                ProjMatrix = ProjMat)
  if(INFO_CONTENT){
    out <- append(out,
                  list(InformationContent= InfoContent))
  }
  if(verbose==2)
    out <- append(out,
                  list(DesignMatrix     = DesMat,
                       CovarianceMatrix = CovMat,
                       VarianceMatrix   = Var))
  return(out)
  }
}

#' @title Print an object of class `wlsPower`
#'
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
  message("Power                                = ", round(x$power,4))
  if(x$Params$dfAdjust!="none"){
    message("ddf adjustment                       = ", x$Params$dfAdjust,"\n",
            "Denominator degrees of freedom       = ", x$Params$denomDF)
  }
  message("Significance level (two sided)       = ", x$Params$sig.level)
  if("N_opt" %in% names(x))
  message("Needed N per cluster per period      = ", x$N_opt,"\n")
}



#' @title plot the information content of a wls object
#'
#' @inheritParams plot.wlsPower
#' @param IC a matrix with information content for each cluster at each time period
#'
#' @return a plotly object
#' @export
#'

plot_InfoContent <- function(IC,
                             annotations=NULL,
                             show_colorbar=TRUE,
                             marginal_plots=TRUE){

  if(is.null(annotations)){
    annotations <- ifelse(length(IC$Cells)<=1e2,TRUE,FALSE)
  }
  mx  <- max(IC$Cells)
  sumCl <- dim(IC$Cells)[1]
  timep <- dim(IC$Cells)[2]
  gaps  <- 20/sqrt(length(IC$Cells))

  dat=cbind(expand.grid(y=seq_len(sumCl),
                        x=seq_len(timep)),
            z=as.numeric(IC$Cells),
            zChar=as.character(round(IC$Cells,3)))
  dat$zChar[is.na(dat$zChar)] <- ""

  PLT <- plot_ly(data=dat, x=~x, y=~y,  z=~z,
          type="heatmap",
          colors=grDevices::colorRamp(c("white","gold","firebrick")),
          xgap=gaps, ygap=gaps, name=" ",
          showscale=show_colorbar,
          hovertemplate="Time: %{x}\nCluster: %{y}\nInfoContent: %{z:.6f}" ) %>%
    colorbar(len=1,limits=c(1-1e-8,mx)) %>%
    layout(xaxis=list(title="Time"),
           yaxis=list(title="Cluster", autorange="reversed"))
  if(annotations) PLT <- PLT %>% add_annotations(text=~zChar, showarrow=FALSE)


  if(marginal_plots){
    PLT <- subplot(
      plot_ly(data=data.frame(time   = seq_along(IC$time),
                              InfoC  = IC$time),
              type="bar", x=~time, y=~(InfoC-1), color=I("grey"),
              name=" ", base=1,
              hovertemplate="Time: %{x}\nInfoContent: %{y:.6f}") %>%
        layout(yaxis=list(title=""),
               xaxis=list(title="", showticklabels=FALSE))
      ,
      plotly_empty(type="scatter",mode="marker")
      ,
      PLT
      ,
      plot_ly(data=data.frame(cluster = seq_along(IC$Cluster),
                              InfoC   = IC$Cluster),
              type="bar", orientation="h",
              y=~cluster, x=~(InfoC-1), color=I("grey"),
              name=" ", base=1,
              hovertemplate="Cluster: %{y}\nInfoContent: %{x:.6f}") %>%
        layout(xaxis=list(title=""),
               yaxis=list(title="", showticklabels=FALSE, autorange="reversed"))
      ,
      nrows=2, heights=c(.2,.8), widths=c(.8,.2), titleX=TRUE, titleY=TRUE
    ) %>% layout(showlegend=FALSE)
  }
  PLT
}


#' @title plot cell contributions (weights) of a wls object
#'
#' @inheritParams plot.wlsPower
#'
#' @return a plotly html widget
#' @export
#'
plot_CellWeights <- function(x,
                             annotations=NULL,
                             show_colorbar=TRUE,
                             marginal_plots=TRUE){

  if(is.null(annotations)){
    annotations <- ifelse(length(x$ProjMatrix)<=1e2,TRUE,FALSE)
  }
  wgt <- x$ProjMatrix
  mx  <- max(abs(wgt))
  sumCl <- dim(wgt)[1]
  timep <- dim(wgt)[2]
  gaps  <- 20/sqrt(length(wgt))

  if(!is.null(x$DesignMatrix$incompMat))
    wgt[x$DesignMatrix$incompMat==0] <- NA

  dat <- cbind(expand.grid(y=seq_len(sumCl),
                           x=seq_len(timep)),
               wgt=as.numeric(wgt),
               wgtChar=as.character(round(wgt,3)) )
  dat$wgtChar[is.na(dat$wgtChar)] <- ""

  PLT <- plot_ly(data=dat, x=~x, y=~y, z=~wgt, type="heatmap",
                 colors=grDevices::colorRamp(c("steelblue","white","firebrick")),
                 xgap=gaps, ygap=gaps, name=" ",
                 showscale=show_colorbar,
                 hovertemplate="Time: %{x}\nCluster: %{y}\nWeight: %{z:.6f}") %>%
    colorbar(len=1,limits=c(-mx,mx)) %>%
    layout(xaxis=list(title="Time"),
           yaxis=list(title="Cluster", autorange="reversed"))
  if(annotations) PLT <- PLT %>% add_annotations(x=dat$x, y=dat$y,
                                                 text=dat$wgtChar, showarrow=FALSE)

  if(marginal_plots){
    PLT <- suppressWarnings(
      subplot(
        plot_ly(data=data.frame(time   = seq_len(timep),
                                weight = colSums(abs(wgt),na.rm=TRUE)),
                type="bar", x=~time, y=~weight, color=I("grey"),
                name=" ",
                hovertemplate="Time: %{x}\nWeight: %{y:.6f}") %>%
          layout(yaxis=list(title="Sum|weights|"),
                 xaxis=list(title="", showticklabels=FALSE))
        ,
        plotly_empty(type="scatter",mode="marker")
        ,
        PLT
        ,
        plot_ly(data=data.frame(cluster=seq_len(sumCl),
                                weight=rowSums(abs(wgt),na.rm=TRUE)),
                type="bar", orientation="h",
                y=~cluster, x=~weight, color=I("grey"),
                name=" ",
                hovertemplate="Cluster: %{y}\nAbsWeight: %{x:.6f}") %>%
          layout(xaxis=list(title="Sum|weights|"),
                 yaxis=list(title="", showticklabels=FALSE, autorange="reversed"))
        ,
        nrows=2, heights=c(.2,.8), widths=c(.8,.2), titleX=TRUE, titleY=TRUE
      ) %>% layout(showlegend=FALSE)
    )
  }
  PLT
}



#' @title plot an object of class `wlsPower`
#'
#' @description Up to four plots (selectable by `which`) that visualise:
#' the contribution of each cluster-period cell to the treatment effect estimator,
#' the information content of each cluster-period cell,
#' the treatment status for each cluster for each time point and
#' the covariance matrix. By default, only the first two plots are returned.
#'
#' @param x object of class wlsPower
#' @param which Specify a subset of the numbers `1:4` to select plots. The default is
#' `1:2` or `1`, depending on whether `x` contains the information content.
#' @param show_colorbar logical, should the colorbars be shown?
#' @param ... Arguments to be passed to methods
#' @param annotations logical, should the cell contributions be annotated in the Plot?
#' @param marginal_plots should the influence of whole periods, clusters also be plotted?
#'
#' @method plot wlsPower
#'
#' @return a list of plotly html widgets
#'
#' @export
#'
plot.wlsPower <- function(x, which=NULL, show_colorbar=NULL,
                          annotations=NULL, marginal_plots=TRUE,...){
  out <- list()
  if(is.null(which))
    which <- if("InformationContent" %in% names(x)) 1:2 else 1

  if (1 %in% which){
    out$WgtPlot <- plot_CellWeights(x, annotations=annotations,
                                show_colorbar=show_colorbar,
                                marginal_plots=marginal_plots)
  }

  if (2 %in% which){
    if(!("InformationContent" %in% names(x)) ) stop("Please rerun wlsPower() with INFO_CONTENT=TRUE")
    out$IMplot <- plot_InfoContent(x$InformationContent,
                                        annotations=annotations,
                                        show_colorbar=show_colorbar,
                                        marginal_plots=marginal_plots)
  }

  if (3 %in% which){
    if(!("DesignMatrix" %in% names(x)) ) stop("Please rerun wlsPower() with verbose=2")
    out$DMplot <- plot(x$DesignMatrix, show_colorbar=show_colorbar)
  }

  if (4 %in% which){
    if(!("CovarianceMatrix" %in% names(x)) ) stop("Please rerun wlsPower() with verbose=2")
    out$CMplot <- plot_CovMat(x$CovarianceMatrix, show_colorbar=show_colorbar)
  }

  return(out)
}

