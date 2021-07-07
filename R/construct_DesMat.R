#' @title Construct the Design Matrix
#'
#' @description
#' Constructs the design matrix with one column for every (fixed)
#' parameter to be estimated and one row for every cluster for every timepoint.
#' This function calls `construct_trtMat` to construct a matrix that indicates
#' treatment status for each cluster at each timepoint.
#' This is then transformed into the first
#' column of the design matrix. `construct_CovMat` further calls
#' `construct_timeAdjust` to get the fixed effect(s) of the timepoints.
#'
#' Note: Unlike the usual notation, the treatment effect is in the first column
#' (for easier access by higher level functions).
#'
#' @inheritParams wlsPower
#' @param trtmatrix an optional user defined matrix
#' to define treatment allocation
#' @param timeBlk an optional user defined matrix that defines
#' the time adjustment in one cluster.
#' Is repeated for every cluster.
#'
#' @return an object of class DesMat
#' @export
#'
#' @examples
#' construct_DesMat(Cl=c(2,0,1))
#' construct_DesMat(Cl=c(2,0,1), N=c(1,3,2))
#'
#' ## manually defined time adjustment (same as above)
#' timeBlock <- matrix(c(1,0,0,0,
#'                       1,1,0,0,
#'                       1,0,1,0,
#'                       1,0,0,1), 4, byrow=TRUE)
#' construct_DesMat(Cl=c(2,0,1), timeBlk=timeBlock)
#'
construct_DesMat <- function(Cl          = NULL,
                             trtDelay    = NULL,
                             dsntype     = "SWD",
                             timepoints  = NULL,
                             timeAdjust  = "factor",
                             period      = NULL,
                             trtmatrix   = NULL,
                             timeBlk     = NULL,
                             N           = NULL,
                             INDIV_LVL   = FALSE){
  if(INDIV_LVL){
    if(length(N)==1){
      N <- rep(N,sum(Cl))
    }
    tmpCl <- sapply(split(N,rep(as.factor(seq_along(Cl)),Cl)),sum)
  }else {
    tmpCl <- Cl
  }

  ## TREATMENT MATRIX ####
  if(!is.null(trtmatrix)){

    if(inherits(trtmatrix,"matrix")){
      trtMat  <- trtmatrix
    }else if (inherits(trtMat, "list") & "swDsn" %in% names(trtmatrix)){
      trtMat <- trtmatrix$swDsn
    } else stop("trtmatrix must be a matrix. It is a ",class(trtMat))

    dsntype <- "userdefined"
    timepoints  <- ncol(trtMat)
    Cl          <- table(do.call(paste,split(trtMat,col(trtMat))))
    tmpCl       <- Cl
  }else{
    trtMat  <- construct_trtMat(Cl            =Cl,
                                trtDelay      =trtDelay,
                                dsntype       =dsntype,
                                timepoints    =timepoints)
    timepoints <- dim(trtMat)[2]  ## construct_trtMat has good heuristics for guessing
                                  ## number of timepoints (if not provided)
  }
  if(INDIV_LVL)  tmpTrtMat <- trtMat[rep(seq_len(sum(Cl)),N),] else
                 tmpTrtMat <- trtMat

  ## TIME ADJUSTMENT ####
  timeBlks <- construct_timeAdjust(Cl          =tmpCl,
                                  timepoints   =timepoints,
                                  timeAdjust   =timeAdjust,
                                  period       =period,
                                  timeBlk      =timeBlk)

  dsnmatrix <- cbind(trt=as.numeric(t(tmpTrtMat)), timeBlks)

  DesMat  <- list(dsnmatrix  = dsnmatrix,
                  timepoints = timepoints,
                  trtDelay   = trtDelay,
                  Cl         = Cl,
                  N          = if(INDIV_LVL) N else NULL,
                  dsntype    = dsntype,
                  timeAdjust = ifelse(is.null(timeBlk),
                                      timeAdjust,
                                      "userdefined"),
                  trtMat     = trtMat)
  class(DesMat) <- append(class(DesMat),"DesMat")

  return(DesMat)
}


## Methods for class DesMat

#'  print.DesMat
#'
#' @param x  An object of class `DesMat
#' @param ... Arguments to be passed to methods
#'
#' @method print DesMat
#'
#' @return Messages with information about the design.
#'
#' @export
#'
print.DesMat <- function(x, ...){

  dsn_out <- switch (x$dsntype,
                    "SWD"               = "stepped wedge" ,
                    "parallel"          = "parallel",
                    "parallel_baseline" = "parallel with baseline period(s)",
                    "userdefined"       = "userdefined")

  message("Timepoints                         = ", x$timepoints,"\n",
          "Number of clusters per seqence     = ", paste(x$Cl, collapse= ", "))
  if(!is.null(x$N)){
  message("Number of subclusters per cluster  = ", paste(x$N, collapse=", "))
  }
  message("Design type                        = ", dsn_out,"\n",
          "Time adjustment                    = ", x$timeAdjust, "\n",
          "Dimension of design matrix         = ", dim(x$dsnmatrix)[1]," x ",
                                                   dim(x$dsnmatrix)[2],"\n",
          "\nTreatment status (clusters x timepoints):")
  print(x$trtMat)
}



#' @title plot.DesMat
#'
#' @param x An object of class `DesMat`
#' @param show_colorbar logical, should the colorbar be shown?
#' @param ... Arguments to be passed to methods
#'
#' @method plot DesMat
#'
#' @return a plotly html widget, displaying the treatment status
#'
#' @export
#'
#' @examples
#' x <- construct_DesMat(C=c(2,2,2,0,2,2,2),.5)

plot.DesMat <- function(x, show_colorbar=FALSE, ...){
  trt <- x$trtMat
  plot_ly(type="heatmap",
          x=~(seq_len(dim(trt)[2])), y=~(seq_len(dim(trt)[1])),
          z=~trt, xgap=5, ygap=5, name=" ",
          showscale=show_colorbar,
          colors=grDevices::colorRamp(c("steelblue","lightgoldenrod1","firebrick")),
          hovertemplate="Time: %{x},   Cluster: %{y} \nTreatment Status: %{z}") %>%
    layout(xaxis = list(title="time", type="category"),
           yaxis = list(title="cluster",autorange="reversed",type="category")) %>%
    colorbar(len=1, title="")
}

#' @title Construct Treatment Matrix
#'
#' @description
#' Constructs a matrix of `#cluster` rows and `#timepoint` columns, indicating
#' treatment status in each cluster at each timepoint.
#'
#' @inheritParams construct_DesMat
#'
#' @return a matrix trtMat, where rows and columns correspond to cluster
#' and timepoints, respectively
#'
#' @export
#'
#'
#' @examples construct_trtMat(Cl=c(1,2,1), trtDelay=c(.2,.8), dsntype="SWD")
#'
#'
construct_trtMat <- function(Cl,
                             trtDelay,
                             dsntype,
                             timepoints=NULL){

  SumCl         <- sum(Cl)
  lenCl         <- length(Cl)

  if(dsntype=="SWD"){
    if(is.null(timepoints)) timepoints <- length(Cl) + 1
    trt    <- matrix(0,lenCl,timepoints)
    trt[upper.tri(trt)] <- 1
    if(!is.null(trtDelay)){
      for(i in seq_along(trtDelay)){
        diag(trt[,-(1:i)]) <- trtDelay[i]
      }
    }
  }else if(dsntype=="parallel"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(is.null(timepoints)){
      if(is.null(trtDelay)){
        timepoints <- 1 ;  warning("timepoints unspecified. Defaults to 1.")
      }else{
        timepoints <- length(trtDelay)+1
        message("timepoints unspecified. Defaults to ", length(trtDelay)+1,
                " (length of trtDelay plus 1)")
      }
    }
    trt     <- matrix(0,nrow=2,ncol=timepoints)
    trt[2,] <- c(trtDelay,rep(1,(timepoints-length(trtDelay))))
  }else if(dsntype=="parallel_baseline"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(length(timepoints)==1){
      timepoints01 <- c(1,timepoints-1)
      message(paste("assumes 1 baseline period and",timepoints-1,
                    "parallel period(s). If intended otherwise,
                    argument timepoints must have length two."))
    }else if(length(timepoints)==2){
      timepoints01 <- timepoints
      timepoints   <- sum(timepoints)
    }else if(is.null(timepoints)){
      timepoints01 <- c(1,length(trtDelay)+1)
      timepoints   <- sum(timepoints01)
      message("timepoints unspecified. Defaults to 1 baseline,",
              length(trtDelay)+1, " parallel period(s).")
    }
    trt     <- matrix(0,nrow=2,ncol=timepoints)
    trt[2,] <- c(rep(0,timepoints01[1]),
                trtDelay,rep(1,(timepoints01[2]-length(trtDelay))))
  }else if (dsntype=="crossover"){
    if(length(Cl)!=2) {stop("In construct_DesMat: Cl must be of length 2.")}
    if(length(timepoints)==1){
      if(timepoints==1) stop("crossover designs must consist of
                             at least 2 timepoints.")
      timepoints01 <- c(floor(timepoints/2),ceiling(timepoints/2))
      message(paste("assumes", floor(timepoints/2) ,"AB period(s) and",
                    ceiling(timepoints/2),"BA period(s). If intended otherwise,
                    argument timepoints must have length two."))
    }else if (length(timepoints)==2){
      lenTp        <- timepoints
      timepoints   <- sum(timepoints)
    }else if(is.null(timepoints)){
      len          <- length(trtDelay)+1
      lenTp        <- c(len,len)
      timepoints   <- sum(lenTp)
      message("timepoints unspecified. Defaults to ", len,
              " AB period(s), and ", len, " BA period(s).")
    }
    trt   <- matrix(0, nrow=2, ncol=timepoints)
    vecAB <- c(trtDelay,rep(1,lenTp[1]-length(trtDelay)))
    vecBA <- c(trtDelay,rep(1,lenTp[2]-length(trtDelay)))

    trt[1,seq_len(lenTp[1])]              <- vecAB
    trt[2,(lenTp[1]+1):timepoints] <- vecBA
  }
  ## force matrix type (needed if timepoints==1)
  trtMat <- as.matrix(trt[rep(seq_len(lenCl),Cl),])

  return(trtMat)
}


#' @title Construct the time period adjustment in the design matrix
#'
#' @description Offers several options to adjust for secular trends.
#'
#' @inheritParams construct_DesMat
#'
#' @return a matrix with one row for every cluster at every timepoint and number of columns
#' depending of adjustment type.
#'
#' @export

construct_timeAdjust <- function(Cl,
                                 timepoints,
                                 timeAdjust = "factor",
                                 period     = NULL,
                                 timeBlk    = NULL){

  SumCl   <- sum(Cl)
  if(!is.null(timeBlk)) {
    timepoints <- dim(timeBlk)[1]
    timeBlks   <- timeBlk[rep(seq_len(timepoints),SumCl),]
    return(timeBlks)
  }

  if(timepoints==1) timeAdjust <- "none"
  if(timeAdjust=="periodic" & is.null(period)) period <- timepoints

  timeBlks <- switch (timeAdjust,
    factor   = cbind(1,rbind(0,diag(timepoints-1))
                     )[rep(seq_len(timepoints),SumCl),]
    ,
    none     = matrix(rep(1,timepoints*SumCl))
    ,
    linear   = cbind(rep(1,timepoints*SumCl),
                     rep(seq_len(timepoints)/timepoints,SumCl))
    ,
    periodic = cbind(rep(1,timepoints),
                     sin(0:(timepoints-1)*(2*pi/period)),
                     cos(0:(timepoints-1)*(2*pi/period))
                     )[rep(seq_len(timepoints),SumCl),]
    ,
    quadratic= cbind(rep(1,timepoints*SumCl),
                     rep(seq_len(timepoints)/timepoints,SumCl),
                     rep(seq_len(timepoints)/timepoints,SumCl)^2)
  )

  return(timeBlks)
}
