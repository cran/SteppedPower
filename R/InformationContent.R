

#' Title  Formula-based calculation of information content
#'
#' @param CovMat #' @param CovMat numeric, a positive-semidefinite matrix with
#' (#Clusters \eqn{\cdot} timepoints) rows and columns.
#' @param W numeric, the inverse of a covariance matrix. If CovMat is specified,
#' input for W is ignored
#' @param dsn a matrix with  (#Clusters \eqn{\cdot} #timepoints) rows and p
#' columns, where p are the degrees of freedom of fixed effects in a gls model.
#' This usually contains the intervention effect and some specification of the
#' time effect.
#' @param sumCl number of clusters
#' @param tp number of time points
#'
#' @return A matrix containing the information content for every cluster-period cell
#'
#' @export

compute_InfoContent <- function(CovMat = NULL,
                                W      = NULL,
                                dsn,
                                sumCl,
                                tp){
  if(is.null(W) | !is.null(CovMat)) W  <- spdinv(as.matrix(CovMat))
  X2 <- dsn[,-1]
  x1 <- dsn[, 1]
  Q  <- W %*% X2 %*% spdinv(t(X2)%*%W%*%X2) %*%t(X2)%*%W
  WQ <- W - Q
  hh <- as.numeric( t(x1)%*%WQ / c(t(x1)%*%WQ%*%x1) )
  return(matrix(1/(1 - hh^2*as.numeric(t(x1)%*%WQ%*%x1) / diag(WQ) ),
                sumCl,tp, byrow=TRUE) )
}
