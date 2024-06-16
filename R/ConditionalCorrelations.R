#' @title Partial Conditional Correlations
#' @description Compute partial conditional correlations
#' @param Q Variance-covariance matrix of dimension
#' @return Get partial conditional correlations
#' @examples
#' \donttest{
#' data("dy2012")
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' pcc = ConditionalCorrelation(fit$Q)
#' }
#' @author David Gabauer
#' @export
ConditionalCorrelation = function (Q) {
  if (length(dim(Q))<=1) {
    stop("Q needs to be at least a 2-dimensional matrix")
  }
  k = dim(Q)[1]
  NAMES = colnames(Q)
  if (length(dim(Q))==2) {
    Q = array(Q, c(k,k,1), dimnames=list(NAMES,NAMES))
  }
  R = Q
  for (i in 1:k) {
    for (j in 1:k) {
      R[i,j,] = Q[i,j,] / (sqrt(Q[i,i,])*sqrt(Q[j,j,]))
    }
  }
  return(R)
}
