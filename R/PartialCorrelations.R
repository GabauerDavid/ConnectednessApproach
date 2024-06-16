#' @title Partial Contemporaneous Correlations
#'
#' @description Get partial contemporaneous correlations
#' @param Q variance-covariance matrix
#' @return Get partial contemporaneous correlations
#' @examples
#' \donttest{
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' pcc = PartialCorrelations(fit$Q)
#' }
#' @references Dahlhaus, R., & Eichler, M. (2003). Causality and graphical models in time series analysis. Oxford Statistical Science Series, 115-137.
#' @author David Gabauer
#' @export
PartialCorrelations = function (Q) {
  message("Partial correlations are computed according to:\n Dahlhaus, R., & Eichler, M. (2003). Causality and graphical models in time series analysis. Oxford Statistical Science Series, 115-137.")
  if (length(dim(Q))<=1) {
    stop("Q needs to be at least a 2-dimensional matrix")
  }
  k = dim(Q)[1]
  NAMES = colnames(Q)
  if (length(dim(Q))==2) {
    Q = array(Q, c(k,k,1), dimnames=list(NAMES,NAMES))
  }
  pcc = Q
  for (l in 1:dim(Q)[3]) {
    precision = MASS::ginv(Q[,,l])
    theta = diag(1/sqrt(diag(precision)))
    pcc[,,l] = -theta %*% precision %*% theta
  }
  return(pcc)
}
