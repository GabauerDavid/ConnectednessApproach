
#' @title Minnesota Prior
#' @description Get Minnesota Prior
#' @param gamma Diagonal value of variance-covariance matrix
#' @param k Number of series
#' @param nlag Lag length
#' @return Get Minnesota Prior
#' @examples
#' \donttest{
#' prior = MinnesotaPrior(0.1, k=4, nlag=1)
#' }
#' @references Koop, G., & Korobilis, D. (2010). Bayesian multivariate time series methods for empirical macroeconomics. Now Publishers Inc.
#' @author David Gabauer
#' @export
MinnesotaPrior = function(gamma=0.1, k, nlag) {
  if (k<=1) {
    stop("k represents the number of series and needs to be an integer larger than 1")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  m = nlag*(k^2)
  bprior = c(cbind(0*diag(k), matrix(0, ncol=(nlag-1)*k, nrow=k)))
  V_i = matrix(0, nrow=(m/k), ncol=k)
  for (i in 1:k) {
    for (j in 1:(m/k)) {
      V_i[j,i] = gamma/(ceiling(j/k)^2)
    }
  }
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  return = list(bprior=bprior, Vprior=Vprior)
}
