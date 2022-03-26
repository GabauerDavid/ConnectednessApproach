
#' @title Uninformative Prior
#' @description Get Uninformative Prior
#' @param k Number of series
#' @param nlag Lag length
#' @return Get Uninformative Prior
#' @examples
#' prior = UninformativePrior(k=4, nlag=1)
#' @references Koop, G., & Korobilis, D. (2010). Bayesian multivariate time series methods for empirical macroeconomics. Now Publishers Inc.
#' @author David Gabauer
#' @export
UninformativePrior = function(k, nlag) {
  if (k<=1) {
    stop("k represents the number of series and needs to be an integer larger than 1")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  m = nlag*(k^2)
  bprior = c(cbind(0*diag(k), matrix(0, ncol=(nlag-1)*k, nrow=k)))
  V_i = matrix(4, nrow=(m/k), ncol=k)
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  return = list(bprior=bprior, Vprior=Vprior)
}
