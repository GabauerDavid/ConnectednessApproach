
#' @title Bivariate DCC-GARCH
#' @description This function multiple Bivariate DCC-GARCH models that captures more accurately conditional covariances and correlations
#' @param x zoo dataset
#' @param spec A cGARCHspec A cGARCHspec object created by calling cgarchspec.
#' @param copula "mvnorm" or "mvt" (see, rmgarch package)
#' @param method "Kendall" or "ML" (see, rmgarch package)
#' @param transformation "parametric", "empirical" or "spd" (see, rmgarch package)
#' @param time.varying Boolean value to either choose DCC-GARCH or CCC-GARCH
#' @param asymmetric Whether to include an asymmetry term to the DCC model (thus estimating the aDCC).
#' @return Estimate Bivariate DCC-GARCH
#' @importFrom rmgarch cgarchspec
#' @importFrom rmgarch cgarchfit
#' @importFrom rmgarch rcor
#' @importFrom rmgarch rcov
#' @references
#' Engle, R. (2002). Dynamic conditional correlation: A simple class of multivariate generalized autoregressive conditional heteroskedasticity models. Journal of Business & Economic Statistics, 20(3), 339-350.
#' @author David Gabauer
#' @export
BivariateDCCGARCH = function(x, spec, copula="mvt", method="Kendall", transformation="parametric", time.varying=TRUE, asymmetric=FALSE) {
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  t = nrow(x)
  k = ncol(x)
  NAMES = colnames(x)
  Z = NULL
  H = R = array(NA, c(k,k,t), dimnames=list(NAMES, NAMES, as.character(rownames(x))))
  for (i in 1:k) {
    for (j in 1:k) {
      if (i>j) {
        mgarch.spec = rmgarch::cgarchspec(uspec=multispec(c(spec[i], spec[j])), dccOrder=c(1,1), asymmetric=asymmetric,
                                          distribution.model=list(copula=copula, method=method, time.varying=time.varying, transformation=transformation))
        copula_fit = rmgarch::cgarchfit(mgarch.spec, data=x[,c(i,j)], solver=c("hybrid", "solnp"), fit.control=list(eval.se=FALSE))
        R[c(i,j),c(i,j),] = rcor(copula_fit)
        H[c(i,j),c(i,j),] = rcov(copula_fit)
      }
    }
    Z = cbind(Z, copula_fit@mfit$Z[,1])
  }
  return = list(H=H, R=R, Z=Z)
}

