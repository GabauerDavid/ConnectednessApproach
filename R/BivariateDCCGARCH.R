
#' @title Bivariate DCC-GARCH
#' @description This function multiple Bivariate DCC-GARCH models that captures more accurately conditional covariances and correlations
#' @param data zoo dataset
#' @param spec A cGARCHspec A cGARCHspec object created by calling cgarchspec.
#' @param copula "mvnorm" or "mvt" (see, rmgarch package)
#' @param method "Kendall" or "ML" (see, rmgarch package)
#' @param transformation "parametric", "empirical" or "spd" (see, rmgarch package)
#' @param time.varying Boolean value to either choose DCC-GARCH or CCC-GARCH
#' @param asymmetric Whether to include an asymmetry term to the DCC model (thus estimating the aDCC).
#' @return Estimate Bivariate DCC-GARCH
#' @examples
#' #data("g2020")
#' #uspec = ugarchspec()
#' #fit = BivariateDCCGARCH(data=g2020, spec=replicate(nol(g2020), uspec))
#' @references
#' Antonakakis, N., Chatziantoniou, I. & Gabauer, D. (2022).
#' @author David Gabauer
#' @export
BivariateDCCGARCH = function(data, spec, copula=c("mvnorm","mvt"), method=c("Kendall","ML"), transformation=c("parametric", "empirical", "spd"), time.varying=TRUE, asymmetric=FALSE) {
  copula = match.arg(copula)
  method = match.arg(method)
  transformation = match.arg(transformation)
  t = nrow(data)
  k = ncol(data)
  Z_t = NULL
  H_t = R_t = array(NA, c(k,k,t))
  for (i in 1:k) {
    for (j in 1:k) {
      if (i>j) {
        mgarch.spec = cgarchspec(uspec=multispec(c(spec[i], spec[j])), dccOrder=c(1,1), asymmetric=FALSE,
                                 distribution.model=list(copula=copula, method=method, time.varying=time.varying, transformation=transformation))
        copula_fit = cgarchfit(mgarch.spec, data=data[,c(i,j)], solver=c("hybrid", "solnp"), fit.control=list(eval.se=FALSE))
        r = rcor(copula_fit)
        R_t[c(i,j),c(i,j),] = r
        h = rcov(copula_fit)
        H_t[c(i,j),c(i,j),] = h
      }
    }
    Z_t = cbind(Z_t, copula_fit@mfit$Z[,1])
  }
  return = list(H_t=H_t, R_t=R_t, Z_t=Z_t)
}

