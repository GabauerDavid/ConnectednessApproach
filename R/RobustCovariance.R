
#' @title Robust Covariance
#' @description Estimation of robust covariance
#' @param x zoo data matrix
#' @param method Either "mcd", "mve", "ogk", "tyler", "sccm" (covariance), "pearson", "spearman" and "kendall" (correlation) or a combination of the two.
#' @return Estimate robust covariance matrix
#' @examples
#' data("dy2012")
#' RobustCovariance(dy2012, method=c("pearson"))
#' RobustCovariance(dy2012, method=c("mcd", "pearson"))
#' @author David Gabauer
#' @importFrom rrcov CovMcd
#' @importFrom rrcov CovMve
#' @importFrom rrcov CovSest
#' @importFrom stats cor
#' @export
RobustCovariance = function(x, method="pearson") {
  
  x = as.matrix(x)
  if ("mcd" %in% method) {
    V = rrcov::CovMcd(x)@cov
  } else if ("mve" %in% method) {
    V = rrcov::CovMve(x)@cov
  } else if ("ogk" %in% method) {
    V = rrcov::CovOgk(x)@cov
  } else if ("sscm" %in% method) {
    V = rrcov::CovSest(x)@cov
  } else if ("sde" %in% method) {
    V = rrcov::CovSde(x)@cov
  } else {
    V = cov(x)
  }
  
  if (length(method)>1 || method %in% c("spearman", "kendall")) {
    if ("spearman" %in% method) {
      R = stats::cor(x, method="spearman")
      V = S %*% R %*% S
    } else if ("kendall" %in% method) {
      R = sin(pi/2*stats::cor(x, method="kendall"))
      V = S %*% R %*% S
    }
  }
  V
}
