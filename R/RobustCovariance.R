
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
RobustCovariance = function (x, method = "pearson") {
  x = as.matrix(x)
  if ("mcd" %in% method) {
    Q = rrcov::CovMcd(x)@cov
  } else if ("mve" %in% method) {
    Q = rrcov::CovMve(x)@cov
  } else if ("ogk" %in% method) {
    Q = rrcov::CovOgk(x)@cov
  } else if ("sscm" %in% method) {
    Q = rrcov::CovSest(x)@cov
  } else if ("sde" %in% method) {
    Q = rrcov::CovSde(x)@cov
  } else {
    Q = cov(x)
  }
  
  R = ConditionalCorrelation(Q)[, , 1]
  S = diag(sqrt(diag(Q)))
  if (length(method) > 1 || method %in% c("spearman", "kendall")) {
    if ("spearman" %in% method) {
      R = 2 * sin(pi/6 * stats::cor(x, method = "spearman"))
      Q = S %*% R %*% S
    } else if ("kendall" %in% method) {
      R = sin(pi/2 * stats::cor(x, method = "kendall"))
      Q = S %*% R %*% S
    }
  }
  return(list(Q = Q, R = R))
}
