
#' @title Vector autoregression
#' @description Estimation of a VAR using equation-by-equation OLS regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param method OLS estimates are based on Pearson, Spearman or Kendall correlation.
#' @param configuration model configuration
#' @return Estimate VAR model
#' @examples
#' data("dy2012")
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' @references Sims, C. A. (1980). Macroeconomics and reality. Econometrica, 1-48.
#' @author David Gabauer
#' @importFrom stats lm
#' @importFrom stats embed
#' @export
VAR = function (x, configuration = list(nlag = 1, method="pearson")) {
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  nlag = as.numeric(configuration$nlag)
  if (nlag <= 0) {
    stop("nlag needs to be a positive integer")
  }
  
  cov_method = function(Z, method) {
    D = diag(apply(Z,2,sd))
    R = cor(Z, method=method)
    if (configuration$method=="kendall") {
      R = sin(pi/2*R)
    }
    D %*% R %*% D
  }
  
  Res = B = se = NULL
  for (i in 1:k) {
    z = stats::embed(x, nlag + 1)
    Z = cbind(z[, -c(1:k)])
    XX = cov_method(Z, method=configuration$method)*nrow(Z)
    yX = matrix(cov_method(cbind(z[,i], Z), method=configuration$method)[1,-1], ncol=1)*nrow(Z)
    XXinv = solve(XX)
    beta = XXinv %*% yX
    resid = (z[,i] - mean(z[,i])) - c(Z%*%beta)
    resid = scale(resid, T, F)
    sigma2 = sum(resid^2) / (nrow(Z)-ncol(Z)-1)
    
    B = rbind(B, c(beta))
    se = rbind(se, sqrt(diag(sigma2*XXinv)))
    Res = cbind(Res, resid)
  }
  Q = array(t(Res) %*% Res/nrow(Res), c(k, k, 1), dimnames = list(NAMES, 
                                                                  NAMES, tail(as.character(zoo::index(x)), 1)))
  results = list(B = B, Q = Q, se = se)
}

