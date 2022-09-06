
#' @title Vector autoregression
#' @description Estimation of a VAR using equation-by-equation OLS regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param configuration model configuration
#' @return Estimate VAR model
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' @references Sims, C. A. (1980). Macroeconomics and reality. Econometrica, 1-48.
#' @author David Gabauer
#' @importFrom stats lm
#' @importFrom stats embed
#' @export
VAR = function(x, configuration=list(nlag=1)) {
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  nlag = as.numeric(configuration$nlag)
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  
  Res = B = se = NULL
  for (i in 1:k) {
    z = stats::embed(x, nlag+1)
    fit = summary(stats::lm(z[,i] ~ z[,-c(1:k)]))
    B = rbind(B, fit$coefficients[-1,1])
    se = rbind(se, fit$coefficients[-1,2])
    Res = cbind(Res, fit$residuals)
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(as.character(zoo::index(x)),1)))
  results = list(B=B, Q=Q, se=se)
}
