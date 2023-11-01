
#' @title Quantile vector autoregression
#' @description Estimation of a QVAR using equation-by-equation quantile regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param tau quantile between 0 and 1
#' @param method See methods for rq in quantreg package. Default is "fn".
#' @param configuration model configuration
#' @return Estimate QVAR model
#' @examples
#' #data(dy2012)
#' #fit = QVAR(dy2012, configuration=list(nlag=1, tau=0.5))
#' @importFrom quantreg rq
#' @references
#' White, H., Kim, T. H., & Manganelli, S. (2015). VAR for VaR: Measuring tail dependence using multivariate regression quantiles. Journal of Econometrics, 187(1), 169-188.
#' 
#' Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters, 204, 109891.
#' @author David Gabauer
#' @export
QVAR = function(x, configuration=list(nlag=1, tau=0.5, method="fn")) {
  tau = as.numeric(configuration$tau)
  nlag = as.numeric(configuration$nlag)
  if (is.null(configuration$method)) {
    configuration$method = "fn"
  }
  method = as.character(configuration$method)

  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if ((sum(tau <= 0 | tau >= 1))>0) {
    stop("tau needs to be within 0 and 1")
  }
  k = ncol(x)
  if (length(tau)!=k) {
    tau = rep(tau,k)
  }
  NAMES = colnames(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }

  Res = B = se = NULL
  for (i in 1:k) {
    z = embed(x, nlag+1)
    fit = rq(z[,i] ~ z[,-c(1:k)], tau=tau[i], method=method)
    B = rbind(B, fit$coefficients[-1])
    Res = cbind(Res, fit$residuals)
    se = rbind(se, summary(fit)$coefficients[-1,2])
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(zoo::index(rownames(x)),1)))
  results = list(B=B, Q=Q, se=se)
}
