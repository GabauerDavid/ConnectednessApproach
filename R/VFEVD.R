
#' @title Generalized volatility forecast error variance decomposition and volatility impulse response functions
#' @description This function provides the volatility impulse responses and the forecast error variance decomposition of DCC-GARCH models.
#' @param fit Fitted DCC-GARCH model
#' @param nfore H-step ahead forecast horizon
#' @param standardize Boolean value whether GIRF should be standardized
#' @return Get volatility impulse response functions and forecast error variance decomposition
#' @references
#' Gabauer, D. (2020). Volatility impulse response analysis for DCC‐GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.
#' @author David Gabauer
#' @importFrom rmgarch rcor
#' @importFrom rmgarch rcov
#' @export
VFEVD = function(fit, nfore=100, standardize=FALSE) {
  if (!is(fit, "DCCfit")) {
    stop("fit needs to be of class DCCfit")
  }
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  NAMES = fit@model$modeldata$asset.names
  H = rmgarch::rcov(fit)
  R = rmgarch::rcor(fit)
  R.bar = apply(R,1:2,mean)
  Q.bar = fit@mfit$Qbar
  t = dim(H)[3]
  k = dim(H)[1]
  alpha = array(0,c(k,k,nfore))
  alpha[,,1] = diag(fit@mfit$matcoef[c(seq(3,(4*k),4)),1])
  beta = diag(fit@mfit$matcoef[c(seq(4,(4*k),4)),1])
  ALPHA = fit@mfit$matcoef[(4*k+1),1]
  BETA = fit@mfit$matcoef[(4*k+2),1]

  H.hat = array(0,c(k,k,nfore+1))
  VIRF = H.hat.shock = H.hat.no_shock = array(0,c(k,k,t,nfore+1))
  e = diag(k)
  for (i in 1:t) {
    H.hat[,,1] = H[,,i]
    Q.hat = H.hat
    Q.hat[,,1] = fit@mfit$Q[[i]]
    for (j in 1:nfore) {
      H.hat[,,j+1] = (alpha[,,j])%*%e^2 + beta%*%H.hat[,,j]
      D = diag(diag(H.hat[,,j+1])^0.5)
      u = D%*%e
      if (j==1) {
        Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + ALPHA*crossprod(u) + BETA*H.hat[,,1]
      } else {
        Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar + (ALPHA+BETA)*Q.hat[,,j]
      }
      R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
      H.hat[,,j+1] = D%*%R.hat%*%D
    }
    H.hat.shock[,,i,] = H.hat
  }
  if (standardize) {
    e = 0*diag(k)
    for (i in 1:t) {
      H.hat[,,1] = H[,,i]
      Q.hat = H.hat
      Q.hat[,,1] = fit@mfit$Q[[i]]
      for (j in 1:nfore) {
        H.hat[,,j+1] = beta%*%H.hat[,,j]
        D = diag(diag(H.hat[,,j+1])^0.5)
        if (j==1) {
          Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + BETA*H.hat[,,1]
        } else {
          Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar+(ALPHA+BETA)*Q.hat[,,j]
        }
        R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
        H.hat[,,j+1] = D%*%R.hat%*%D
      }
      H.hat.no_shock[,,i,] = H.hat
    }
  }
  for (i in 1:t) {
    VIRF[,,i,] = H.hat.shock[,,i,] - H.hat.no_shock[,,i,]
  }
  date = dimnames(H)[[3]]
  VFEVD = array(NA, c(k,k,t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:t) {
    num = apply(VIRF[,,i,]^2,1:2,sum)
    den = c(apply(num,1,sum))
    fevd = t(num)/den
    VFEVD[,,i] = (fevd/apply(fevd, 1, sum))
  }
  return = list(IRF=VIRF, FEVD=VFEVD)
}
