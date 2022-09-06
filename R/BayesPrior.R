
#' @title Bayes Prior
#' @description Get Bayes prior
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param size Sample size used to calculate prior parameters
#' @return Get Bayes Prior
#' @examples
#' data(dy2012)
#' prior = BayesPrior(dy2012, nlag=1)
#' @importFrom methods is
#' @references Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy. The Review of Economic Studies, 72(3), 821-852.
#' @author David Gabauer
#' @export
BayesPrior = function(x, size=NULL, nlag) {
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (is.null(size)) {
    size = nrow(x)
  }
  if ((size-nlag)<=0) {
    stop("size needs to be larger than nlag")
  }
  
  p = ncol(x)
  size = size - nlag
  yt = t(x[(nlag+1):(size+nlag),])
  m = p + nlag*(p^2)
  Zt = NULL
  for (i in (nlag+1):(size+nlag)) {
    ztemp = diag(p)
    for (j in 1:nlag) {
      xlag = x[(i-j),1:p]
      xtemp = matrix(0,p,p*p)
      for (jj in 1:p) {
        xtemp[jj,((jj-1)*p+1):(jj*p)] = xlag
      }
      ztemp = cbind(ztemp, xtemp)
    }
    Zt = rbind(Zt, ztemp)
  }
  vbar = matrix(0,m,m)
  xhy = matrix(0,m,1)
  for (i in 1:size) {
    zhat1 = Zt[((i-1)*p+1):(i*p),]
    vbar = vbar + t(zhat1)%*%zhat1
    xhy = xhy + t(zhat1)%*%as.matrix(yt[,i])
  }
  vbar = MASS::ginv(vbar)
  aols = vbar%*%xhy
  
  sse2 = matrix(0,p,p)
  for (i in 1:size) {
    zhat1 = Zt[((i-1)*p+1):(i*p),]
    sse2 = sse2 + (yt[,i] - zhat1%*%aols)%*%t(yt[,i] - zhat1%*%aols)
  }
  hbar = sse2/size
  return(list(bprior=aols[-c(1:p)],Vprior=vbar[-c(1:p),-c(1:p)],Q=hbar))
}
