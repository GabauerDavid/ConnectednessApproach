#' @title Impulse response functions
#' @description This function calculates orthorgonalized/generalized impulse response functions of time or frequency domain.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual Variance-Covariance Matrix
#' @param nfore H-step ahead forecast horizon
#' @param orth Boolean
#' @return Orthorgonal/generalized time/frequency impulse response functions
#' @examples
#' \donttest{
#' data("dy2012")
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' irf = IRF(Phi=fit$B, Sigma=fit$Q, nfore=10, orth=TRUE)
#' }
#' @references
#' Stiassny, A. (1996). A spectral decomposition for structural VAR models. Empirical Economics, 21(4), 535-555.
#' 
#' Koop, G., Pesaran, M. H., & Potter, S. M. (1996). Impulse response analysis in nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. Economics Letters, 58(1), 17-29.
#' @author David Gabauer
#' @export
IRF = function (Phi, Sigma, nfore=10, orth=TRUE) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  
  if (length(dim(Sigma))>2) {
    Sigma = Sigma[,,1]
  }
  p = 0
  k = 0
  if (length(Phi) > 0) {
    k = dim(Phi)[1]
    k1 = dim(Phi)[2]
    p = floor(k1/k)
  }
  if (is.null(Sigma)) {
    Sigma = diag(rep(1, k))
  }
  if (orth) {
    m1 = eigen(Sigma)
    v1 = sqrt(m1$values)
    vv = diag(v1)
    Pmtx = m1$vectors
    Sh = Pmtx %*% vv %*% t(Pmtx)
  }
  if (k < 1) {
    k = 1
  }
  PSI = diag(rep(1, k))
  if (orth) {
    WGT = c(PSI %*% Sh)
  } else {
    WGT = c(PSI)
  }
  for (il in 1:nfore) {
    ilk = il * k
    tmp = matrix(0, k, k)
    if (p > 0) {
      iend = min(il, p)
      for (j in 1:iend) {
        jdx = (il - j)
        kdx = j * k
        tmp = tmp + Phi[, (kdx - k + 1):kdx] %*% PSI[,(jdx * k + 1):(jdx * k + k)]
      }
    }
    PSI = cbind(PSI, tmp)
    if (orth) {
      WGT = cbind(WGT, c(tmp %*% Sh))
    } else {
      WGT = cbind(WGT, c(tmp))
    }
  }
  wk1 = WGT
  for (i in 1:k^2) {
    wk1[i, ] = cumsum(WGT[i, ])
  }
  tdx = c(1:(nfore + 1)) - 1
  if (orth) {
    gmax = max(WGT)
    gmin = min(WGT)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
  } else {
    gmax = max(WGT)
    gmin = min(WGT)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
  }
  k = sqrt(nrow(WGT))
  irf = list()
  for (i in 1:k) {
    irf[[i]] = t(WGT[((i-1)*k+1):(i*k),])
  }
  return = list(irf=irf)
}
