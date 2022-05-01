
#' @title Forecast error variance decomposition
#' @description This function computes the orthogonalized/generalized forecast error variance decomposition
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param type Time or Frequency connectedness approach
#' @param range Partition range for frequency approach only.
#' @param generalized Generalized or orthogonalized FEVD
#' @return Orthogonalized/generalized time/frequency forecast error variance decomposition
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' fevd = FEVD(Phi=fit$B, Sigma=fit$Q, nfore=10, type="time", generalized=TRUE)$FEVD
#' @references
#' Stiassny, A. (1996). A spectral decomposition for structural VAR models. Empirical Economics, 21(4), 535-555.\\
#' Koop, G., Pesaran, M. H., & Potter, S. M. (1996). Impulse response analysis in nonlinear multivariate models. Journal of Econometrics, 74(1), 119-147.\\
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. Economics Letters, 58(1), 17-29.
#' @importFrom stats fft
#' @export
FEVD = function (Phi, Sigma, nfore=100, type=c("time", "frequency"), generalized=TRUE, range=NULL) {
  type = match.arg(type)
  if (nfore <= 0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma)) <= 1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi)) <= 1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Sigma)) > 2) {
    Sigma = Sigma[, , 1]
  }
  k = ncol(Sigma)
  irf = IRF(Phi=Phi, Sigma=Sigma, nfore=nfore, orth=FALSE)$irf
  if (type == "time") {
    Phi = lapply(1:nfore, function(j) sapply(irf, function(i) i[j, ])) # change representation
    denom = diag(Reduce("+", lapply(Phi, function(i) i %*% Sigma %*% t(i)))) # denominator of FEVD
    if (generalized) {
      irf_ = lapply(Phi, function(i) t(i %*% Sigma %*% solve(diag(sqrt(diag(Sigma))))))
      IRF = array(unlist(irf_), c(k,k,nfore))
      fevd = IRF^2
      enum = apply(fevd,1:2,sum)
      den = c(apply(enum,1,sum))
      nfevd = t(enum)/den
      FEVD = nfevd / apply(nfevd, 1, sum)
    } else {
      irf_ = lapply(Phi, function(i) chol(Sigma) %*%  t(i)) # orthorgonalized impulse response function
      IRF = array(unlist(irf_), c(k,k,nfore))
      fevd = IRF^2
      tab_ = array(0,c(k,k,nfore))
      for (i in 1:nfore) {
        tab_[,,i] = t(fevd[,,i] %*% solve(diag(denom)))
      }
      FEVD = apply(tab_,1:2,sum)
    }
  } else if (type == "frequency") {
    fftir = lapply(irf, function(i) apply(i, 2, fft))
    fftir = lapply(1:(nfore + 1), function(j) sapply(fftir, function(i) i[j, ]))
    denom = diag(Re(Reduce("+", lapply(fftir, function(i) i %*% Sigma %*% t(Conj(i))/nfore)[range])))
    if (generalized) {
      irf_ = lapply(fftir, function(i) t(abs(i %*% Sigma) %*% solve(diag(sqrt(diag(Sigma))))))
      IRF = array(unlist(irf_), c(k,k,nfore+1))
      irf_ = lapply(fftir, function(i) (abs(i %*% Sigma)))
      IRF_ = array(unlist(irf_), c(k,k,nfore+1))
      enum_ = IRF_^2/(nfore + 1)
      tab = list()
      for (i in 1:dim(enum_)[3]) {
        fevd_ = enum_[,,1]
        for (j in 1:dim(enum_)[2]) {
          fevd_[j,] = enum_[j,,i] / (denom[j] * diag(Sigma))
        }
        tab[[i]] = t(fevd_)
      }
      tot = apply(Reduce("+", tab[range]), 2, sum)
      FEVD = lapply(tab, function(i) t(i)/tot)
    } else {
      irf_ = lapply(fftir, function(i) (abs(i %*% t(chol(Sigma)))))
      IRF = array(unlist(irf_), c(k,k,nfore+1))
      enum = IRF^2/(nfore + 1)
      FEVD = list()
      for (i in 1:dim(enum)[3]) {
        fevd_ = enum[,,1]
        for (j in 1:dim(enum)[2]) {
          fevd_[j,] = enum[j,,i] / denom[j]
        }
        FEVD[[i]] = fevd_
      }
    }
  }
  return = list(IRF=IRF, FEVD=FEVD)
}
