
#' @title Forecast error variance decomposition
#' @description This function computes the orthogonalized/generalized forecast error variance decomposition
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param type Time or Frequency connectedness approach
#' @param range Partition range for frequency approach only.
#' @param generalized Generalized or orthogonalized FEVD
#' @param orth Orthogonalized shocks
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
FEVD = function (Phi, Sigma, nfore=100, type=c("time","frequency"), generalized=TRUE, orth=FALSE, range=NULL) {
  type = match.arg(type)
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
  irf = IRF(Phi=Phi, Sigma=Sigma, nfore=nfore, orth=orth)$irf
  if (type=="time") {
    Phi = lapply(1:nfore, function(j) sapply(irf, function(i) i[j,]))
    denom = diag(Reduce("+", lapply(Phi, function(i) i %*% Sigma %*% t(i))))
    if (generalized) {
      enum = Reduce("+", lapply(Phi, function(i) (i %*% Sigma)^2))
      tab = sapply(1:nrow(enum), function(j) enum[j,]/(denom[j] * diag(Sigma)))
      FEVD = t(apply(tab, 2, function(i) i/sum(i)))
    } else {
      enum = Reduce("+", lapply(Phi, function(i) (chol(Sigma) %*% t(i))^2))
      FEVD = t(sapply(1:ncol(enum), function(i) enum[,i]/denom[i]))
    }
  } else if (type=="frequency") {
    fftir = lapply(irf, function(i) apply(i, 2, fft))
    fftir = lapply(1:(nfore+1), function(j) sapply(fftir, function(i) i[j,]))
    denom = diag(Re(Reduce("+", lapply(fftir, function(i) i %*% Sigma %*% t(Conj(i))/nfore)[range])))
    if (generalized) {
      enum = lapply(fftir, function(i) (abs(i %*% Sigma))^2/(nfore + 1))
      tab = lapply(enum, function(i) sapply(1:nrow(i), function(j) i[j, ]/(denom[j] * diag(Sigma))))
      tot = apply(Reduce("+", tab[range]), 2, sum)
      FEVD = lapply(tab, function(i) t(i)/tot)
    } else {
      enum = lapply(fftir, function(i) (abs(i %*% t(chol(Sigma))))^2/(nfore + 1))
      FEVD = lapply(enum, function(i) t(sapply(1:nrow(i), function(j) i[j,]/(denom[j]))))
    }
  }
  return = list(IRF=irf, FEVD=FEVD)
}
