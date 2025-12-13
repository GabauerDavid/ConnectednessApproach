
#' @title Quantile correlation
#' @description This function computes the quantile correlation coefficient proposed by Choi and Shin (2022).
#' @param x zoo data matrix
#' @param tau quantile between 0 and 1
#' @param method Either "lasso", "br", "fn" or "sfn". Default is "lasso"
#' @return Get quantile correlations
#' @importFrom zoo zoo
#' @importFrom zoo index
#' @examples
#' \donttest{
#' data("g2020")
#' fit = QuantileCorrelation(g2020, tau=0.5)
#' }
#' @references
#' Choi, J. E., & Shin, D. W. (2022). Quantile correlation coefficient: A new tail dependence measure. Statistical Papers, 63(4), 1075-1104.
#' @author David Gabauer
#' @export
QuantileCorrelation = function(x, tau=0.5, method="lasso") {
  stopifnot(tau > 0, tau < 1)

  # coerce to numeric matrix
  if (zoo::is.zoo(x) || xts::is.xts(x)) x <- zoo::coredata(x)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  n  <- ncol(x)
  nm <- colnames(x)
  if (is.null(nm)) nm <- paste0("V", seq_len(n))
  
  # result matrix
  corQ <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  diag(corQ) <- 1
  
  # choose faster quantile regression solver
  rq_fit <- switch(
    method,
    br  = function(y, X) quantreg::rq.fit.br (X, y, tau = tau)$coefficients,
    sfn = function(y, X) quantreg::rq.fit.sfn(X, y, tau = tau)$coefficients,
    fn  = function(y, X) quantreg::rq.fit.fnb(X, y, tau = tau)$coefficients,
    lasso  = function(y, X) quantreg::rq.fit.lasso(X, y, tau = tau)$coefficients,
    stop("Unknown 'method': use one of 'lasso', 'br', 'sfn', 'fn'", call. = FALSE)
  )
  
  # compute pairwise quantile correlations (upper triangle)
  inds <- which(upper.tri(corQ), arr.ind = TRUE)
  vals <- vapply(seq_len(nrow(inds)), function(k) {
    i <- inds[k, 1]; j <- inds[k, 2]
    xi <- x[, i]; xj <- x[, j]
    
    ok <- is.finite(xi) & is.finite(xj)
    xi <- xi[ok]; xj <- xj[ok]
    if (length(xi) < 3L) return(NA_real_)
    
    # y ~ 1 + x
    X1 <- cbind(1, xi)
    b2.1 <- rq_fit(y = xj, X = X1)[2]
    
    # x ~ 1 + y
    X2 <- cbind(1, xj)
    b1.2 <- rq_fit(y = xi, X = X2)[2]
    
    qcor <- sign(b2.1) * (if ((b2.1 * b1.2) > 0) sqrt(b2.1 * b1.2) else 0)
    if (qcor > 1) qcor <- 0.999999
    if (qcor < -1) qcor <- -0.999999
    qcor
  }, numeric(1))
  
  # fill symmetric matrix
  corQ[upper.tri(corQ)] = vals
  corQ[lower.tri(corQ)] = t(corQ)[lower.tri(corQ)]
  corQ
}
