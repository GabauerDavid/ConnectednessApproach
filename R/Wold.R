#' @title Wold representation theorem
#' @description Transform VAR to VMA coefficients
#' @param x VAR coefficients
#' @param nfore H-step ahead forecast horizon
#' @return Get VMA coefficients
#' @examples
#' data(dy2012)
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' wold = Wold(fit$B, nfore=10)
#' @author David Gabauer
#' @export
Wold = function (x, nfore=10) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  nstep = abs(as.integer(nfore))
  K = nrow(x)
  p = floor(ncol(x)/K)
  A = array(0, c(K,K,nstep))
  for (i in 1:p){
    A[,,i]=x[,((i-1)*K+1):(i*K)]
  }
  Phi = array(0, dim = c(K, K, nstep + 1))
  Phi[,,1]=diag(K)
  Phi[,,2]=Phi[,,1] %*% A[,,1]
  if (nstep > 1) {
    for (i in 3:(nstep + 1)) {
      tmp1 = Phi[,,1] %*% A[,,i-1]
      tmp2 = matrix(0, nrow=K, ncol=K)
      idx = (i-2):1
      for (j in 1:(i-2)) {
        tmp2 = tmp2 + Phi[,,j+1] %*% A[,,idx[j]]
      }
      Phi[,,i] = tmp1 + tmp2
    }
  }
  return(Phi)
}
