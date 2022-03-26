
#' @title Lastrapes and Wiesen (2021) joint connectedness approach
#' @description This function calculates the Lastrapes and Wiesen (2021) joint connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data(lw2021)
#' fit = VAR(lw2021, configuration=list(nlag=2))
#' dca = JointConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=30)
#' dca$TABLE
#' }
#' @references
#' Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.
#' @author David Gabauer
#' @export
JointConnectedness = function(Phi, Sigma, nfore) {
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  if (length(dim(Sigma))<=1) {
    stop("Sigma needs to be at least a 2-dimensional matrix")
  }
  if (length(dim(Phi))<=1) {
    stop("Phi needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(Sigma)
  if (length(dim(Phi))==2) {
    Phi = array(Phi, c(nrow(Phi),ncol(Phi),1))
  }
  if (length(dim(Sigma))==2) {
    Sigma = array(Sigma, c(nrow(Sigma),ncol(Sigma),1))
  }

  k = ncol(Sigma)
  t = dim(Sigma)[3]

  if (is.null(NAMES)) {
    NAMES = 1:k
  }

  date = as.character(dimnames(Sigma)[[3]])
  TCI = array(NA, c(t,1), dimnames=list(date, "TCI"))
  NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
  CT = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
  CT_ = array(NA, c(k+3, k+1, t))
  pb = progress_bar$new(total=t)
  for (ij in 1:t) {
    # calculate the gFEVD
    gSOT = 100*FEVD(Phi[,,ij], Sigma[,,ij], nfore=nfore, type="time", generalized=TRUE)$FEVD
    gSOI = mean(rowSums(gSOT-diag(diag(gSOT))))

    # calculate Xi (the forecast error covariance matrix)
    A = Wold(Phi[,,ij], nfore)  # the VMA coefficient matrices
    Xi_h = array(0,dim=c(k,k,nfore))
    for (h in 1:nfore) {
      Xi_h[,,h] = A[,,h]%*%Sigma[,,ij]%*%t(A[,,h]) # calculated Xi at each h
    }
    Xi = rowSums(Xi_h, dims=2) # sum them along THIRD dimension to form Xi  (note: because this is a row sum, dims=2, actually sums along the third dimension)
    I_K = diag(1,nrow=k,ncol=k)

    #Calculate the elimination matrices.
    #Usually denoted as a KxK-1 matrix M_i, here it is an array where M[,,1]=M_1, and in general M[,,i]=M_i.
    M = array(0,dim=c(k,k-1,k))
    for (i in 1:k){
      M[,,i] = I_K[,-i]  #calculate the elimination matrices
    }

    #Calculate the joint total spillovers from all others to variable i (S_jnt)
    #calculate the numerator of S_jnt
    S_jnt_numerator_h = array(0,dim=c(k,nfore))
    for (i in 1:k){
      for (h in 1:nfore){
        S_jnt_numerator_h[i,h] = I_K[i,]%*%A[,,h]%*%Sigma[,,ij]%*%M[,,i]%*%(solve(t(M[,,i])%*%Sigma[,,ij]%*%M[,,i]))%*%t(M[,,i])%*%Sigma[,,ij]%*%t(A[,,h])%*%I_K[,i]     #calculate the numerator of S_jnt at each h
      }
    }
    S_jnt_from = array(0,dim=c(k))
    for (i in 1:k){
      S_jnt_from[i] = 100*sum(S_jnt_numerator_h[i,])/Xi[i,i]
    }

    jSOI = mean(S_jnt_from)
    lambda = gSOI/jSOI
    TCI[ij,] = jSOI
    FROM[ij,] = S_jnt_from
    TO[ij,] = colSums(gSOT-diag(diag(gSOT)))/lambda
    NET[ij,] = TO[ij,] - FROM[ij,]
    CT[,,ij] = gSOT
    CT_[,,ij] = rbind(rbind(rbind(cbind(gSOT, FROM[ij,]), c(TO[ij,], sum(TO[ij,]))), c(colSums(gSOT), 0)), c(NET[ij,], TCI[ij,]))
    pb$tick()
  }
  TABLE = apply(CT_,1:2,mean)
  rownames(TABLE) = c(NAMES, "TO", "Inc.Own", "NET")
  colnames(TABLE) = c(NAMES, "FROM")
  TABLE = format(round(TABLE, 2), nsmall=2)
  TABLE[k+2,k+1] = "TCI"
  return = list(TABLE=TABLE, CT=CT/100, TCI=TCI, TO=TO, FROM=FROM, NET=NET,
                NPDC=NULL, NPT=NULL, PCI=NULL, INFLUENCE=NULL, approach="Joint")
}
