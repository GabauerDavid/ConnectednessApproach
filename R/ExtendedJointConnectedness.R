
#' @title Balcilar et al. (2021) extended joint connectedness approach
#' @description This function provides extended joint connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Balcilar et al. (2021)
#' data("bgu2021")
#' prior = MinnesotaPrior(0.1, k=ncol(bgu2021), nlag=1)
#' fit = TVPVAR(bgu2021, configuration=list(l=c(0.99,0.99), nlag=1, prior=prior))
#' dca = ExtendedJointConnectedness(Phi=fit$B_t, Sigma=fit$Q_t, nfore=20)
#' dca$TABLE
#' }
#' @references
#' Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy, 73, 102219.
#' @author David Gabauer
#' @export
ExtendedJointConnectedness = function(Phi, Sigma, nfore) {
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

  date = dimnames(Sigma)[[3]]
  TCI = array(NA, c(t,1), dimnames=list(as.character(date), "TCI"))
  NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(as.character(date), NAMES))
  CT = PCI = NPDC = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, as.character(date)))

  pb = progress_bar$new(total=t)
  for (ij in 1:t) {
    # calculate the gFEVD
    gSOT = 100*FEVD(Phi[,,ij], Sigma[,,ij], nfore=nfore, type="time",
                    generalized=TRUE)$FEVD
    gSOI = mean(rowSums(gSOT-diag(diag(gSOT))))

    # calculate Xi (the forecast error covariance matrix)
    A = Wold(Phi[,,ij], nfore)  # the VMA coefficient matrices
    Xi_h = array(0,dim=c(k,k,nfore))
    for (h in 1:nfore) {
      Xi_h[,,h] = A[,,h]%*%Sigma[,,ij]%*%t(A[,,h]) # calculated Xi at each h
    }
    Xi = rowSums(Xi_h, dims=2) # sum them along THIRD dimension to form Xi  (note: because this is a row sum, dims=2, actually sums along the third dimension)
    I_K = diag(1,nrow=k,ncol=k)

    # Calculate the elimination matrices.
    M = array(0,dim=c(k,k-1,k))
    for (i in 1:k){
      M[,,i] = I_K[,-i] # calculate the elimination matrices
    }
    S_jnt_numerator_h = array(0,dim=c(k,nfore))
    for (i in 1:k) {
      for (h in 1:nfore){
        S_jnt_numerator_h[i,h] = I_K[i,]%*%A[,,h]%*%Sigma[,,ij]%*%M[,,i]%*%(ginv(t(M[,,i])%*%Sigma[,,ij]%*%M[,,i]))%*%t(M[,,i])%*%Sigma[,,ij]%*%t(A[,,h])%*%I_K[,i] #calculate the numerator of S_jnt at each h
      }
    }

    S_jnt_numerator = array(0,dim=c(k))
    for (i in 1:k) {
      S_jnt_numerator[i] = sum(S_jnt_numerator_h[i,]) # calculate the numerator of j_jnt  (sum over h)
    }

    S_jnt=array(0,dim=c(k))
    for (i in 1:k) {
      S_jnt[i] = (100*S_jnt_numerator[i])/Xi[i,i]
    }

    # calculate the joint spillover index (jSOI)
    gSOT_diag = gSOT
    diag(gSOT_diag) = 0
    jSOI = mean(S_jnt)
    lambda = S_jnt / apply(gSOT_diag, 1, sum)
    jSOT = gSOT
    colnames(jSOT)=rownames(jSOT)=NAMES
    for (i in 1:k) {
      jSOT[i,] = gSOT[i,]*lambda[i]
    }
    jSOT_diag = jSOT
    diag(jSOT_diag) = 0
    from_jnt = rowSums(jSOT_diag)
    to_jnt = colSums(jSOT_diag)
    jSOI = mean(to_jnt)
    diag(jSOT_diag) = 100 - from_jnt

    dca = ConnectednessTable(jSOT_diag/100)
    CT[,,ij] = dca$FEVD
    TO[ij,] = dca$TO
    FROM[ij,] = dca$FROM
    NET[ij,] = dca$NET
    NPDC[,,ij] = dca$NPDC
    TCI[ij,] = dca$TCI
    PCI[,,ij] = dca$PCI
    NPT[ij,] = dca$NPT
    INFLUENCE[,,ij] = dca$INFLUENCE
    pb$tick()
  }

  TABLE = ConnectednessTable(CT/100)$TABLE
  config = list(nfore=nfore, approach="Extended Joint", generalized=TRUE, corrected=FALSE)
  return = list(TABLE=TABLE, CT=CT/100, TCI=TCI, TO=TO, FROM=FROM,
                NET=NET, NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
}
