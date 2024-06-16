
#' @title R2 decomposed connectedness from correlations
#' @description This function computes the R2 decomposed connectedness measures from correlations
#' @param R zoo correlation data matrix
#' @return Get R2 connectedness measures from correlation matrix
#' @import progress
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @references
#' Naeem, M. A., Chatziantoniou, I., Gabauer, D., & Karim, S. (2023). Measuring the G20 Stock Market Return Transmission Mechanism: Evidence From the R2 Connectedness Approach. International Review of Financial Analysis.
#' 
#' Balli, F., Balli, H. O., Dang, T. H. N., & Gabauer, D. (2023). Contemporaneous and lagged R2 decomposed connectedness approach: New evidence from the energy futures market. Finance Research Letters, 57, 104168.
#' @author David Gabauer
#' @export
R2Correlations = function(R) {
  date = as.character(dimnames(R)[[3]])
  
  k = ncol(R)
  if (length(dim(R))==2) {
    R = array(R, c(k,k,1))
  }
  t = dim(R)[3]
  NAMES = colnames(R)
  
  R2dec = array(1, c(k,k,t), dimnames=list(NAMES,NAMES,date))
  pb = txtProgressBar(max=t,style=3)
  for (j in 1:t) {
    setTxtProgressBar(pb, j)
    for (i in 1:k) {
      Rc = R[,,j]
      ryx = Rc[-i,i,drop=F]
      rxx = Rc[-i,-i]
      
      eigcovx = eigen(rxx, TRUE)
      rootcovx = eigcovx$vectors%*%diag(sqrt(abs(eigcovx$values)))%*%t(eigcovx$vectors)
      cd = rootcovx^2 %*% (solve(rootcovx)%*%ryx)^2
      if (sum(cd)>1) {
        R2dec[i,-i,j] = R2dec[i,-i,c(j-1)]
      } else {
        R2dec[i,-i,j] = cd
      }
    }
  }
  TABLE = ConnectednessTable(apply(R2dec,1:2,mean))$TABLE
  
  t0 = dim(R2dec)[3]
  TCI = array(0, c(t0, 1), dimnames=list(date, "TCI"))
  TO = FROM = NET = array(0, c(t0, k), dimnames=list(date, NAMES))
  FEVD = NPDC = array(0, c(k, k, t0), dimnames=list(NAMES, NAMES, date))
  for (i in 1:dim(R2dec)[3]) {
    ct = ConnectednessTable(R2dec[,,i])
    FEVD[,,i] = ct$FEVD
    NPDC[,,i] = ct$NPDC
    TO[i,] = ct$TO
    FROM[i,] = ct$FROM
    NET[i,] = ct$NET
    TCI[i,] = mean(ct$TO)
  }
  
  PCI = R2dec
  for (l in 1:dim(R2dec)[3]) {
    for (i in 1:k) {
      for (j in 1:k) {
        PCI[i,j,l] = 2*(R2dec[i,j,l]+R2dec[j,i,l])/(R2dec[i,j,l]+R2dec[j,i,l]+R2dec[i,i,l]+R2dec[j,j,l])
      }
    }
  }
  
  config = list(approach="R2")
  return = list(CT=R2dec, TO=TO, FROM=FROM, NET=NET, TCI=TCI, PCI=PCI, NPDC=NPDC, TABLE=TABLE, config=config)
}