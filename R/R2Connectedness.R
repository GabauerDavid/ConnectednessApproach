
#' @title R2 connectedness approach
#' @description This function computes the R2 connectedness measures
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param window.size Rolling-window size or Bayes Prior sample size
#' @param method "pearson", "spearman", or "kendall". "pearson" is default.
#' @param relative Boolean whether relative or absolute R2 should be used
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @return Get R2 connectedness measures
#' @examples
#' \donttest{
#' data("dy2012")
#' dca = R2Connectedness(dy2012, window.size=NULL, nlag=0, method="pearson")
#' dca$TABLE
#' }
#' @import progress
#' @references
#' Naeem, M. A., Chatziantoniou, I., Gabauer, D., & Karim, S. (2023). Measuring the G20 Stock Market Return Transmission Mechanism: Evidence From the R2 Connectedness Approach. International Review of Financial Analysis.
#' 
#' Balli, F., Balli, H. O., Dang, T. H. N., & Gabauer, D. (2023). Contemporaneous and lagged R2 decomposed connectedness approach: New evidence from the energy futures market. Finance Research Letters, 57, 104168.
#' @author David Gabauer
#' @export
R2Connectedness = function(x, window.size=NULL, nlag=0, method="pearson", relative=FALSE, corrected=FALSE) {
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  
  DATE = as.character(index(x))
  x = as.matrix(x)
  k = ncol(x)
  NAMES = colnames(x)
  Z = embed(x, nlag+1)
  
  if (is.null(window.size)) {
    window.size = nrow(Z)
    t0 = 1
  } else {
    window.size = window.size - nlag
    t0 = nrow(Z) - window.size + 1
  }
  date = tail(DATE, t0)
  
  CT = array(0, c(k, k, t0, nlag+1), dimnames=list(NAMES, NAMES, date, 0:nlag))
  pb = txtProgressBar(max=t0,style=3)
  for (j in 1:t0) {
    setTxtProgressBar(pb, j)
    for (i in 1:k) {
      R = cor(Z[j:(j+window.size-1),], method=method)
      if (method=="kendall") {
        R = sin(0.5*pi*R)
      }
      ryx = R[-i,i,drop=F]
      rxx = R[-i,-i]
      
      eigcovx = eigen(rxx, TRUE)
      rootcovx = eigcovx$vectors%*%diag(sqrt(eigcovx$values))%*%t(eigcovx$vectors)
      cd = rootcovx^2 %*% (solve(rootcovx)%*%ryx)^2
      
      CT[i,-i,j,1] = cd[c(1:(k-1))]
      if (nlag>0) {
        CT[i,  ,j,2] = apply(array(cd[-c(1:(k-1))], c(1,k,nlag)), 1:2, sum)
      }
    }
  }
  
  kl = 1
  dimensions = "TCI"
  if (nlag>0) {
    kl = 3
    dimensions = c("Overall", "Contemporaneous", "Lagged")
  }
  TCI = array(0, c(t0, kl), dimnames=list(date, dimensions))
  TO = FROM = NET = array(0, c(t0, k, kl), dimnames=list(date, NAMES, dimensions))
  NPDC = array(0, c(k, k, t0, kl), dimnames=list(NAMES, NAMES, date, dimensions))
  pb = txtProgressBar(max=t0, style=3)
  for (i in 1:t0) {
    setTxtProgressBar(pb, i)
    
    if (nlag>0) {
      ct = ConnectednessTable(CT[,,i,1])
      lt = ConnectednessTable(CT[,,i,2])
      at = ConnectednessTable(CT[,,i,2] + CT[,,i,1])
      
      TO[i,,1] = at$TO
      FROM[i,,1] = at$FROM
      NET[i,,1] = at$NET
      NPDC[,,i,1] = at$NPDC
      TCI[i,1] = at$TCI*(1+(corrected*(k/(k-1)-1)))
      
      TO[i,,2] = ct$TO
      FROM[i,,2] = ct$FROM
      NET[i,,2] = ct$NET
      NPDC[,,i,2] = ct$NPDC
      TCI[i,2] = ct$TCI*(1+(corrected*(k/(k-1)-1)))
      
      TO[i,,3] = lt$TO
      FROM[i,,3] = lt$FROM
      NET[i,,3] = lt$NET
      NPDC[,,i,3] = lt$NPDC
      TCI[i,3] = lt$TCI*(1+(corrected*(k/(k-1)-1)))
      
    } else {
      
      ct = ConnectednessTable(CT[,,i,1])
      TO[i,,1] = ct$TO
      FROM[i,,1] = ct$FROM
      NET[i,,1] = ct$NET
      NPDC[,,i,1] = ct$NPDC
      TCI[i,1] = ct$TCI*(1+(corrected*(k/(k-1)-1)))
    }
  }
  
  TABLE = ConnectednessTable(CT[,,,1])$TABLE
  if (nlag>0) {
    lt = ConnectednessTable(CT[,,,2])$TABLE
    at = ConnectednessTable(CT[,,,1] + CT[,,,2])$TABLE
    TABLE = list("Overall"=at, "Contemporaneous"=TABLE, "Lagged"=lt)
  }
  
  if (nlag==0) {
    TO = TO[,,1]
    FROM = FROM[,,1]
    NET = NET[,,1]
    NPDC = NPDC[,,,1]
  }
  
  config = list(nlag=nlag, approach="R2", window.size=window.size, method=method, relative=relative, corrected=corrected)
  return = list(CT=CT, TO=TO, FROM=FROM, NET=NET, TCI=TCI, NPDC=NPDC, TABLE=TABLE, config=config)
}
