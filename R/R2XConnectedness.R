
#' @title R2 connectedness approach with exogenous input
#' @description This function computes the R2 connectedness measures with exogenous input
#' @param x zoo data matrix
#' @param z zoo data matrix for exogenous input
#' @param nlag Lag length
#' @param xlag Lag length of exogenous variables
#' @param window.size Rolling-window size or Bayes Prior sample size
#' @param method Methods for R2 connectedness are:"pearson", "spearman", or "kendall". "pearson" is default. Methods for pseudo R2 quantile connectedness are: "lasso", "br", "fn", or "sfn". "lasso" is default.
#' @param relative Boolean whether relative or absolute R2 should be used
#' @return Get R2 connectedness measures
#' @examples
#' \donttest{
#' data("dy2012")
#' dca = R2XConnectedness(x=dy2012[,1:2, drop=F], z=dy2012[,3:4, drop=F], window.size=NULL, nlag=1, xlag=1, method="pearson")
#' dca$TABLE
#' }
#' @import progress
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @references
#' Naeem, M. A., Chatziantoniou, I., Gabauer, D., & Karim, S. (2023). Measuring the G20 Stock Market Return Transmission Mechanism: Evidence From the R2 Connectedness Approach. International Review of Financial Analysis.
#' 
#' Balli, F., Balli, H. O., Dang, T. H. N., & Gabauer, D. (2023). Contemporaneous and lagged R2 decomposed connectedness approach: New evidence from the energy futures market. Finance Research Letters, 57, 104168.
#' 
#' Ferrer, R., Shahzad, S. J. H., Furió, D., & Benammar, R. (2025). Systemic Risk in the Tails: Contemporaneous Transmission and Spillover Dynamics in European Renewable Energy Equities.
#' @author David Gabauer
#' @export
R2XConnectedness = function(x, z=NULL, window.size=NULL, nlag=1, xlag=0, method="pearson", relative=FALSE) {
  
  k = ncol(x)
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  
  if (is(z, "list")) {
    if (length(z)==k) {
      exog_data = z
    } else if (length(z)==1) {
      exog_data = list()
      for (idx in 1:k) {
        exog_data[[idx]] = z[[1]]
      }
    }
  }
  
  if (length(dim(z))==2) {
    exog_data = list()
    for (idx in 1:k) {
      exog_data[[idx]] = z
    }
  }
  
  DATE = as.character(index(x))
  namesX = colnames(x)
  namesZ = NULL
  if (!is.null(exog_data)) {
    namesZ = colnames(exog_data[[1]])
  }
  NAMES = c(namesX, namesZ)
  
  # number of covariates
  kX = length(namesX)
  kZ = length(namesZ)
  X = embed(x, nlag + 1)
  
  if (is.null(window.size)) {
    window.size = nrow(X)
    t0 = 1
  } else {
    window.size = window.size - nlag
    t0 = nrow(X) - window.size + 1
  }
  date = tail(DATE, t0)
  
  betaX = array(0, c(kX*(nlag + 1), t0), dimnames = list(rep(namesX, nlag+1), date))
  betaZ = array(0, c(kX, (xlag + 1)*kZ, t0), dimnames = list(namesX, rep(namesZ, xlag+1), date))
  
  ctX = array(0, c(kX, kX+kZ, t0, nlag + 1), dimnames = list(namesX, NAMES, date, 0:nlag))
  pb = txtProgressBar(max = t0, style = 3)
  for (j in 1:t0) {
    setTxtProgressBar(pb, j)
    for (i in 1:kX) {
      
      Z = NULL
      if (!is.null(exog_data)) {
        Z = embed(exog_data[[i]], xlag+1)
        Z = tail(Z, nrow(X))
        Z = as.matrix(Z[j:(j + window.size - 1), ])
      }
      ZZ = cbind(X[j:(j + window.size - 1), ], Z)
      
      R = cor(ZZ, method=method)
      if (method == "kendall") {
        R = sin(0.5 * pi * R)
      } else if (method == "spearman") {
        R = 2 * sin(pi/6 * R)
      }
      
      beta = summary(lm(ZZ[,i]~ZZ[,-i]))$coefficients[,1]
      kz = ifelse(is.null(ncol(Z)), 0, ncol(Z))
      betaX[-i,j] = head(beta[-1], length(beta)-kz-1)
      if (kz!=0) {
        betaZ[i,,j] = tail(beta, ncol(Z))
      }
      
      ryx = R[-i, i, drop = F]
      rxx = R[-i, -i]
      eigcovx = eigen(rxx, TRUE)
      rootcovx = eigcovx$vectors %*% diag(sqrt(eigcovx$values)) %*% t(eigcovx$vectors)
      R2d = rootcovx^2 %*% (solve(rootcovx) %*% ryx)^2
      R2X = head(R2d, length(R2d)-kz)
      R2Z = array(tail(R2d, kz), c(xlag+1, kZ))
      
      ctX[i, -i, j, 1] = c(R2X[c(1:(kX - 1)),1], R2Z[1,]) # contemporaneous
      if (nlag > 0) {
        ctX[i, , j, 2] = c(apply(
          array(R2X[-c(1:(kX - 1)),1], c(1, kX, nlag)), 1:2, sum), 
          apply(R2Z[-1,,drop=F],2,sum)
        )
      }
    }
  }
  
  kl = 1
  dimensions = "TCI"
  if (nlag > 0) {
    kl = 3
    dimensions = c("Overall", "Contemporaneous", "Lagged")
  }
  TCIX = TCI = array(0, c(t0, kl), dimnames = list(date, dimensions))
  TO = FROMX = FROM = NET = array(0, c(t0, kX, kl), dimnames = list(date, namesX, dimensions))
  NPDC = array(0, c(kX, kX, t0, kl), dimnames = list(namesX, namesX, date, dimensions))
  
  pb = txtProgressBar(max = t0, style = 3)
  for (i in 1:t0) {
    setTxtProgressBar(pb, i)
    if (nlag > 0) {
      ct = ConnectednessTable(ctX[, , i, 1])
      lt = ConnectednessTable(ctX[, , i, 2])
      at = ConnectednessTable(ctX[, , i, 2] + ctX[,,i,1])
      
      TO[i, , 1] = at$TO
      TO[i, , 2] = ct$TO
      TO[i, , 3] = lt$TO
      
      FROM[i, , 1] = at$FROM
      FROM[i, , 2] = ct$FROM
      FROM[i, , 3] = lt$FROM
      
      FROMX[i, , 1] = at$FROMX
      FROMX[i, , 2] = ct$FROMX
      FROMX[i, , 3] = lt$FROMX
      
      NET[i, , 1] = at$NET
      NET[i, , 2] = ct$NET
      NET[i, , 3] = lt$NET
      
      TCI[i, 1] = at$TCI
      TCI[i, 2] = ct$TCI
      TCI[i, 3] = lt$TCI
      
      TCIX[i, 1] = at$TCIX
      TCIX[i, 2] = ct$TCIX
      TCIX[i, 3] = lt$TCIX
      
      NPDC[,,i,1] = at$NPDC
      NPDC[,,i,2] = ct$NPDC
      NPDC[,,i,3] = lt$NPDC
      
    } else {
      ct = ConnectednessTable(ctX[, , i, 1])
      TO[i, , 1] = ct$TO
      FROM[i, , 1] = ct$FROM
      NET[i, , 1] = ct$NET
      NPDC[, , i, 1] = ct$NPDC
      TCI[i, 1] = ct$TCI
    }
  }
  
  TABLE = ConnectednessTable(ctX[, , , 1])$TABLE
  if (nlag > 0) {
    lt = ConnectednessTable(ctX[, , , 2])$TABLE
    at = ConnectednessTable(ctX[, , , 1] + ctX[, , , 2])$TABLE
    TABLE = list(Overall = at, Contemporaneous = TABLE, Lagged = lt)
  }
  
  if (nlag == 0) {
    TO = TO[, , 1]
    FROM = FROM[, , 1]
    NET = NET[, , 1]
    NPDC = NPDC[, , , 1]
  }
  config = list(nlag = nlag, approach = "R2", window.size = window.size, 
                method = method, relative = relative)
  return = list(CT = ctX, TO = TO, FROM = FROM, R2=FROMX, NET = NET, NPDC=NPDC,
                TCI = TCI, TCIX = TCIX, TABLE = TABLE, config = config, betaX=betaX, betaZ=betaZ)
}
