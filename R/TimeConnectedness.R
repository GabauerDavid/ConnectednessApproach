
#' @title Diebold and Yilmaz (2009, 2012) connectedness approach
#' @description This function allows to calculate the Diebold and Yilmaz (2009, 2012) connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param generalized Orthorgonalized/generalized FEVD
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param FEVD Alternatively, to provide Phi and Sigma it is also possible to use FEVD directly.
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of DY2012
#' data("dy2012")
#' fit = VAR(dy2012, configuration=list(nlag=4))
#' dca = TimeConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=10, generalized=TRUE)
#' dca$TABLE
#' }
#' @references
#' Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal, 119(534), 158-171.\\
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting, 28(1), 57-66.
#' @author David Gabauer
#' @export
TimeConnectedness = function(Phi=NULL, Sigma=NULL, nfore=10, generalized=TRUE, corrected=FALSE, FEVD=NULL) {
  if ((is.null(Phi) || is.null(Phi)) && is.null(FEVD)) {
    stop("Either Sigma and Phi need to be given or FEVD")
  }
  if (is.null(FEVD)) {
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
    k = ncol(Sigma)
    if (length(dim(Phi))==2) {
      Phi = array(Phi, c(nrow(Phi),ncol(Phi),1))
    }
    if (length(dim(Sigma))==2) {
      Sigma = array(Sigma, c(nrow(Sigma),ncol(Sigma),1))
    }
    t = dim(Sigma)[3]
  
    if (is.null(NAMES)) {
      NAMES = 1:k
    }
    date = as.character(as.Date(dimnames(Sigma)[[3]]))
    CT = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, as.character(date)))
    for (i in 1:t) {
      CT[,,i] = FEVD(Phi=Phi[,,i], Sigma=Sigma[,,i], nfore=nfore, generalized=generalized, type="time")$FEVD
    }
  } else {
    CT = FEVD
    t = dim(CT)[3]
    NAMES = dimnames(CT)[[1]]
    k = length(NAMES)
    date = as.character(as.Date(dimnames(CT)[[3]]))
  }
  TCI = array(NA, c(t,1), dimnames=list(date, "TCI"))
  NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
  PCI = NPDC = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
  pb = progress::progress_bar$new(total=t)
  for (i in 1:t) {
    ct = ConnectednessTable(CT[,,i])
    TO[i,] = ct$TO
    FROM[i,] = ct$FROM
    NET[i,] = ct$NET
    NPT[i,] = ct$NPT
    NPDC[,,i] = ct$NPDC
    INFLUENCE[,,i] = ct$INFLUENCE
    PCI[,,i] = ct$PCI
    if (corrected) {
      TCI[i,] = ct$cTCI
    } else {
      TCI[i,] = ct$TCI
    }
    pb$tick()
  }
  TABLE = ConnectednessTable(CT)$TABLE
  config = list(nfore=nfore, approach="Time", generalized=generalized, corrected=corrected)
  return = list(TABLE=TABLE, CT=CT, TCI=TCI, TO=TO, FROM=FROM,
                NET=NET, NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
}
