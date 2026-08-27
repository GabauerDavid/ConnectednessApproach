#' @title Connectedness table
#' @description This function provides a standard connectedness table.
#' @param FEVD Forecast error variance decomposition
#' @param digit Number of decimal places
#' @return Get connectedness table
#' @examples
#' \donttest{
#' data("dy2012")
#' fit = VAR(dy2012, configuration=list(nlag=1))
#' fevd = FEVD(Phi=fit$B, Sigma=fit$Q, nfore=10, type="time", generalized=TRUE)$FEVD
#' dca = ConnectednessTable(fevd)
#' }
#' @export
ConnectednessTable = function(FEVD, digit = 3) {
  if (length(dim(FEVD))<=1) {
    stop("FEVD needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(FEVD)
  namesX = rownames(FEVD)
  kX = length(namesX)
  namesZ = NAMES[-c(1:kX)]
  kZ = length(namesZ)
  
  CT = apply(FEVD, 1:2, mean) * 100
  ctX = CT[1:kX,1:kX]
  ctZ = CT[1:kX,-c(1:kX)]
  
  OWN = diag(diag(ctX))
  TO = colSums(ctX - OWN)
  FROM = rowSums(ctX) - rowSums(OWN)
  FROMX = rowSums(CT) - rowSums(OWN)
  NET = TO - FROM
  TCI = mean(FROM)
  TCIX = mean(FROMX)
  cTCI = TCI*k/(k-1)
  
  NPDC = ctX - t(ctX)
  NPT = rowSums(NPDC<0)
  INFLUENCE = 100*abs(NPDC/t(t(ctX)+ctX))
  
  table = format(round(cbind(CT, FROM, FROMX), digit), nsmall = digit)
  to = c(format(round(c(TO, sum(TO)), digit), nsmall = digit))
  to = c(to, rep(NA, kZ+1))
  inc = c(format(round(colSums(ctX), digit), nsmall = digit))
  inc = c(inc, rep(NA, kZ), 'TCI', 'R2')
  tci = paste0(format(round(TCI, digit), nsmall = digit))
  tcix = paste0(format(round(TCIX, digit), nsmall = digit))
  net = c(format(round(NET, digit), nsmall = digit))
  net = c(net, rep(NA, kZ), tci, tcix)
  npt = c(format(round(c(NPT),digit),nsmall=digit),  rep(NA, kZ+2))
  
  TABLE = rbind(table, to, inc, net, npt)
  colnames(TABLE) = c(NAMES, "FROM", "R2")
  rownames(TABLE) = c(namesX, "TO", "Inc.Own", "NET", "NPT")
  
  if (kZ==0) {
    TABLE = TABLE[,-ncol(TABLE)]
  }
  TABLE[which(is.na(TABLE), arr.ind=TRUE)] = ""
  
  PCI = matrix(NA, kX, kX)
  for (i in 1:kX) {
    for (j in 1:kX) {
      PCI[i,j] = 200*(CT[i,j]+CT[j,i])/(CT[i,i]+CT[i,j]+CT[j,i]+CT[j,j])
    }
  }
  
  return = list(FEVD = CT, TCI = TCI, TCIX = TCIX, cTCI=cTCI, PCI=PCI, NPDC=NPDC,
                TO = TO, FROM = FROM, FROMX = FROMX, NET = NET, TABLE = TABLE,
                NPDC=NPDC,  NPT=NPT, INFLUENCE=INFLUENCE)
}
