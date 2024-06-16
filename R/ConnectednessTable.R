
#' @title Connectedness table
#' @description This function provides standard connectedness table.
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
ConnectednessTable = function(FEVD, digit=2) {
  if (length(dim(FEVD))<=1) {
    stop("FEVD needs to be at least a 2-dimensional matrix")
  }
  NAMES = colnames(FEVD)
  k = dim(FEVD)[1]
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  CT = apply(FEVD,1:2,mean)*100 # spillover from others to one specific
  OWN = diag(diag(CT))
  TO = colSums(CT-OWN)
  FROM = rowSums(CT-OWN)
  NET = TO-FROM
  TCI = mean(TO)
  cTCI = TCI*k/(k-1)
  NPDC = CT-t(CT)
  NPT = rowSums(NPDC<0)
  INFLUENCE = 100*abs(NPDC/t(t(CT)+CT))
  table = format(round(cbind(CT,FROM),digit),nsmall=digit)
  to = c(format(round(c(TO,sum(TO)),digit),nsmall=digit))
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "cTCI/TCI")
  tci = paste0(format(round(cTCI,digit),nsmall=digit),"/",format(round(TCI,digit),nsmall=digit))
  net = c(format(round(NET,digit),nsmall=digit))
  net = c(net, tci) 
  npt = c(format(round(NPT,digit),nsmall=digit), "")
  
  TABLE = rbind(table,to,inc,net,npt)
  colnames(TABLE) = c(NAMES,"FROM")
  rownames(TABLE) = c(NAMES,"TO","Inc.Own","NET","NPT")
  PCI = matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      PCI[i,j] = 200*(CT[i,j]+CT[j,i])/(CT[i,i]+CT[i,j]+CT[j,i]+CT[j,j])
    }
  }
  return = list(FEVD=CT, TCI=TCI, cTCI=cTCI, PCI=PCI,
                TO=TO, FROM=FROM, NET=NET, NPDC=NPDC, TABLE=TABLE,
                NPT=NPT, INFLUENCE=INFLUENCE)
}
