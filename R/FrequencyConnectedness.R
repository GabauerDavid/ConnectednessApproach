#' @title Baruník and Křehlík (2018) frequency connectedness approach
#' @description This function calculates the Baruník and Křehlík (2018) frequency connectedness measures.
#' @param Phi VAR coefficient matrix
#' @param Sigma Residual variance-covariance matrix
#' @param nfore H-step ahead forecast horizon
#' @param partition Frequency spectrum
#' @param generalized Orthorgonalized/generalized FEVD
#' @param scenario ABS or WTH
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param orth Orthorgonalized shocks
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("dy2012")
#' partition = c(pi+0.00001, pi/4, 0)
#' fit = VAR(dy2012, configuration=list(nlag=4))
#' dca = FrequencyConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=100, partition=partition)
#' }
#' @import frequencyConnectedness
#' @references
#' Baruník, J., & Křehlík, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics, 16(2), 271-296.
#' @author David Gabauer
#' @export
FrequencyConnectedness = function(Phi, Sigma, nfore=100, partition=c(pi,pi/2,0), generalized=TRUE, orth=FALSE, scenario="ABS", corrected=FALSE) {
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

  k = dim(Sigma)[1]
  t = dim(Sigma)[3]
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  
  periods = round(pi/partition)
  period_names = NULL
  for (i in 1:(length(periods)-1)) {
    period_names = c(period_names, paste0(periods[i], "-", periods[i+1]))
  }
  period_names = c("Total",period_names)
  date = as.character(dimnames(Sigma)[[3]])
  interval = length(period_names)
  new_p = frequencyConnectedness::getPartition(partition, nfore)
  range = sort(unique(do.call(c, new_p)))
  
  TCI = array(0, c(t,interval), dimnames=list(date, period_names))
  PCI = INFLUENCE = CT = NPDC = array(0, c(k, k, t, interval), dimnames=list(NAMES, NAMES, date, period_names))
  NET = FROM = TO = array(0, c(t, k, interval), dimnames=list(date, NAMES, period_names))
  NPT = array(0, c(t, k, interval-1), dimnames=list(date, NAMES, period_names[-1]))
  PCI = INFLUENCE = array(0, c(k, k, t, interval-1), dimnames=list(NAMES, NAMES, date, period_names[-1]))
  pb = progress_bar$new(total=t)
  for (i in 1:t) {
    decomp = FEVD(Phi=Phi[,,i], Sigma=Sigma[,,i], nfore=nfore, generalized=generalized, type="frequency", range=range)$FEVD
    for (ij in 1:length(decomp)) {
      rownames(decomp[[ij]]) = colnames(decomp[[ij]]) = 1:ncol(Sigma)
    }
    tables = lapply(new_p, function(j) Reduce('+', decomp[j]))
    for (j in 2:(interval)) {
      if (scenario=="ABS") {
        dca = ConnectednessTable(tables[[j-1]])
        CT[,,i,j] = dca$FEVD
        TO[i,,j] = dca$TO
        FROM[i,,j] = dca$FROM
        NET[i,,j] = dca$NET
        NPDC[,,i,j] = dca$NPDC
        PCI[,,i,j-1] = dca$PCI
        INFLUENCE[,,i,j-1] = dca$INFLUENCE
        NPT[i,,j-1] = dca$NPT
        if (corrected) {
          TCI[i,j] = dca$cTCI
        } else {
          TCI[i,j] = dca$TCI
        }
      } else if (scenario=="WTH") {
        dca = ConnectednessTable(tables[[j-1]]/sum(sum(tables[[j-1]]))*k)
        CT[,,i,j] = dca$FEVD
        TO[i,,j] = dca$TO
        FROM[i,,j] = dca$FROM
        NET[i,,j] = dca$NET
        NPDC[,,i,j] = dca$NPDC
        PCI[,,i,j-1] = dca$PCI
        INFLUENCE[,,i,j-1] = dca$INFLUENCE
        NPT[i,,j-1] = dca$NPT
        if (corrected) {
          TCI[i,j] = dca$cTCI
        } else {
          TCI[i,j] = dca$TCI
        }
      }
    }
    pb$tick()
  }
  CT[,,,1] = apply(CT,1:3,sum)
  TCI[,1] = apply(TCI,1,sum)
  TO[,,1] = apply(TO,1:2,sum)
  FROM[,,1] = apply(FROM,1:2,sum)
  NET[,,1] = apply(NET,1:2,sum)
  NPDC[,,,1] = apply(NPDC,1:3,sum)

  TABLE = array(NA,c(k+4,k+1,interval), dimnames=list(c(NAMES, "TO", "Inc.Own", "Net", "NPDC"), c(NAMES, "FROM"), period_names))
  for (i in 1:interval) {
    TABLE[,,i] = ConnectednessTable(CT[,,,i]/100)$TABLE
  }
  config = list(partition=partition, nfore=nfore, generalized=generalized, orth=orth, scenario=scenario, corrected=corrected, approach="Frequency")
  return = list(TABLE=TABLE, CT=CT/100, TCI=TCI, TO=TO, FROM=FROM,
                NET=NET, NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
}
