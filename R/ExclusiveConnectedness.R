
#' @title Exclusive Connectedness Measures
#' @description This function results in exclusive connectedness measures
#' @param dca Dynamic connectedness object
#' @param group Vector of group indices
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Chatziantoniou, et al. (2022)
#' data("cegg2022")
#' dca = ConnectednessApproach(cegg2022,
#'                             nlag=1,
#'                             nfore=20,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             corrected=TRUE,
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, 
#'                             kappa2=0.99, prior="BayesPrior")))
#' exc = ExclusiveConnectedness(dca, group=c(1,2,3))
#' }
#' @references Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.
#' @author David Gabauer
#' @export
ExclusiveConnectedness = function(dca, group=c(1,2), start=NULL, end=NULL) {
  message("Partial connectedness measures are implemented according to:\n Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.")
  corrected = dca$config$corrected
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Partial connectedness measures are not implemented for",dca$config$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = dca$CT[,,start:end]
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    k = dim(ct)[1]
    t = dim(ct)[3]
    
    CT = ct
    for (i in group) {
      CT[,i,] = ct[,i,]*0
      CT[i,,] = ct[i,,]*0
    }
  
    TCI = array(NA, c(t,1), dimnames=list(as.character(date), "TCI"))
    NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
    NPDC = PCI = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
    for (i in 1:t) {
      dca_ = ConnectednessTable(CT[,,i])
      NPDC[,,i] = dca_$NPDC
      PCI[,,i] = dca_$PCI
      infl = dca_$INFLUENCE
      infl[which(is.nan(infl), arr.ind=TRUE)] = 0
      INFLUENCE[,,i] = infl
      TO[i,] = dca_$TO
      FROM[i,] = dca_$FROM
      NET[i,] = dca_$NET
      NPT[i,] = dca_$NPT
      if (corrected) {
        TCI[i,] = dca_$cTCI
      } else {
        TCI[i,] = dca_$TCI
      }
    }
    ind = which(is.nan(PCI),arr.ind=TRUE)
    if (length(ind)>0) {
      PCI[ind] = 0
    }
    ind = which(is.nan(INFLUENCE),arr.ind=TRUE)
    if (length(ind)>0) {
      INFLUENCE[ind] = 0
    }
    TABLE = ConnectednessTable(CT)$TABLE
    config = list(approach="Exclusive")
    return = list(TABLE=TABLE, TCI=TCI, NET=NET, TO=TO, FROM=FROM, NPT=NPT,
                  NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=config)
  }
}
