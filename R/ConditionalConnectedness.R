
#' @title ConditionalConnectedness
#' @description This function computes the conditional connectedness measures.
#' @param dca Dynamic connectedness object
#' @param group Group vector
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Chatzianzoniou, Gabauer and Stenfors (2022)
#' data(cgs2022)
#' dca = ConnectednessApproach(cgs2022, 
#'                             nlag=1, 
#'                             nfore=10, 
#'                             window.size=250,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.99, 
#'                             prior="BayesPrior")))
#' cc = ConditionalConnectedness(dca, group=c(1,4,7,10,13,16))
#' }
#' @references Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.
#' @author David Gabauer
#' @export
ConditionalConnectedness = function(dca, group=c(1,2,3), start=NULL, end=NULL, corrected=FALSE) {
  message("Conditional connectedness measures are implemented according to:\n Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.")
  
  if (dca$approach=="Frequency" | dca$approach=="Joint") {
    stop(paste("Conditional connectedness measures are not implemented for",dca$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = dca$CT[group,group,start:end,drop=FALSE]
    k = length(group)
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    t = length(date)
    
    TCI = array(NA, c(t,1), dimnames=list(date,"TCI"))
    NPT = TO = FROM = NET = array(NA, c(t,k), dimnames=list(date,NAMES))
    INFLUENCE = PCI = FEVD = NPDC = array(NA, c(k,k,t), dimnames=list(NAMES,NAMES,date))
    for (i in 1:t) {
      cc = ConnectednessTable(ct[,,i]/rowSums(ct[,,i]))
      FEVD[,,i] = cc$FEVD
      TO[i,] = cc$TO
      FROM[i,] = cc$FROM
      NET[i,] = cc$NET
      NPT[i,] = cc$NPT
      NPDC[,,i] = cc$NPDC
      PCI[,,i] = cc$PCI
      INFLUENCE[,,i] = cc$INFLUENCE
      if (corrected) {
        TCI[i,] = cc$cTCI
      } else {
        TCI[i,] = cc$TCI
      }
    }
    TABLE = ConnectednessTable(FEVD/100)$TABLE
    return = list(TABLE=TABLE, FEVD=FEVD, TCI=TCI, NET=NET, TO=TO, FROM=FROM, 
                  NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, approach="Conditional")
  }
}
