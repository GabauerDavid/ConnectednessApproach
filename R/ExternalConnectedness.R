
#' @title External Connectedness Measures
#' @description This function provides external connectedness measures
#' @param dca Dynamic connectedness object
#' @param groups List of at least two group vectors
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("gg2018")
#' dca = ConnectednessApproach(gg2018, model="TVP-VAR",
#'                             connectedness="Time",
#'                             nlag=1, nfore=10, window.size=200,
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, 
#'                             kappa2=0.99, prior="BayesPrior")))
#' ext = ExternalConnectedness(dca, groups=list("US"=c(1,2,3,4), "JP"=c(5,6,7,8)))
#' }
#' @references Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.
#' @author David Gabauer
#' @export
ExternalConnectedness = function(dca, groups=list(c(1), c(2:ncol(dca$NET))), start=NULL, end=NULL) {
  message("The decomposed connectedness measures are implemented according to:\n Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.")
  corrected = dca$config$corrected
  if (length(groups)<=1) {
    stop("groups need to consist of at least 2 vectors")
  }
  if (dca$config$approach=="Frequency" | dca$config$approach=="Joint") {
    stop(paste("Decomposed connectedness measures are not implemented for",dca$approach, "connectedness"))
  } else {
    if (is.null(start)) {
      start = 1
    }
    if (is.null(end)) {
      end = dim(dca$CT)[3]
    }
    ct = 100*dca$CT[,,start:end]
    NAMES = colnames(ct)
    k = dim(ct)[2]
    if (length(dim(ct))==2) {
      ct = array(ct, c(k,k,1),dimnames=list(NAMES,NAMES))
    }
    ct_inter = ct_wo = ct
    date = as.character(dimnames(ct)[[3]])
    t = dim(ct)[3]
    
    m = length(groups)
    NAMES_group = names(groups)
    if (is.null(NAMES_group)) {
      NAMES_group = paste0("GROUP", 1:m)
    }
    
    for (i in 1:m) {
      group_1 = groups[[i]]
      ct_wo[group_1,group_1,] = 0
    }
    
    TCI_wo = array(NA, c(t, 1), dimnames=list(date, c("TCI")))
    INFLUENCE_wo = PCI_wo = NPDC_wo = array(NA, c(k, k, t), dimnames=list(NAMES,NAMES,date))
    TO_wo = FROM_wo = NET_wo = NPT_wo = array(NA, c(t, k), dimnames=list(date, NAMES))
    for (i in 1:t) {
      dca_ = ConnectednessTable(ct_wo[,,i]/100)
      TO_wo[i,] = dca_$TO
      FROM_wo[i,] = dca_$FROM
      NET_wo[i,] = dca_$NET
      NPT_wo[i,] = dca_$NPT
      NPDC_wo[,,i] = dca_$NPDC
      PCI_wo[,,i] = dca_$PCI
      infl = dca_$INFLUENCE
      infl[which(is.nan(infl), arr.ind=TRUE)] = 0
      INFLUENCE_wo[,,i] = infl
      if (corrected) {
        TCI_wo[i,] = dca_$cTCI
      } else {
        TCI_wo[i,] = dca_$TCI
      }
    }
    if (corrected) {
      m_ = (k-1)
    } else {
      m_ = k
    }
    
    TCI_group = array(NA, c(t,m), dimnames=list(date, NAMES_group))
    for (i in 1:m) {
      group = groups[i][[1]]
      TCI_group[,i] = rowSums(TO_wo[,group,drop=FALSE])/m_
    }
    
    TABLE = ConnectednessTable(ct_wo/100)$TABLE
    config = list(approach="External")
    return = list(TABLE=TABLE, gTCI=TCI_group, TCI=TCI_wo, TO=TO_wo, FROM=FROM_wo, NPT=NPT_wo,
                  NET=NET_wo, NPDC=NPDC_wo, PCI=PCI_wo, INFLUENCE=INFLUENCE_wo, config=config)
  }
}
