
#' @title Inclusive Connectedness Measures
#' @description This function results in inclusive connectedness measures
#' @param dca Dynamic connectedness object
#' @param group Vector of group indices
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("cegg2022")
#' dca = ConnectednessApproach(cegg2022,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             nlag=1,
#'                             nfore=20,
#'                             corrected=TRUE,
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, 
#'                             kappa2=0.99, prior="BayesPrior")))
#' inc = InclusiveConnectedness(dca, group=c(1,2,3))
#' }
#' @references Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.
#' @author David Gabauer
#' @export
InclusiveConnectedness = function(dca, group=c(1,2), start=NULL, end=NULL) {
  message("The partial connectedness measures are implemented according to:\n Chatziantoniou, I., Elsayed, A., Gabauer, D., & Gozgor, G. (2022). Oil price shocks and exchange rate dynamics: New evidence from decomposed and partial connectedness measures for oil importing and exporting economies.")
  corrected = dca$config$corrected
  if (is.null(start)) {
    start = 1
  }
  if (is.null(end)) {
    end = dim(dca$CT)[3]
  }

  if (dca$config$approach=="Joint") {
    stop(paste("Partial connectedness measures are not implemented for",dca$config$approach, "connectedness"))
  } else if (dca$config$approach=="Frequency") {
    ct = dca$CT[,,start:end,]
    NAMES = dimnames(ct)[[1]]
    date = as.character(dimnames(ct)[[3]])
    k = dim(ct)[1]
    t = dim(ct)[3]
    mn = dim(CT)[4]
    TABLE = list()
    horizons = dimnames(CT)[[4]]
    
    CT = ct*0
    for (ij in 1:mn) {
      for (i in group) {
        CT[,i,,ij] = ct[,i,,ij]
        CT[i,,,ij] = ct[i,,,ij]
      }
    }
    
    TCI = array(0, c(t,mn), dimnames=list(date,horizons))
    FROM = TO = NPT = NET = array(0, c(t,k,mn), dimnames=list(date, NAMES, horizons))
    NPDC = PCI = INFLUENCE = array(0, c(k,k,t,mn), dimnames=list(NAMES, NAMES, date, horizons))
    for (ij in 2:mn) {
      for (i in 1:t) {
        dca_ = ConnectednessTable(CT[,,i,ij])
        NPDC[,,i,ij] = dca_$NPDC
        TO[i,,ij] = dca_$TO
        FROM[i,,ij] = dca_$FROM
        NET[i,,ij] = dca_$NET
        NPT[i,,ij] = dca_$NPT
        if (corrected) {
          TCI[i,ij] = dca_$cTCI
        } else {
          TCI[i,ij] = dca_$TCI
        }
      }
      TABLE[[ij]] = ConnectednessTable(CT[,,,ij])$TABLE
    }
    TABLE[[1]] = ConnectednessTable(CT[,,,1])$TABLE
    names(TABLE) = horizons
    
    TCI[,1] = apply(TCI,1,sum)
    TO[,,1] = apply(TO,1:2,sum)
    FROM[,,1] = apply(FROM,1:2,sum)
    NET[,,1] = apply(NET,1:2,sum)
    NPDC[,,,1] = apply(NPDC,1:3,sum)
    for (ij in 1:t) {
      for (jl in length(horizons):1) {
        for (i in 1:k) {
          for (j in 1:k) {
            PCI[i,j,ij,jl] = 200*(CT[i,j,ij,jl]+CT[j,i,ij,jl])/(CT[i,i,ij,1]+CT[i,j,ij,1]+CT[j,i,ij,1]+CT[j,j,ij,1])
          }
        }
        INFLUENCE[,,ij,jl] = abs(NPDC[,,ij,jl]/t(t(CT[,,ij,1])+CT[,,ij,1]))
      }
      NPT[ij,,1] = rowSums(NPDC[,,ij,1]<0)
    }
    
    ind = which(is.nan(PCI),arr.ind=TRUE)
    if (length(ind)>0) {
      PCI[ind] = 0
    }
    ind = which(is.nan(INFLUENCE),arr.ind=TRUE)
    if (length(ind)>0) {
      INFLUENCE[ind] = 0
    }
    
  } else {
    ct = dca$CT[,,start:end]
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    k = dim(ct)[1]
    t = dim(ct)[3]
    
    CT = ct*0
    for (i in group) {
      CT[,i,] = ct[,i,]
      CT[i,,] = ct[i,,]
    }
    
    TCI = array(NA, c(t,1), dimnames=list(as.character(date), "TCI"))
    NPT = NET = FROM = TO = array(NA, c(t, k), dimnames=list(date, NAMES))
    NPDC = PCI = INFLUENCE = array(NA, c(k, k, t), dimnames=list(NAMES, NAMES, date))
    for (i in 1:t) {
      dca_ = ConnectednessTable(CT[,,i])
      NPDC[,,i] = dca_$NPDC
      PCI[,,i] = dca_$PCI
      INFLUENCE[,,i] = dca_$INFLUENCE
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
  }
  return = list(TABLE=TABLE, TCI=TCI, NET=NET, TO=TO, FROM=FROM, NPT=NPT,
                NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=list(approach="Inclusive"))
}
