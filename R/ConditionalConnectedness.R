
#' @title ConditionalConnectedness
#' @description This function computes the conditional connectedness measures.
#' @param dca Dynamic connectedness object
#' @param group Group vector
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
ConditionalConnectedness = function(dca, group=c(1,2,3), start=NULL, end=NULL) {
  corrected = dca$config$corrected
  message("Conditional connectedness measures are implemented according to:\n Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.")
  k = length(group)
  if (is.null(start)) {
    start = 1
  }
  if (is.null(end)) {
    end = dim(dca$CT)[3]
  }
  
  if (dca$config$approach=="Joint") {
    stop(paste("Conditional connectedness measures are not implemented for",dca$config$approach, "connectedness"))
  } else if (dca$config$approach=="Frequency" ) {
    ct = dca$CT[group,group,start:end,,drop=FALSE]
    NAMES = dimnames(ct)[[1]]
    date = dimnames(ct)[[3]]
    t = length(date)
    mn = dim(CT)[4]
    TABLE = list()
    horizons = dimnames(CT)[[4]]
    
    TCI = array(0, c(t,mn), dimnames=list(date,horizons))
    FROM = TO = NPT = NET = array(0, c(t,k,mn), dimnames=list(date, NAMES, horizons))
    CT_ = NPDC = PCI = INFLUENCE = array(0, c(k,k,t,mn), dimnames=list(NAMES, NAMES, date, horizons))
    for (jl in 2:mn) {
      for (i in 1:t) {
        cc = ConnectednessTable(ct[,,i,jl]/rowSums(ct[,,i,1]))
        CT_[,,i,jl] = cc$FEVD/100
        TO[i,,jl] = cc$TO
        FROM[i,,jl] = cc$FROM
        NET[i,,jl] = cc$NET
        NPT[i,,jl] = cc$NPT
        NPDC[,,i,jl] = cc$NPDC
        if (corrected) {
          TCI[i,jl] = cc$cTCI
        } else {
          TCI[i,jl] = cc$TCI
        }
      }
      TABLE[[jl]] = ConnectednessTable(CT_[,,,jl])$TABLE
    }
    CT_[,,,1] = apply(CT_,1:3,sum)
    TABLE[[1]] = ConnectednessTable(CT_[,,,1])$TABLE
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
            PCI[i,j,ij,jl] = 200*(CT_[i,j,ij,jl]+CT_[j,i,ij,jl])/(CT_[i,i,ij,1]+CT_[i,j,ij,1]+CT_[j,i,ij,1]+CT_[j,j,ij,1])
          }
        }
        INFLUENCE[,,ij,jl] = abs(NPDC[,,ij,jl]/t(t(CT_[,,ij,1])+CT_[,,ij,1]))
      }
      NPT[ij,,1] = rowSums(NPDC[,,ij,1]<0)
    }
  } else {
    approach = dca$config$approach=="Extended Joint"
    ct = dca$CT[group,group,start:end,drop=FALSE]
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
      PCI[,,i] = cc$PCI / (2*approach)
      INFLUENCE[,,i] = cc$INFLUENCE
      if (corrected) {
        TCI[i,] = cc$cTCI
      } else {
        TCI[i,] = cc$TCI
      }
    }
    TABLE = ConnectednessTable(FEVD/100)$TABLE
    if (approach) {
      k = dim(NET)[2]
      TABLE[k+2,k+1] = "TCI"
      TABLE[k+3,k+1] = format(round(mean(TCI),2), nsmall=2)
    }
  }
  return = list(TABLE=TABLE, FEVD=FEVD, TCI=TCI, NET=NET, TO=TO, FROM=FROM, 
                NPT=NPT, NPDC=NPDC, PCI=PCI, INFLUENCE=INFLUENCE, config=list(approach="Conditional"))
  
}
