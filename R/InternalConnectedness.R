
#' @title Internal Connectedness Measures
#' @description This function provides internal connectedness measures
#' @param dca Dynamic connectedness object
#' @param groups List of at least two group vectors
#' @param start Start index
#' @param end End index 
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("gg2018")
#' dca = ConnectednessApproach(gg2018, 
#'                             nlag=1, 
#'                             nfore=10, 
#'                             window.size=200,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.99, 
#'                             prior="BayesPrior")))
#' int = InternalConnectedness(dca, groups=list("US"=c(1,2,3,4), "JP"=c(5,6,7,8)))
#' }
#' @references Gabauer, D., & Gupta, R. (2018). On the transmission mechanism of country-specific and international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness decomposition approach. Economics Letters, 171, 63-71.
#' @author David Gabauer
#' @export
InternalConnectedness = function(dca, groups=list(c(1), c(2:ncol(dca$NET))), start=NULL, end=NULL) {
  corrected = dca$config$corrected
  if (is.null(start)) {
    start = 1
  }
  if (is.null(end)) {
    end = dim(dca$CT)[3]
  }
  
  if (length(groups)<=1) {
    stop("groups need to consist of at least 2 vectors")
  }
  if (dca$config$approach=="Joint") {
    stop(paste("Decomposed connectedness measures are not implemented for",dca$approach, "connectedness"))
  } else if (dca$config$approach=="Frequency") {
    ct = dca$CT[,,start:end,]
    NAMES = colnames(ct)
    mn = dim(CT)[4]
    TABLE = list()
    horizons = dimnames(CT)[[4]]
    k = dim(ct)[2]
    ct_inter = ct_wo = ct
    date = as.character(dimnames(ct)[[3]])
    t = dim(ct)[3]
    
    m = length(groups)
    NAMES_group = names(groups)
    if (is.null(NAMES_group)) {
      NAMES_group = paste0("GROUP", 1:m)
    }
    
    for (jl in 1:mn) {
      for (i in 1:m) {
        for (j in 1:m) {
          if (i>j) {
            group_1 = groups[[i]]
            group_2 = groups[[j]]
            ct_wo[group_1,group_2,,jl] = 0
            ct_wo[group_2,group_1,,jl] = 0
          }
        }
      }
    }
    
    TCI_wo = array(0, c(t,mn), dimnames=list(date,horizons))
    FROM_wo = TO_wo = NPT_wo = NET_wo = array(0, c(t,k,mn), dimnames=list(date, NAMES,horizons))
    NPDC_wo = PCI_wo = INFLUENCE_wo = array(0, c(k,k,t,mn), dimnames=list(NAMES, NAMES, date,horizons))
    for (jl in 2:mn) {
      for (i in 1:t) {
        dca_ = ConnectednessTable(ct_wo[,,i,jl])
        TO_wo[i,,jl] = dca_$TO
        FROM_wo[i,,jl] = dca_$FROM
        NET_wo[i,,jl] = dca_$NET
        NPT_wo[i,,jl] = dca_$NPT
        NPDC_wo[,,i,jl] = dca_$NPDC
        if (corrected) {
          TCI_wo[i,jl] = dca_$cTCI
        } else {
          TCI_wo[i,jl] = dca_$TCI
        }
      }
      if (corrected) {
        m_ = (k-1)
      } else {
        m_ = k
      }
      
      TCI_group = array(NA, c(t,m,mn), dimnames=list(date, NAMES_group, horizons))
      for (i in 1:m) {
        group = groups[i][[1]]
        TCI_group[,i,jl] = rowSums(TO_wo[,group,jl,drop=FALSE])/m_
      }
      TABLE[[jl]] = ConnectednessTable(ct_wo[,,,jl])$TABLE
    }
    TABLE[[1]] = ConnectednessTable(ct_wo[,,,1])$TABLE
    names(TABLE) = horizons
    TCI_wo[,1] = apply(TCI_wo,1,sum)
    TO_wo[,,1] = apply(TO_wo,1:2,sum)
    FROM_wo[,,1] = apply(FROM_wo,1:2,sum)
    NET_wo[,,1] = apply(NET_wo,1:2,sum)
    NPDC_wo[,,,1] = apply(NPDC_wo,1:3,sum)
    TCI_group[,i,jl] = rowSums(TO_wo[,group,jl,drop=FALSE])/m_
    for (ij in 1:t) {
      for (jl in length(horizons):1) {
        for (i in 1:k) {
          for (j in 1:k) {
            PCI_wo[i,j,ij,jl] = 200*(ct_wo[i,j,ij,jl]+ct_wo[j,i,ij,jl])/(ct_wo[i,i,ij,1]+ct_wo[i,j,ij,1]+ct_wo[j,i,ij,1]+ct_wo[j,j,ij,1])
          }
        }
        INFLUENCE_wo[,,ij,jl] = abs(NPDC_wo[,,ij,jl]/t(t(ct_wo[,,ij,1])+ct_wo[,,ij,1]))
      }
      NPT_wo[ij,,1] = rowSums(NPDC_wo[,,ij,1]<0)
    }
    ind = which(is.nan(PCI_wo) | is.infinite(PCI_wo),arr.ind=TRUE)
    if (length(ind)>0) {
      PCI_wo[ind] = 0
    }
    ind = which(is.nan(INFLUENCE_wo) | is.infinite(INFLUENCE_wo),arr.ind=TRUE)
    if (length(ind)>0) {
      INFLUENCE_wo[ind] = 0
    }
  } else {
    ct = dca$CT[,,start:end]
    NAMES = colnames(ct)
    k = dim(ct)[2]
    ct_wo = ct
    date = as.character(dimnames(ct)[[3]])
    t = dim(ct)[3]

    m = length(groups)
    NAMES_group = names(groups)
    if (is.null(NAMES_group)) {
      NAMES_group = paste0("GROUP", 1:m)
    }
  
    for (i in 1:m) {
      for (j in 1:m) {
        if (i>j) {
          group_1 = groups[[i]]
          group_2 = groups[[j]]
          ct_wo[group_1,group_2,] = 0
          ct_wo[group_2,group_1,] = 0
        }
      }
    }
  
    TCI_wo = array(NA, c(t, 1), dimnames=list(date, c("TCI")))
    PCI_wo = INFLUENCE_wo = NPDC_wo = array(NA, c(k, k, t), dimnames=list(NAMES,NAMES,date))
    TO_wo = FROM_wo = NET_wo = NPT_wo = array(NA, c(t, k), dimnames=list(date, NAMES))
    for (i in 1:t) {
      dca_ = ConnectednessTable(ct_wo[,,i])
      TO_wo[i,] = dca_$TO
      FROM_wo[i,] = dca_$FROM
      NET_wo[i,] = dca_$NET
      NPT_wo[i,] = dca_$NPT
      NPDC_wo[,,i] = dca_$NPDC
      inf = dca_$INFLUENCE
      inf[which(is.nan(inf),arr.ind=TRUE)] = 0
      INFLUENCE_wo[,,i] = inf
      PCI_wo[,,i] = dca_$PCI
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
    TABLE = ConnectednessTable(ct_wo)$TABLE
  }
  return = list(TABLE=TABLE, gTCI=TCI_group, TCI=TCI_wo, TO=TO_wo, FROM=FROM_wo, 
                NET=NET_wo, NPDC=NPDC_wo, PCI=PCI_wo, INFLUENCE=INFLUENCE_wo, NPT=NPT_wo, config=list(approach="Internal"))
  
}
