
#' @title Aggregated Connectedness Measures
#' @description This function results in aggregated connectedness measures.
#' @param dca Dynamic connectedness object
#' @param groups List of at least two group vectors
#' @param start Start index
#' @param end End index
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' #Replication of Gabauer and Gupta (2018)
#' data("gg2018")
#' dca = ConnectednessApproach(gg2018, 
#'                             nlag=1, 
#'                             nfore=10, 
#'                             window.size=200,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.99, 
#'                             prior="BayesPrior")))
#' ac = AggregatedConnectedness(dca, groups=list("US"=c(1,2,3,4), "JP"=c(5,6,7,8)))
#' }
#' @references Chatziantoniou, I., Gabauer, D., & Stenfor, A. (2021). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves (No. 2021-06). University of Portsmouth, Portsmouth Business School, Economics and Finance Subject Group.
#' @author David Gabauer
#' @export
AggregatedConnectedness = function (dca, groups, start = NULL, end = NULL)  {
  corrected = dca$config$corrected
  message("Aggregated connectedness measures are introduced accoring to:\n Stenfors, A., Chatziantoniou, I., & Gabauer, D. (2022). Independent Policy, Dependent Outcomes: A Game of Cross-Country Dominoes across European Yield Curves. Journal of International Financial Markets, Institutions and Money.")
  if (is.null(start)) {
    start = 1
  }
  if (is.null(end)) {
    end = dim(dca$CT)[3]
  }
  NAMES = dimnames(dca$NET)[[2]]
  k = length(NAMES)
  m = length(groups)
  CT = dca$CT
  t = dim(CT)[3]
  weights = NULL
  for (i in 1:m) {
    weights[i] = length(groups[i][[1]])
  }
  if (is.null(names(groups))) {
    NAMES_group = paste0("GROUP", 1:m)
  }else {
    NAMES_group = names(groups)
  }
  date = as.character(as.Date(dimnames(CT)[[3]]))
  if (length(groups) <= 1) {
    stop("groups need to consist of at least 2 vectors")
  }
  if (dca$config$approach == "Joint") {
    stop(paste("Aggregated connectedness measures are not implemented for", 
               dca$config$approach, "connectedness"))
  }else if (dca$config$approach == "Frequency") {
    mn = dim(CT)[4]
    TABLE = list()
    horizons = dimnames(CT)[[4]]
    TCI_ = TCI = array(0, c(t, mn), dimnames = list(date, 
                                                    horizons))
    FROM = TO = NPT = NET = array(0, c(t, m, mn), dimnames = list(date, 
                                                                  NAMES_group, horizons))
    CT_ = NPDC = PCI = INFLUENCE = array(0, c(m, m, t, mn), 
                                         dimnames = list(NAMES_group, NAMES_group, date, 
                                                         horizons))
    for (jl in 2:length(horizons)) {
      for (il in 1:t) {
        ct0 = ct = CT[, , il, jl]
        for (i in 1:m) {
          for (j in 1:m) {
            x = ct0[groups[i][[1]], groups[j][[1]], drop=FALSE]
            CT_[i, j, il, jl] = sum(rowSums(x)/nrow(x))
          }
        }
        dca_ = ConnectednessTable(CT_[, , il, jl])
        if (corrected) {
          TCI[il, jl] = dca_$cTCI
          TCI_[il, jl] = sum(dca_$TO * (k - weights)/(k - 
                                                        1))
        }else {
          TCI[il, jl] = dca_$TCI
          TCI_[il, jl] = sum(dca_$TO * (k - weights)/k)
        }
        TO[il, , jl] = dca_$TO
        FROM[il, , jl] = dca_$FROM
        NET[il, , jl] = dca_$NET
        NPT[il, , jl] = dca_$NPT
        NPDC[, , il, jl] = dca_$NPDC
      }
      TABLE[[jl]] = ConnectednessTable(CT_[, , , jl])$TABLE
    }
    CT_[, , , 1] = apply(CT_, 1:3, sum)
    TABLE[[1]] = ConnectednessTable(CT_[, , , 1])$TABLE
    names(TABLE) = horizons
    TCI[, 1] = apply(TCI, 1, sum)
    TO[, , 1] = apply(TO, 1:2, sum)
    FROM[, , 1] = apply(FROM, 1:2, sum)
    NET[, , 1] = apply(NET, 1:2, sum)
    NPDC[, , , 1] = apply(NPDC, 1:3, sum)
    for (ij in 1:t) {
      for (jl in length(horizons):1) {
        for (i in 1:m) {
          for (j in 1:m) {
            PCI[i, j, ij, jl] = 200 * (CT_[i, j, ij, 
                                           jl] + CT_[j, i, ij, jl])/(CT_[i, i, ij, 
                                                                         1] + CT_[i, j, ij, 1] + CT_[j, i, ij, 
                                                                                                     1] + CT_[j, j, ij, 1])
          }
        }
        INFLUENCE[, , ij, jl] = abs(NPDC[, , ij, jl]/t(t(CT_[, 
                                                             , ij, 1]) + CT_[, , ij, 1]))
      }
      NPT[ij, , 1] = rowSums(NPDC[, , ij, 1] < 0)
    }
  }else {
    approach = dca$config$approach == "Extended Joint"
    TCI_ = TCI = array(NA, c(t, 1), dimnames = list(date, 
                                                    "TCI"))
    NPT = FROM = TO = NET = array(NA, c(t, m), dimnames = list(date, 
                                                               NAMES_group))
    CT_ = PCI = NPDC = INFLUENCE = array(NA, c(m, m, t), 
                                         dimnames = list(NAMES_group, NAMES_group, date))
    for (il in 1:t) {
      ct0 = ct = CT[, , il]
      for (i in 1:m) {
        for (j in 1:m) {
          x = ct0[groups[i][[1]], groups[j][[1]],drop=F]
          CT_[i, j, il] = sum(rowSums(x)/nrow(x))
        }
      }
      dca_ = ConnectednessTable(CT_[, , il])
      if (corrected) {
        TCI[il, ] = dca_$cTCI
        TCI_[il, ] = sum(dca_$TO * (k - weights)/(k - 
                                                    1))
      }
      else {
        TCI[il, ] = dca_$TCI
        TCI_[il, ] = sum(dca_$TO * (k - weights)/k)
      }
      TO[il, ] = dca_$TO
      FROM[il, ] = dca_$FROM
      NET[il, ] = dca_$NET
      NPT[il, ] = dca_$NPT
      NPDC[, , il] = dca_$NPDC
      PCI[, , il] = dca_$PCI/(2 * approach)
      INFLUENCE[, , il] = dca_$INFLUENCE
    }
    TABLE = ConnectednessTable(CT_)$TABLE
    if (approach) {
      k = dim(NET)[2]
      TABLE[k + 2, k + 1] = "TCI"
      TABLE[k + 3, k + 1] = format(round(mean(TCI), 2), 
                                   nsmall = 2)
    }
  }
  return = list(TABLE = TABLE, TCI_ext = TCI_, TCI = TCI, 
                TO = TO, FROM = FROM, NPT = NPT, NET = NET, NPDC = NPDC, 
                INFLUENCE = INFLUENCE, PCI = PCI, config = list(approach = "Aggregated"))
}