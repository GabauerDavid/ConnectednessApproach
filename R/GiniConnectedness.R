
#' @title Gini-Based Connectedness Measures
#' @description This function results in Gini-based connectedness measures
#' @param dca Dynamic connectedness object
#' @return Get Gini-based connectedness measures
#' @examples 
#' \donttest{
#' data("dy2012")
#' dca = ConnectednessApproach(dy2012,
#'                             nlag=1,
#'                             nfore=20,
#'                             model="VAR",
#'                             connectedness="Time",
#'                             corrected=TRUE)
#' gca = GiniConnectedness(dca)
#' }
#' @references Gabauer, D., Stenfors, A., & Vinco, V. (2026). Eurozone sectoral inflation networks. International Review of Economics and Finance, 109, 105405.
#' @author David Gabauer
#' @export
GiniConnectedness = function(dca) {
  Gini = function (x) {
    x = as.numeric(na.omit(x))
    n = length(x)
    x = sort(x)
    G = sum(x * 1L:n)
    G = 2 * G/sum(x) - (n + 1L)
    G/n
  }
  CT_iq = CT_eq = dca$CT
  for (idx in 1:dim(dca$CT)[3]) {
    for (jdx in 1:dim(dca$CT)[2]) {
      g = Gini(dca$CT[jdx,,idx])
      CT_eq[jdx,,idx] = (1-g)*dca$CT[jdx,,idx]
      CT_iq[jdx,,idx] = g*dca$CT[jdx,,idx]
    }
  }
  return(list(dca_eq=TimeConnectedness(FEVD=CT_eq), dca_iq=TimeConnectedness(FEVD=CT_iq)))
}
