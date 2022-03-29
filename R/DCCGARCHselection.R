
#' @title DCC-GARCH selection specification
#' @description This function calculates the optimal DCC-GARCH specification
#' @param x zoo data matrix
#' @param distributions Vector of distributions
#' @param models Vector of GARCH models
#' @param ar AR(p)
#' @param ma MA(q)
#' @param prob The quantile (coverage) used for the VaR.
#' @param conf.level Confidence level of VaR test statistics
#' @param lag Lag length of weighted Portmanteau statistics
#' @return Get best DCC-GARCH
#' @importFrom stats bartlett.test coef fitted fligner.test integrate qnorm quantile residuals sd sigma var.test
#' @references
#' Ghalanos, A. (2014). rugarch: Univariate GARCH models, R package version 1.3-3.
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.
#' @author David Gabauer
#' @export
DCCGARCHselection = function(x, distributions=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH"), prob=0.05, conf.level=0.90, lag=20, ar=0, ma=0) {
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  for (i in 1:k) {
    print(colnames(x)[i])
    sel = GARCHselection(x=x[,i], distributions=distributions, models=models, 
                         prob=prob, conf.level=conf.level, lag=lag, ar=ar, ma=ma)
    if (i==1) {
      mspec = sel$best_ugarch
    } else {
      mspec = c(mspec, sel$best_ugarch)
    }
  }
  mspec
}
