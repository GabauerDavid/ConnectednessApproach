
#' @title Univariate GARCH selection criterion
#' @description This function estimates and evaluates a combination of GARCH models with different distributions and suggests the best GARCH models among all alternatives given some test statistics
#' @param x zoo data matrix
#' @param distributions Vector of distributions
#' @param models Vector of GARCH models
#' @param ar AR(p)
#' @param ma MA(q)
#' @param prob The quantile (coverage) used for the VaR.
#' @param conf.level Confidence level of VaR test statistics
#' @param lag Lag length of weighted Portmanteau statistics
#' @return Get optimal univariate GARCH model specification
#' @importFrom stats bartlett.test coef fitted fligner.test integrate qnorm quantile residuals sd sigma var.test
#' @references
#' Ghalanos, A. (2014). rugarch: Univariate GARCH models, R package version 1.3-3.
#' 
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.
#' @author David Gabauer
#' @export
GARCHselection = function(x, distributions=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH"), prob=0.05, conf.level=0.90, lag=20, ar=0, ma=0) {
  message("A dynamic version of the optimal univariate GARCH selection procedure is implemented according to:\n Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.")
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  GARCH_IC = matrix(Inf, nrow=length(distributions), ncol=length(models))
  colnames(GARCH_IC) = models
  rownames(GARCH_IC) = distributions
  spec_list = list()
  table_list = list()
  for (i in 1:length(models)) {
    spec_list[[i]] = list()
    table_list[[i]] = list()
    for (j in 1:length(distributions)) {
      spec_list[[i]][[j]] = list()
    }
    names(spec_list[[i]]) = distributions
  }
  names(spec_list) = names(table_list) = models
  
  for (j in 1:length(models)) {
    message(paste0("-",models[j]))
    for (i in 1:length(distributions)) {
      message(paste0("--",distributions[i]))
      if (models[j] %in% c("AVGARCH","TGARCH","APARCH","NAGARCH","NGARCH","ALLGARCH")) {
        ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                          variance.model=list(model="fGARCH", submodel=models[j], garchOrder=c(1,1)),
                                          distribution.model=distributions[i])
      } else {
        ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                          variance.model=list(model=models[j], garchOrder=c(1,1)),
                                          distribution.model=distributions[i])
      }
      ugarch.fit = rugarch::ugarchfit(ugarch.spec, data=x, solver="hybrid", solver.list=list(outer.iter=10, inner.iter=1000, eval.se=FALSE, tol=1e-12))
      if (ugarch.fit@fit$convergence==0) {
        fit = GARCHtests(ugarch.fit, prob=prob, conf.level=conf.level, lag=lag)
        GARCH_IC[i,j] = fit$InformationCriterion
        spec_list[[models[j]]][[distributions[i]]] = ugarch.spec
        table_list[[j]][[distributions[i]]] = fit
      }
    }
  }
  GARCH_selection = which(GARCH_IC==min(GARCH_IC),arr.ind=TRUE)
  best_ugarch = spec_list[[GARCH_selection[2]]][[GARCH_selection[1]]]
  best_table = table_list[[GARCH_selection[2]]][[GARCH_selection[1]]]
  return = list(best_ugarch=best_ugarch, best_table=best_table, GARCH_IC=GARCH_IC, spec_list=spec_list, table_list=table_list)
}
