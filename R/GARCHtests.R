
#' @title Univariate GARCH test statistics
#' @description This function provides the results of multiple univariate GARCH test statistics
#' @param fit Fitted univariate GARCH
#' @param prob The quantile (coverage) used for the VaR.
#' @param conf.level Confidence level of VaR test statistics
#' @param lag Lag length of weighted Portmanteau statistics
#' @return Get best univariate GARCH
#' @references
#' Ghalanos, A. (2014). rugarch: Univariate GARCH models, R package version 1.3-3.
#' 
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics, 26(1), 1375-1408.
#' @author David Gabauer
#' @export
#' @importFrom rugarch qdist
#' @importFrom rugarch VaRDurTest
#' @importFrom stats qnorm
#' @importFrom utils tail
#' @importFrom xts as.xts
#' @importFrom zoo as.zoo
GARCHtests = function(fit, lag=20, prob=0.05, conf.level=0.90){
  distribution = fit@model$modeldesc$distribution
  model = fit@model$modeldesc$vmodel
  submodel = fit@model$modeldesc$vsubmodel
  ar = fit@model$modelinc['ar']
  ma = fit@model$modelinc['ma']
  if (is.null(submodel)) {
    ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                             variance.model=list(model=model, garchOrder=c(1,1)),
                             distribution.model=distribution)
  } else {
    ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                             variance.model=list(model=model, submodel=submodel, garchOrder=c(1,1)),
                             distribution.model=distribution)
  }
  x = xts::xts(as.numeric(fit@model$modeldata$data), as.Date(fit@model$modeldata$index))
  t = length(x)

  VaR = rugarch::fitted(fit) + rugarch::sigma(fit)*rugarch::qdist(fit@model$modeldesc$distribution, p=prob, mu=fit@fit$matcoef[1,1], sigma=1, skew=ifelse(is.na(rugarch::coef(fit)["skew"]),0,rugarch::coef(fit)["skew"]), shape=rugarch::coef(fit)["shape"])
  var_test = rugarch::VaRTest(prob, actual=x, VaR=as.numeric(VaR), conf.level=conf.level)
  i = 99
  while (is.nan(var_test$uc.LRp)) {
    n = round(t*i/100)
    var_test = rugarch::VaRTest(prob, actual=utils::tail(x,n), VaR=utils::tail(as.numeric(VaR),n), conf.level=conf.level)
    i = i - 1
  }
  vardur_test = rugarch::VaRDurTest(prob, x, VaR,conf.level=conf.level)
  f = function(x) {
    rugarch::qdist(fit@model$modeldesc$distribution, p=x, mu=0,sigma=1,skew=ifelse(is.na(rugarch::coef(fit)["skew"]),0,rugarch::coef(fit)["skew"]), shape=rugarch::coef(fit)["shape"])
  }
  ES = rugarch::fitted(fit) + rugarch::sigma(fit)*stats::integrate(f,0,prob)$value/prob
  ES = rugarch::ESTest(prob, x, ES, VaR, boot=TRUE, n.boot=1000, conf.level=conf.level)
  sign.bias = rugarch::signbias(fit)[1,][1:2]
  warch = WeightedBoxTest(rugarch::residuals(fit), type="Ljung-Box", lag=lag, sqrd.res=TRUE)

  statistics = c(sign.bias[[1]], warch$statistic, var_test$uc.LRstat, vardur_test$rLL, ES$actual.exceed/ES$expected.exceed)
  pvalues = c(sign.bias[2]$prob,warch$p.value, var_test$uc.LRp, round(ES$boot.p.value,3), vardur_test$LRp)
  if (length(which(is.nan(pvalues)))>0) {
    pvalues[which(is.nan(pvalues))] = 1.00
  }
  TABLE = rbind(statistics, pvalues)
  colnames(TABLE) = c("SignBias", paste0("WARCH(",lag,")"),"VaR","CVaR","VaR Dur.")

  qprob = stats::qnorm(1-prob)
  loss = sum(abs(fit@fit$robust.tval[-c(1:2)])<=qprob)
  if (is.na(fit@fit$robust.matcoef[1,2])==FALSE) {
    if ("skew" %in% rownames(fit@fit$robust.matcoef)) {
      upper = fit@fit$robust.matcoef["skew",1] + qprob*fit@fit$robust.matcoef["skew",2]
      lower = fit@fit$robust.matcoef["skew",1] - qprob*fit@fit$robust.matcoef["skew",2]
      if (upper>1 && lower<1) {
        loss = loss + 100
      }
    }
  }
  IC = -2*rugarch::likelihood(fit) + loss*log(t)
  IC = IC + sum(TABLE[2,]<0.10)*10^5
  IC = ifelse(is.na(IC), Inf, IC)
  
  return = list(InformationCriterion=IC, TABLE=TABLE)
}
