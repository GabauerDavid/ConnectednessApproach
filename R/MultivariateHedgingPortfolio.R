
#' @title Multivariate Hedging Portfolio
#' @description This function calculates the hedge ratios of Kroner and Sultan (1993)
#' @param x zoo return matrix (in percentage)
#' @param H Residual variance-covariance, correlation or pairwise connectedness matrix
#' @param method Cumulative sum or cumulative product
#' @param statistics Hedging effectiveness statistic
#' @param metric Risk measure of Sharpe Ratio (StdDev, VaR, or CVaR)
#' @param digit Number of decimal places
#' @return Get hedge ratios
#' @importFrom zoo zoo
#' @importFrom zoo index
#' @importFrom PerformanceAnalytics SharpeRatio
#' @importFrom PerformanceAnalytics Return.annualized
#' @importFrom PerformanceAnalytics StdDev.annualized
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' hr = MultivariateHedgingPortfolio(g2020/100, fit$Q)
#' hr$TABLE
#' @author David Gabauer
#' @export
MultivariateHedgingPortfolio = function (x, H, method = c("cumsum", "cumprod"), statistics = c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), metric="StdDev", digit = 2) {
  method = match.arg(method)
  statistics = match.arg(statistics)
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  t = nrow(x)
  date = as.character(rownames(x))
  NAMES = colnames(x)
  if (length(dim(H)) == 2) {
    H = array(H, c(k, k, 1))
  }
  if (dim(H)[[3]] == 1) {
    H = array(H, c(k, k, t), dimnames = list(NAMES, NAMES, date))
  }
  
  BETA = H*0
  for (i in 1:t) {
    for (j in 1:k) {
      COV = H[,,i]
      covxy = COV[-j,j]
      covxx = COV[-j,-j]
      beta = solve(covxx) %*% covxy
      BETA[j,-j,i] = beta
    }
  }
  
  portfolio_return = cumulative_portfolio_return = array(NA, c(t, k), dimnames = list(date))
  ret = risk = SR = HE = pvalue = array(NA, c(k, 1), dimnames = list(NAMES))
  summary = NULL
  for (i in 1:k) {
    beta = t(BETA[i,,])
    summary_ = cbind(apply(beta,2,mean)[-i],
                     apply(beta,2,sd)[-i],
                     apply(beta,2,quantile,0.05)[-i],
                     apply(beta,2,quantile,0.95)[-i])
    rownames(summary_) = paste0(NAMES[i], "/", NAMES[-i])
    summary = rbind(summary, summary_)
    
    portfolio_return[,i] = x[,i] - rowSums(x * beta)
    HE[i,] = 1 - var(portfolio_return[,i])/var(x[, i])
    z = zoo(portfolio_return[,i], order.by=index(x))
    SR[i,] = SharpeRatio(z, FUN=(metric), annualize=TRUE)
    ret[i,] = Return.annualized(z)
    risk[i,] = StdDev.annualized(z)
    df = rbind(data.frame(val = x[,i], group = "A"), data.frame(val = portfolio_return[,i], group = "B"))
    pvalue[i,] = VarianceTest(val ~ as.character(group), data = df, method = statistics)$p.value

    if (method == "cumsum") {
      cumulative_portfolio_return[,i] = cumsum(portfolio_return[,i])
    } else if (method == "cumprod") {
      cumulative_portfolio_return[,i] = cumprod(1 + portfolio_return[,i]) - 1
    }
  }
  TABLE = cbind(summary, rep(HE,each=k-1), rep(pvalue,each=k-1), rep(ret,each=k-1), rep(risk,each=k-1), rep(SR,each=k-1))
  colnames(TABLE) = c("Mean", "Std.Dev.", "5%", "95%", "HE", 
                      "p-value", "Return", "Risk","SR")
  return = list(TABLE = format(round(TABLE, digit), nsmall = digit), Beta=BETA,
                portfolio_return = portfolio_return, 
                cumulative_portfolio_return = cumulative_portfolio_return)
}
