
#' @title Equally weighted portfolio
#' @description This function calculates the equality weighted portfolio
#' @param x zoo return matrix (in percentage)
#' @param method Cumulative sum or cumulative product 
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @param metric Risk measure of Sharpe Ratio (StdDev, VaR, or CVaR)
#' @return Get portfolio weights
#' @importFrom zoo zoo
#' @importFrom zoo index
#' @examples
#' data("g2020")
#' mcp = EquallyWeightedPortfolio(g2020/100, statistics="Fisher")
#' mcp$TABLE
#' @references
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.
#' 
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
EquallyWeightedPortfolio = function (x, method = c("cumsum", "cumprod"), statistics = c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), metric="StdDev", digit = 2) {
  message("Hedging effectiveness is calculated according to:\n Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.\n\n          Statistics of the hedging effectiveness measure are implemented according to:\n Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.")
  method = match.arg(method)
  statistics = match.arg(statistics)
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  t = nrow(x)
  date = as.character(rownames(x))
  NAMES = colnames(x)
  portfolio_weights = array(1/k, c(t, k), dimnames = list(date, NAMES))
  
  summary = NULL
  for (i in 1:k) {
    x_ = as.matrix(portfolio_weights[, i])
    summary_ = matrix(NA, nrow = ncol(x_), ncol = 4)
    for (ij in 1:ncol(x_)) {
      summary_[ij, ] = matrix(c(mean(x_[, ij]), stats::sd(x_[, ij]), stats::quantile(x_[, ij], 0.05), stats::quantile(x_[, ij], 0.95)), nrow = 1)
    }
    colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
    summary = rbind(summary, summary_)
  }
  rownames(summary) = NAMES
  portfolio_return = array(NA, c(t, 1), dimnames = list(date))
  for (i in 1:t) {
    portfolio_return[i, ] = sum(portfolio_weights[i, ] * 
                                  as.numeric(x[i, ]))
  }
  if (method == "cumsum") {
    cumulative_portfolio_return = cumsum(portfolio_return)
  } else if (method == "cumprod") {
    cumulative_portfolio_return = cumprod(1 + portfolio_return) - 1
  }
  SR = HE = pvalue = array(NA, c(k, 1), dimnames = list(NAMES))
  for (i in 1:k) {
    HE[i, ] = 1 - var(portfolio_return)/var(x[, i])
    z = zoo::zoo(portfolio_return, order.by=index(x))
    SR[i,] = PerformanceAnalytics::SharpeRatio(z, FUN=(metric), annualize=TRUE)
    df = rbind(data.frame(val = portfolio_return, group = "A"), 
               data.frame(val = x[, i], group = "B"))
    pvalue[i, ] = VarianceTest(val ~ as.character(group), 
                               data = df, method = statistics)$p.value

  }
  TABLE = cbind(summary, HE, pvalue, SR)
  colnames(TABLE) = c("Mean", "Std.Dev.", "5%", 
                      "95%", "HE", "p-value", "SR")
  return = list(TABLE = format(round(TABLE, digit), nsmall = digit), 
                portfolio_weights = portfolio_weights, HE = HE, pvalue = pvalue, 
                portfolio_return = portfolio_return, cumulative_portfolio_return = cumulative_portfolio_return)
}
