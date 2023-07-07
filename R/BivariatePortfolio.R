
#' @title Kroner and Ng (1998) optimal bivariate portfolio weights
#' @description This function calculates the optimal portfolio weights according to Kroner and Ng (1998)
#' @param x zoo return matrix (in percentage)
#' @param H Residual variance-covariance, correlation or pairwise connectedness matrix
#' @param method Cumulative sum or cumulative product
#' @param long Allow only long portfolio position
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @param metric Risk measure of Sharpe Ratio (StdDev, VaR, or CVaR)
#' @return Get bivariate portfolio weights
#' @importFrom stats var.test
#' @importFrom zoo zoo
#' @importFrom zoo index
#' @importFrom PerformanceAnalytics SharpeRatio
#' @importFrom rmgarch cgarchspec
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' bpw = BivariatePortfolio(g2020/100, fit$Q, method="cumsum", statistics="Fisher")
#' bpw$TABLE
#' @references
#' Kroner, K. F., & Ng, V. K. (1998). Modeling asymmetric comovements of asset returns. The Review of Financial Studies, 11(4), 817-844.
#' 
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.
#' 
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
BivariatePortfolio = function (x, H, method = c("cumsum", "cumprod"), long = TRUE, statistics = c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), metric="StdDev", digit = 2) {
  message("The optimal bivariate portfolios are computed according to:\n Kroner, K. F., & Ng, V. K. (1998). Modeling asymmetric comovements of asset returns. The Review of Financial Studies, 11(4), 817-844.\n\n          Hedging effectiveness is calculated according to:\n Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.\n\n          Statistics of the hedging effectiveness measure are implemented according to:\n Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.")
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
  summary = NULL
  portfolio_weights = array(0.5, c(k, k, t), dimnames = list(NAMES, NAMES, date))
  for (i in 1:k) {
    for (j in 1:k) {
      pw = (H[j, j, ] - H[i, j, ])/(H[i, i, ] - 2 * H[i, j, ] + H[j, j, ])
      pw[which(is.na(pw))] = 0.5
      if (long) {
        pw = ifelse(pw > 1, 1, pw)
        pw = ifelse(pw < 0, 0, pw)
      }
      portfolio_weights[j, i, ] = pw
      portfolio_weights[i, j, ] = 1 - pw
      x_ = as.matrix(portfolio_weights[j, i, ])
      summary_ = matrix(NA, nrow = ncol(x_), ncol = 4)
      for (ij in 1:ncol(x_)) {
        summary_[ij, ] = matrix(c(mean(x_[, ij]), 
                                  stats::sd(x_[, ij]), 
                                  stats::quantile(x_[, ij], 0.05), 
                                  stats::quantile(x_[, ij], 0.95)), nrow = 1)
      }
      colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
      rownames(summary_) = paste0(NAMES[i], "/", NAMES[j])
      summary = rbind(summary, summary_)
    }
  }
  SR = pvalue = HE = array(NA, c(k, k), dimnames = list(NAMES, NAMES))
  portfolio_return = cumulative_portfolio_return = cumulative_asset_return = array(NA, c(k, k, t), dimnames = list(NAMES, NAMES, date))
  for (i in 1:k) {
    for (j in 1:k) {
      portfolio_return[j, i, ] = portfolio_weights[j, i, ] * x[, i] + (1 - portfolio_weights[j, i,]) * x[, j]
      HE[i, j] = 1 - var(portfolio_return[j, i, ])/var(x[,i])
      SR[i, j] = SR[j, i] = PerformanceAnalytics::SharpeRatio(zoo::zoo(portfolio_return[i, j, ], order.by=index(x)), FUN=(metric), annualize=TRUE)
      df = rbind(data.frame(val=x[, i], group = "A"), data.frame(val=portfolio_return[j,i,], group = "B"))
      pvalue[i, j] = VarianceTest(formula=val ~ as.character(group), data=df, method=statistics)$p.value
      if (method == "cumsum") {
        cumulative_asset_return[j, i, ] = cumsum(x[, i])
        cumulative_portfolio_return[j, i, ] = cumsum(portfolio_return[j, i, ])
      } else if (method == "cumprod") {
        cumulative_asset_return[j, i, ] = cumprod(1 + x[, i])
        cumulative_portfolio_return[j, i, ] = cumprod(1 + portfolio_return[j, i, ])
      }
    }
  }
  TABLE = cbind(summary, c(t(HE)), c(t(pvalue)), c(t(SR)))
  TABLE = TABLE[-which(TABLE[, 1] == 0.5), ]
  colnames(TABLE) = c("Mean", "Std.Dev.", "5%", "95%", "HE", "p-value", "SR")
  return = list(TABLE = format(round(TABLE, digit), nsmall = digit), 
                portfolio_weights = portfolio_weights, portfolio_return = portfolio_return, 
                cumulative_portfolio_return = cumulative_portfolio_return)
}
