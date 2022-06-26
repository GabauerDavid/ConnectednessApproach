
#' @title Kroner and Sultan (1993) hedge ratios
#' @description This function calculates the hedge ratios of Kroner and Sultan (1993)
#' @param x zoo return matrix (in percentage)
#' @param H Residual variance-covariance, correlation or pairwise connectedness matrix
#' @param method Cumulative sum or cumulative product
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @return Get hedge ratios
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' hr = HedgeRatio(g2020, fit$Q)
#' hr$TABLE
#' @references
#' Kroner, K. F., & Sultan, J. (1993). Time-varying distributions and dynamic hedging with foreign currency futures. Journal of Financial and Quantitative Analysis, 28(4), 535-551.
#' 
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.
#' 
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
HedgeRatio = function(x, H, method=c("cumsum","cumprod"), statistics=c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), digit=2) {
  message("Hedge ratios are implemented according to:\n Kroner, K. F., & Sultan, J. (1993). Time-varying distributions and dynamic hedging with foreign currency futures. Journal of Financial and Quantitative Analysis, 28(4), 535-551.")
  message("Hedging effectiveness is calculated according to:\n Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.")
  message("Statistics of the hedging effectiveness measure are implemented according to:\n Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.")
  
  method = match.arg(method)
  statistics = match.arg(statistics)
  x = x / 100
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  k = ncol(x)
  t = nrow(x)
  date = as.character(rownames(x))
  NAMES = colnames(x)

  HR = array(NA, c(k, k, t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:k) {
    for (j in 1:k) {
      HR[i,j,] = H[i,j,] / H[j,j,]
    }
  }

  summary = NULL
  for (i in 1:k) {
    for (j in 1:k) {
      x_ = as.matrix(HR[i,j,])
      summary_ = matrix(NA, nrow=ncol(x_), ncol=4)
      for (ij in 1:ncol(x_)){
        summary_[ij,] = matrix(c(mean(x_[,ij]), stats::sd(x_[,ij]), stats::quantile(x_[,ij],0.05), stats::quantile(x_[,ij],0.95)), nrow=1)
      }
      colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
      rownames(summary_) = paste0(NAMES[i], "/", NAMES[j])
      summary = rbind(summary, summary_)
    }
  }
  colnames(summary) = c("Mean","Std.Dev.","5%","95%")
  
  statistics = "Bartlett"
  HE = pvalue = array(NA,c(k,k), dimnames=list(NAMES,NAMES))
  portfolio_return = cumulative_portfolio_return = array(NA,c(k,k,t), dimnames=list(NAMES,NAMES,date))
  for (i in 1:k) {
    for (j in 1:k) {
      portfolio_return[i,j,] = as.numeric(x[,i] - HR[i,j,]*x[,j])
      HE[j,i] = 1 - var(portfolio_return[i,j,])/var(x[,i])
      df = rbind(data.frame(val=portfolio_return[i,j,], group="A"), data.frame(val=x[,i], group="B"))
      if (statistics=="Fisher") {
        pvalue[i,] = VarianceTest(val~as.character(group), data=df, method="Fisher")$p.value
      } else if (statistics=="Bartlett") {
        pvalue[i,] = VarianceTest(val~as.character(group), data=df, method="Bartlett")$p.value
      } else if (statistics=="Fligner-Killeen") {
        pvalue[i,] = VarianceTest(val~as.character(group), data=df, method="Fligner-Killeen")$p.value
      } else if (statistics=="Levene") {
        pvalue[i,] = VarianceTest(val~as.character(group), data=df, method="Levene")$p.value
      } else if (statistics=="Brown-Forsythe") {
        pvalue[i,] = VarianceTest(val~as.character(group), data=df, method="Brown-Forsythe")$p.value
      } else {
        stop("No valid hedging effectiveness statistics have been chosen.")
      }

      if (method=="cumsum") {
        cumulative_portfolio_return[i,j,] = cumsum(portfolio_return[i,j,])
      } else if (method=="cumprod") {
        cumulative_portfolio_return[i,j,] = cumprod(1+portfolio_return[i,j,])-1
      }
    }
  }
  TABLE = cbind(summary,c(HE),c(pvalue))
  TABLE = TABLE[-which(TABLE[,1]==1),]
  colnames(TABLE)=c("Mean","Std.Dev.","5%","95%","HE","p-value")

  return = list(TABLE=format(round(TABLE,digit),nsmall=digit), hedge_ratio=HR, portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return)
}
