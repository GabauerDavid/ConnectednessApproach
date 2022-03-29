
#' @title Minimum connectedness portfolio
#' @description This function calculates the minimum connectedness portfolio
#' @param x zoo return matrix (in percentage)
#' @param H Pairwise connectedness matrix or alternatively variance-covariance or correlation matrix
#' @param method Cumulative sum or cumulative product 
#' @param long Allow only long portfolio position
#' @param statistics Hedging effectiveness statistic
#' @param digit Number of decimal places
#' @return Get portfolio weights
#' @examples
#' data("g2020")
#' fit = VAR(g2020, configuration=list(nlag=1))
#' dca = TimeConnectedness(Phi=fit$B, Sigma=fit$Q, nfore=10, generalized=TRUE)
#' mcp = MinimumConnectednessPortfolio(g2020, dca$PCI, statistics="Fisher")
#' mcp$TABLE
#' @references
#' Broadstock, D. C., Chatziantoniou, I., & Gabauer, D. (2020). Minimum Connectedness Portfolios and the Market for Green Bonds: Advocating Socially Responsible Investment (SRI) Activity. Available at SSRN 3793771.\\
#' Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.\\
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
MinimumConnectednessPortfolio = function(x, H, method=c("cumsum","cumprod"), long=TRUE, statistics=c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"), digit=2) {
  message("The minimum connectedness portfolio is implemented according to:\n Broadstock, D. C., Chatziantoniou, I., & Gabauer, D. (2020). Minimum Connectedness Portfolios and the Market for Green Bonds: Advocating Socially Responsible Investment (SRI) Activity. Available at SSRN 3793771.")
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
  I = matrix(1, k, 1)
  if (dim(H)[[3]]==1) {
    H = array(H, c(k,k,t), dimnames=list(NAMES,NAMES,date))
  }

  portfolio_weights = array(NA,c(t,k), dimnames=list(date, NAMES))
  for (i in 1:t) {
    V_inv = MASS::ginv(H[,,i])
    pw = (V_inv%*%I) / c(t(I)%*%V_inv%*%I)
    if (long) {
      pw = ifelse(pw<0,0,pw)
      pw = ifelse(pw>1,1,pw)
      pw = pw/sum(pw)
    }
    portfolio_weights[i,] = pw
  }
  
  summary = NULL
  for (i in 1:k) {
    x_ = as.matrix(portfolio_weights[,i])
    summary_ = matrix(NA, nrow=ncol(x_), ncol=4)
    for (ij in 1:ncol(x_)){
      summary_[ij,] = matrix(c(mean(x_[,ij]), stats::sd(x_[,ij]), stats::quantile(x_[,ij],0.05), stats::quantile(x_[,ij],0.95)), nrow=1)
    }
    colnames(summary_) = c("Mean", "Std.Dev.", "5%", "95%")
    summary = rbind(summary, summary_)
  }
  rownames(summary) = NAMES
  
  portfolio_return = array(NA,c(t,1), dimnames=list(date))
  for (i in 1:t) {
    portfolio_return[i,] = sum(portfolio_weights[i,]*as.numeric(x[i,]))
  }
  
  if (method=="cumsum") {
    cumulative_portfolio_return = cumsum(portfolio_return)
  } else if (method=="cumprod") {
    cumulative_portfolio_return = cumprod(1+portfolio_return)-1
  }
  
  HE = pvalue = array(NA,c(k,1), dimnames=list(NAMES))
  for (i in 1:k) {
    HE[i,] = 1 - var(portfolio_return)/var(x[,i])
    df = rbind(data.frame(val=portfolio_return, group="A"), data.frame(val=x[,i], group="B"))
    if (statistics=="Fisher") {
      pvalue[i,] = var.test(x=portfolio_return, y=x[,i],ratio=1)$p.value
    } else if (statistics=="Bartlett") {
      pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Bartlett", verbose=F)$p.value
    } else if (statistics=="Fligner-Killeen") {
      pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Fligner", verbose=F)$p.value
    } else if (statistics=="Levene") {
      pvalue[i,] = onewaytests::homog.test(val~as.character(group), data=df, method="Levene", verbose=F)$p.value
    } else if (statistics=="Brown-Forsythe") {
      pvalue[i,] = onewaytests::bf.test(val~as.character(group), data=df, verbose=F)$p.value
    } else {
      stop("No valid hedging effectiveness statistics have been chosen.")
    }
  }
  
  TABLE = cbind(summary,HE,pvalue)
  colnames(TABLE)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
  
  return = list(TABLE=format(round(TABLE,digit),nsmall=digit), portfolio_weights=portfolio_weights, HE=HE, pvalue=pvalue, portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return)
}
