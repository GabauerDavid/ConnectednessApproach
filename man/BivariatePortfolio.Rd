% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BivariatePortfolio.R
\name{BivariatePortfolio}
\alias{BivariatePortfolio}
\title{Kroner and Ng (1998) optimal bivariate portfolio weights}
\usage{
BivariatePortfolio(
  x,
  H,
  method = c("cumsum", "cumprod"),
  long = TRUE,
  statistics = c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"),
  metric = "StdDev",
  digit = 2
)
}
\arguments{
\item{x}{zoo return matrix (in percentage)}

\item{H}{Residual variance-covariance, correlation or pairwise connectedness matrix}

\item{method}{Cumulative sum or cumulative product}

\item{long}{Allow only long portfolio position}

\item{statistics}{Hedging effectiveness statistic}

\item{metric}{Risk measure of Sharpe Ratio (StdDev, VaR, or CVaR)}

\item{digit}{Number of decimal places}
}
\value{
Get bivariate portfolio weights
}
\description{
This function calculates the optimal portfolio weights according to Kroner and Ng (1998)
}
\examples{
\donttest{
data("g2020")
fit = VAR(g2020, configuration=list(nlag=1))
bpw = BivariatePortfolio(g2020/100, fit$Q, method="cumsum", statistics="Fisher")
bpw$TABLE
}
}
\references{
Kroner, K. F., & Ng, V. K. (1998). Modeling asymmetric comovements of asset returns. The Review of Financial Studies, 11(4), 817-844.

Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.

Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
}
\author{
David Gabauer
}
