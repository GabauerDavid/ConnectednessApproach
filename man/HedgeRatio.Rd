% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HedgeRatio.R
\name{HedgeRatio}
\alias{HedgeRatio}
\title{Kroner and Sultan (1993) hedge ratios}
\usage{
HedgeRatio(
  x,
  H,
  method = c("cumsum", "cumprod"),
  statistics = c("Fisher", "Bartlett", "Fligner-Killeen", "Levene", "Brown-Forsythe"),
  metric = "StdDev",
  digit = 2
)
}
\arguments{
\item{x}{zoo return matrix (in percentage)}

\item{H}{Residual variance-covariance, correlation or pairwise connectedness matrix}

\item{method}{Cumulative sum or cumulative product}

\item{statistics}{Hedging effectiveness statistic}

\item{metric}{Risk measure of Sharpe Ratio (StdDev, VaR, or CVaR)}

\item{digit}{Number of decimal places}
}
\value{
Get hedge ratios
}
\description{
This function calculates the hedge ratios of Kroner and Sultan (1993)
}
\examples{
\donttest{
data("g2020")
fit = VAR(g2020, configuration=list(nlag=1))
hr = HedgeRatio(g2020/100, fit$Q)
hr$TABLE
}
}
\references{
Kroner, K. F., & Sultan, J. (1993). Time-varying distributions and dynamic hedging with foreign currency futures. Journal of Financial and Quantitative Analysis, 28(4), 535-551.

Ederington, L. H. (1979). The hedging performance of the new futures markets. The Journal of Finance, 34(1), 157-170.

Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
}
\author{
David Gabauer
}
