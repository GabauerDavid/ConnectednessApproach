% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VFEVD.R
\name{VFEVD}
\alias{VFEVD}
\title{Generalized volatility forecast error variance decomposition and volatility impulse response functions}
\usage{
VFEVD(fit, nfore = 100, standardize = FALSE)
}
\arguments{
\item{fit}{Fitted DCC-GARCH model}

\item{nfore}{H-step ahead forecast horizon}

\item{standardize}{Boolean value whether GIRF should be standardized}
}
\value{
Get volatility impulse response functions and forecast error variance decomposition
}
\description{
This function provides the volatility impulse responses and the forecast error variance decomposition of DCC-GARCH models.
}
\references{
Gabauer, D. (2020). Volatility impulse response analysis for DCC‐GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.
}
\author{
David Gabauer
}
