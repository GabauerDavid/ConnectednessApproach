
#' @title WeightedBoxTest
#' @description Weighted portmanteau tests for testing the null hypothesis of adequate ARMA fit and/or for detecting nonlinear processes. Written in the style of Box.test() and is capable of performing the traditional Box Pierce (1970), Ljung Box (1978) or Monti (1994) tests.
#' @param x a numeric vector or univariate time series, or residuals of a fitted time series
#' @param lag the statistic will be based on lag autocorrelation coefficients. lag=1 by default
#' @param type test to be performed, partial matching is used. "Box-Pierce" by default
#' @param fitdf number of degrees of freedom to be subtracted if x is a series of residuals, set at 0 by default
#' @param sqrd.res A flag, should the series/residuals be squared to detect for nonlinear effects?, FALSE by default
#' @param log.sqrd.res A flag, should a log of the squared series/residuals be used to detect for nonlinear effects? FALSE by default
#' @param abs.res A flag, should the absolute series or residuals be used to detect for nonlinear effects? FALSE by default
#' @param weighted A flag determining if the weighting scheme should be utilized. TRUE by default. If set to FALSE, the traditional test is performed with no weights
#' @return Get Uninformative Prior
#' @references Box, G. E. P. and Pierce, D. A. (1970), Distribution of residual correlations in autoregressive-integrated moving average time series models. Journal of the American Statistical Association, 65, 1509-1526.
#' 
#' Fisher, T. J. and Gallagher, C. M. (2012), New Weighted Portmanteau Statistics for Time Series Goodness-of-Fit Testing. Journal of the American Statistical Association, accepted.
#' 
#' Ljung, G. M. and Box, G. E. P. (1978), On a measure of lack of fit in time series models. Biometrika 65, 297-303.
#' 
#' Mahdi, E. and McLeod, A. I. (2012), Improved multivariate portmanteau test. Journal of Time Series Analysis 65(2), 297-303.
#' 
#' Monti, A. C. (1994), A proposal for a residual autocorrelation test in linear models. Biometrika 81(4), 776-780.
#' 
#' Pena, D. and Rodriguez, J. (2002) A powerful portmanteau test of lack of fit for time series. Journal of the American Statistical Association 97(458), 601-610.
#' @author David Gabauer
#' @importFrom stats acf
#' @importFrom stats na.pass
#' @importFrom stats pchisq
#' @importFrom stats pgamma
#' @export
WeightedBoxTest = function(x, lag=1, type=c("Box-Pierce", "Ljung-Box", "Monti"),
                           fitdf=0, sqrd.res=FALSE, log.sqrd.res=FALSE, 
                           abs.res=FALSE, weighted=TRUE) {
  if (NCOL(x) > 1) 
    stop("x is not a vector or univariate time series")
  if (lag < 1) 
    stop("Lag must be positive")
  if (fitdf < 0) 
    stop("Fitdf cannot be negative")
  if (fitdf >= lag) 
    stop("Lag must exceed fitted degrees of freedom")
  if ((sqrd.res && log.sqrd.res) || (sqrd.res && abs.res) || 
      (log.sqrd.res && abs.res)) 
    stop("Only one option of: sqrd.res, log.sqrd.res or abs.res can be selected")
  DNAME <- deparse(substitute(x))
  type <- match.arg(type)
  if (abs.res) {
    x <- abs(x)
  }
  if (sqrd.res || log.sqrd.res) {
    x <- x^2
  }
  if (log.sqrd.res) {
    x <- log(x)
  }
  if (weighted) {
    if (type == "Monti") {
      METHOD <- "Weighted Monti test (Gamma Approximation)"
      cor <- acf(x, lag.max = lag, type = "partial", 
                 plot = FALSE, na.action = na.pass)
      obs <- cor$acf[1:lag]
    }
    else {
      cor <- acf(x, lag.max = lag, type = "correlation", 
                 plot = FALSE, na.action = na.pass)
      obs <- cor$acf[2:(lag + 1)]
    }
    if (type == "Ljung-Box") {
      METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
    }
    n <- sum(!is.na(x))
    weights <- (lag - 1:lag + 1)/(lag)
    if (type == "Box-Pierce") {
      METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
      STATISTIC <- n * sum(weights * obs^2)
    }
    else {
      STATISTIC <- n * (n + 2) * sum(weights * (1/seq.int(n - 
                                                            1, n - lag) * obs^2))
    }
    if (sqrd.res) {
      fitdf <- 0
      names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
    } else if (log.sqrd.res) {
      fitdf <- 0
      names(STATISTIC) <- "Weighted X-squared on Log-Squared Residuals for detecting nonlinear processes"
    } else if (abs.res) {
      fitdf <- 0
      names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes"
    } else {
      names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
    }
    shape <- (3/4) * (lag + 1)^2 * lag/(2 * lag^2 + 3 * lag + 
                                          1 - 6 * lag * fitdf)
    scale <- (2/3) * (2 * lag^2 + 3 * lag + 1 - 6 * lag * 
                        fitdf)/lag/(lag + 1)
    PARAMETER <- c(shape, scale)
    names(PARAMETER) <- c("Shape", "Scale")
    PVAL <- 1 - pgamma(STATISTIC, shape = shape, scale = scale)
    names(PVAL) <- "Approximate p-value"
  } else {
    if (type == "Monti") {
      METHOD <- "Monti test"
      cor <- acf(x, lag.max = lag, type = "partial", 
                 plot = FALSE, na.action = na.pass)
      obs <- cor$acf[1:lag]
    } else {
      cor <- acf(x, lag.max = lag, type = "correlation", 
                 plot = FALSE, na.action = na.pass)
      obs <- cor$acf[2:(lag + 1)]
    }
    if (type == "Ljung-Box") {
      METHOD <- "Ljung-Box test"
    }
    n <- sum(!is.na(x))
    if (type == "Box-Pierce") {
      METHOD <- "Box-Pierce test"
      STATISTIC <- n * sum(obs^2)
    } else {
      STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1, 
                                                n - lag) * obs^2))
    }
    if (sqrd.res) {
      fitdf <- 0
      names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
    } else if (log.sqrd.res) {
      fitdf <- 0
      names(STATISTIC) <- "X-squared on Log-Squared Residuals for detecting nonlinear processes"
    } else if (abs.res) {
      fitdf <- 0
      names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes"
    } else {
      names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
    }
    mydf <- lag - fitdf
    PARAMETER <- c(mydf)
    names(PARAMETER) <- c("df")
    PVAL <- 1 - pchisq(STATISTIC, df = mydf)
    names(PVAL) <- "p-value"
  }
  structure(list(statistic=STATISTIC, parameter=PARAMETER, 
                 p.value=PVAL, method=METHOD, data.name=DNAME))
}
