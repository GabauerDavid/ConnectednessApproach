#' @title Connectedness Approach
#' @description This function provides a modular framework combining various models and connectedness frameworks.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param nfore H-step ahead forecast horizon
#' @param window.size Rolling-window size or Bayes Prior sample size
#' @param model Estimation model
#' @param corrected Boolean value whether corrected or standard TCI should be computed
#' @param connectedness Type of connectedness approach
#' @param VAR_config Config for VAR model
#' @param DCC_config Config for DCC-GARCH model
#' @param Connectedness_config Config for connectedness approach
#' @return Get connectedness measures
#' @examples
#' \donttest{
#' data("acg2020")
#' dca = ConnectednessApproach(acg2020, 
#'                             nlag=1, 
#'                             nfore=12,
#'                             model="TVP-VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.96,
#'                                             prior="MinnesotaPrior", gamma=0.1)))
#' dca$TABLE
#' }
#' @import progress
#' @importFrom rugarch ugarchspec
#' @importFrom rugarch multispec
#' @importFrom rmgarch dccspec
#' @importFrom rmgarch dccfit
#' @references
#' Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal, 119(534), 158-171.
#' 
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting, 28(1), 57-66.
#' 
#' Barunik, J., & Krehlik, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics, 16(2), 271-296.
#' 
#' Gabauer, D. (2020). Volatility impulse response analysis for DCC-GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.
#' 
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management, 13(4), 84.
#' 
#' Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.
#' 
#' Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy, 73, 102219.
#' 
#' Chatziantoniou, I., & Gabauer, D. (2021). EMU risk-synchronisation and financial fragility through the prism of dynamic connectedness. The Quarterly Review of Economics and Finance, 79, 1-14.
#' 
#' Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters, 204, 109891.
#' 
#' Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management, 60, 100680.
#' 
#' Gabauer, D., Gupta, R., Marfatia, H., & Miller, S. (2020). Estimating US Housing Price Network Connectedness: Evidence from Dynamic Elastic Net, Lasso, and Ridge Vector Autoregressive Models (No. 202065). University of Pretoria, Department of Economics.
#' 
#' Chatziantoniou, I., Gabauer, D., & Gupta, R. (2021). Integration and Risk Transmission in the Market for Crude Oil: A Time-Varying Parameter Frequency Connectedness Approach (No. 202147).
#' 
#' Chatziantoniou, I., Aikins Abakah, E. J., Gabauer, D., & Tiwari, A. K. (2022). Quantile time-frequency price connectedness between green bond, green equity, sustainable investments and clean energy markets. Journal of Cleaner Production.
#' 
#' Cunado, J, Chatziantoniou, I., Gabauer, D., Hardik, M., & de Garcia, F.P. (2022). Dynamic spillovers across precious metals and energy realized volatilities: Evidence from quantile extended joint connectedness measures.
#' @author David Gabauer
#' @export
ConnectednessApproach = function(x,
                                 nlag=1, 
                                 nfore=10, 
                                 window.size=NULL, 
                                 corrected=FALSE,
                                 model=c("VAR", "QVAR", "LASSO", "Ridge", "Elastic", "TVP-VAR", "DCC-GARCH"),
                                 connectedness=c("Time","Frequency", "Joint", "Extended Joint"),
                                 VAR_config=list(
                                   QVAR=list(tau=0.5),
                                   ElasticNet=list(nfolds=10, alpha=NULL, loss="mae", delta_alpha=0.1),
                                   TVPVAR=list(kappa1=0.99, kappa2=0.99, prior="BayesPrior", gamma=0.01)),
                                 DCC_config=list(standardize=FALSE),
                                 Connectedness_config = list(
                                   TimeConnectedness=list(generalized=TRUE),
                                   FrequencyConnectedness=list(partition=c(pi,pi/2,0), generalized=TRUE, scenario="ABS")
                                 )) {
  if (class(x)!="zoo") {
    stop("Data needs to be of type 'zoo'")
  }
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (nfore<=0) {
    stop("nfore needs to be a positive integer")
  }
  model = match.arg(model)
  connectedness = match.arg(connectedness)
  
  NAMES = colnames(x)
  k = ncol(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  t = nrow(x)
  if (is.null(window.size)) {
    window.size = nrow(x)
    t0 = 1
  } else {
    window.size = window.size - nlag
    t0 = t - window.size + 1
  }

  if (model=="VAR") {
    var_model = VAR
    configuration = list(nlag=nlag)
  } else if (model=="QVAR") {
    var_model = QVAR
    configuration = list(nlag=nlag, tau=VAR_config$QVAR$tau)
  } else if (model=="LASSO") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=1, nfolds=VAR_config$ElasticNet$nfolds, loss=VAR_config$ElasticNet$loss)
  } else if (model=="Ridge") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=0, nfolds=VAR_config$ElasticNet$nfolds, loss=VAR_config$ElasticNet$loss)
  } else if (model=="Elastic") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=VAR_config$ElasticNet$alpha, nfolds=VAR_config$ElasticNet$nfolds,
                         loss=VAR_config$ElasticNet$loss, delta_alpha=VAR_config$ElasticNet$delta_alpha)
  } else if (model=="TVP-VAR") {
    prior_ = VAR_config$TVPVAR$prior
    if (prior_=="MinnesotaPrior") {
      prior = MinnesotaPrior(gamma=VAR_config$TVPVAR$gamma, k=k, nlag=nlag)
    } else if (prior_=="UninformativePrior") {
      prior = UninformativePrior(k=k, nlag=nlag)
    } else if (prior_=="BayesPrior") {
      prior = BayesPrior(x=x, size=window.size, nlag=nlag)
    } else {
      stop("Error Prior choice")
    }
    configuration = list(l=c(VAR_config$TVPVAR$kappa1, VAR_config$TVPVAR$kappa2), nlag=nlag, prior=prior)
    var_model = TVPVAR
  } else if (model=="DCC-GARCH") {
    ugarch.spec = rugarch::ugarchspec(mean.model=list(armaOrder=c(0,0)),
                             variance.model=list(garchOrder=c(1,1), model="sGARCH"),
                             distribution.model="norm")
    mgarch.spec = rmgarch::dccspec(uspec=rugarch::multispec(replicate(k, ugarch.spec)),
                          dccOrder=c(1,1), distribution="mvnorm")
  } else {
    stop("Model does not exist")
  }

  message("Estimating model")
  if (model=="TVP-VAR") {
    fit = var_model(x, configuration=configuration)
    B_t = fit$B_t
    Q_t = fit$Q_t
  } else if (model=="DCC-GARCH") {
    fit = rmgarch::dccfit(mgarch.spec, data=x)
    fevd = VFEVD(fit, nfore=nfore, standardize=DCC_config$standardize)$FEVD
    Q_t = fevd
  } else {
    Q_t = array(NA, c(k, k, t0))
    B_t = array(NA, c(k, k*nlag, t0))
    pb = progress_bar$new(total=t0)
    for (i in 1:t0) {
      fit = var_model(x[i:(i+window.size-1),], configuration=configuration)
      B_t[,,i] = fit$B
      Q_t[,,i] = fit$Q
      pb$tick()
    }
  }
  DATE = as.character(zoo::index(x))
  date = DATE[(length(DATE)-dim(Q_t)[3]+1):length(DATE)]
  dimnames(Q_t)[[1]] = dimnames(Q_t)[[2]] = NAMES
  dimnames(Q_t)[[3]] = as.character(date)
  
  message("Computing connectedness measures")
  if (connectedness=="Time") {
    generalized = Connectedness_config$TimeConnectedness$generalized
    if (model=="DCC-GARCH") {
      dca = TimeConnectedness(FEVD=fevd, corrected=corrected)
      message("The DCC-GARCH connectedness approach is implemented according to:\n Gabauer, D. (2020). Volatility impulse response analysis for DCC-GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting, 39(5), 788-796.")
    } else {
      dca = TimeConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore,
                              generalized=generalized,
                              corrected=corrected)
      if (model=="VAR" && !generalized) {
        message("The (orthogonalized) VAR connectedness approach is implemented according to:\n Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal, 119(534), 158-171.")
      } else if (model=="VAR" && generalized) {
        message("The (generalized) VAR connectedness approach is implemented according to:\n Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting, 28(1), 57-66.")
      } else if (model=="TVP-VAR") {
        message("The TVP-VAR connectedness approach is implemented according to:\n Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management, 13(4), 84.")
      } else if (model=="QVAR") {
        message("The QVAR connectedness approach is implemented according to:\n Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters, 204, 109891.")
      } else if (model=="LASSO" || model=="Ridge" || model=="Elastic") {
        message("The Elastic Net and its restricted models, namely, the LASSO and Ridge VAR connectedness approach are implemented according to:\n Gabauer, D., Gupta, R., Marfatia, H., & Miller, S. (2020). Estimating US Housing Price Network Connectedness: Evidence from Dynamic Elastic Net, Lasso, and Ridge Vector Autoregressive Models (No. 202065). University of Pretoria, Department of Economics.")
      }
    }
  } else if (connectedness=="Frequency") {
    dca = FrequencyConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore,
                                 partition=Connectedness_config$FrequencyConnectedness$partition,
                                 generalized=Connectedness_config$FrequencyConnectedness$generalized,
                                 scenario=Connectedness_config$FrequencyConnectedness$scenario, 
                                 corrected=corrected)
    if (model=="VAR") {
      message("The VAR frequency connectedness approach is implemented according to:\n Barunik, J., & Krehlik, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics, 16(2), 271-296.")
    } else if (model=="TVP-VAR") {
      message("The TVP-VAR frequency connectedness approach is implemented according to:\n Chatziantoniou, I., Gabauer, D., & Gupta, R. (2021). Integration and Risk Transmission in the Market for Crude Oil: A Time-Varying Parameter Frequency Connectedness Approach (No. 202147).")
    } else if (model=="QVAR") {
      message("The QVAR frequency connectedness approach is implemented according to:\n Chatziantoniou, I., Aikins Abakah, E. J., Gabauer, D., & Tiwari, A. K. (2021). Quantile time-frequency price connectedness between green bond, green equity, sustainable investments and clean energy markets: Implications for eco-friendly investors. Available at SSRN 3970746.")
    }
  } else if (connectedness=="Joint") {
    dca = JointConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore)
    if (model=="VAR") {
      message("The VAR joint connectedness approach is implemented according to:\n Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.")
    }
  } else if (connectedness=="Extended Joint") {
    dca = ExtendedJointConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore)
    if (model=="VAR" || model=="TVP-VAR") {
      message("The VAR and TVP-VAR extended joint connectedness approach is implemented according to:\n Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy, 73, 102219.")
    } else if (model=="QVAR") {
      message("The QVAR extended joint connectedness approach is implemented according to:\n Cunado, J, Chatziantoniou, I., Gabauer, D., Hardik, M., & de Garcia, F.P. (2022). Dynamic spillovers across precious metals and energy realized volatilities: Evidence from quantile extended joint connectedness measures.")
    }
  } else {
    stop("Connectedness approach does not exist")
  }
  dca
}

