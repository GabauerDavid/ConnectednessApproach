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
#'                             model="VAR",
#'                             connectedness="Time",
#'                             VAR_config=list(TVPVAR=list(kappa1=0.99, kappa2=0.96,
#'                                             prior="MinnesotaPrior", gamma=0.1)))
#' dca$TABLE
#' }
#' @import progress
#' @importFrom stats cor
#' @importFrom rugarch ugarchspec
#' @importFrom rugarch multispec
#' @importFrom rmgarch dccspec
#' @importFrom rmgarch dccfit
#' @references
#' Adekoya, O. B., Akinseye, A. B., Antonakakis, N., Chatziantoniou, I., Gabauer, D., & Oliyide, J. (2022). Crude oil and Islamic sectoral stocks: Asymmetric TVP-VAR connectedness and investment strategies. Resources Policy.
#' 
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management.
#' 
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics.
#' 
#' Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2021). The impact of Euro through time: Exchange rate dynamics under different regimes. International Journal of Finance & Economics.
#' 
#' Balcilar, M., Gabauer, D., & Umar, Z. (2021). Crude Oil futures contracts and commodity markets: New evidence from a TVP-VAR extended joint connectedness approach. Resources Policy.
#' 
#' Balli, F., Balli, H. O., Dang, T. H. N., & Gabauer, D. (2023). Contemporaneous and lagged R2 decomposed connectedness approach: New evidence from the energy futures market. Finance Research Letters.
#' 
#' Barunik, J., & Krehlik, T. (2018). Measuring the frequency dynamics of financial connectedness and systemic risk. Journal of Financial Econometrics.
#' 
#' Broadstock, D. C., Chatziantoniou, I., & Gabauer, D. (2022). Minimum connectedness portfolios and the market for green bonds: Advocating socially responsible investment (SRI) activity. In Applications in energy finance: The energy sector, economic activity, financial markets and the environment. Cham: Springer International Publishing.
#' 
#' Chatziantoniou, I., & Gabauer, D. (2021). EMU risk-synchronisation and financial fragility through the prism of dynamic connectedness. The Quarterly Review of Economics and Finance.
#' 
#' Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters.
#' 
#' Chatziantoniou, I., Gabauer, D., & Gupta, R. (2023). Integration and risk transmission in the market for crude oil: New evidence from a time-varying parameter frequency connectedness approach. Resources Policy.
#' 
#' Chatziantoniou, I., Aikins Abakah, E. J., Gabauer, D., & Tiwari, A. K. (2022). Quantile time-frequency price connectedness between green bond, green equity, sustainable investments and clean energy markets. Journal of Cleaner Production.
#' 
#' Chatziantoniou, I., Elsayed, A. H., Gabauer, D., & Gozgor, G. (2023). Oil price shocks and exchange rate dynamics: Evidence from decomposed and partial connectedness measures for oil importing and exporting economies. Energy Economics.
#' 
#' Cocca, T., Gabauer, D., & Pomberger, S. (2024). Clean energy market connectedness and investment strategies: New evidence from DCC-GARCH R2 decomposed connectedness measures. Energy Economics.
#' 
#' Cunado, J., Chatziantoniou, I., Gabauer, D., de Gracia, F. P., & Hardik, M. (2023). Dynamic spillovers across precious metals and oil realized volatilities: Evidence from quantile extended joint connectedness measures. Journal of Commodity Markets.
#' 
#' Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal.
#' 
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting.
#' 
#' Gabauer, D. (2020). Volatility impulse response analysis for DCC-GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting.
#' 
#' Gabauer, D. (2021). Dynamic measures of asymmetric & pairwise connectedness within an optimal currency area: Evidence from the ERM I system. Journal of Multinational Financial Management.
#' 
#' Gabauer, D., Chatziantoniou, I., & Stenfors, A. (2023). Model-free connectedness measures. Finance Research Letters.
#' 
#' Gabauer, D., Gupta, R., Marfatia, H. A., & Miller, S. M. (2024). Estimating US housing price network connectedness: Evidence from dynamic Elastic Net, Lasso, and ridge vector autoregressive models. International Review of Economics & Finance.
#' 
#' Gabauer, D., & Stenfors, A. (2024). Quantile-on-quantile connectedness measures: Evidence from the US treasury yield curve. Finance Research Letters, 60, 104852.
#' 
#' Lastrapes, W. D., & Wiesen, T. F. (2021). The joint spillover index. Economic Modelling, 94, 681-691.
#' 
#' Naeem, M. A., Chatziantoniou, I., Gabauer, D., & Karim, S. (2024). Measuring the G20 stock market return transmission mechanism: Evidence from the R2 connectedness approach. International Review of Financial Analysis.
#' 
#' Stenfors, A., Chatziantoniou, I., & Gabauer, D. (2022). Independent policy, dependent outcomes: A game of cross-country dominoes across European yield curves. Journal of International Financial Markets, Institutions and Money.
#' 
#' Zhang, Y., Gabauer, D., Gupta, R., & Ji, Q. (2024). How connected is the oil-bank network? Firm-level and high-frequency evidence. Energy Economics.
#' 
#' @author David Gabauer
#' @export
ConnectednessApproach = function(x,
                                 nlag=1, 
                                 nfore=10, 
                                 window.size=NULL, 
                                 corrected=FALSE,
                                 model=c("VAR", "QVAR", "LAD", "LASSO", "Ridge", "Elastic", "TVP-VAR", "DCC-GARCH"),
                                 connectedness=c("Time","Frequency", "Joint", "Extended Joint", "R2"),
                                 VAR_config=list(
                                   VAR=list(method="pearson"),
                                   QVAR=list(tau=0.5, method="fn"),
                                   ElasticNet=list(nfolds=10, alpha=NULL, loss="mae", n_alpha=10),
                                   TVPVAR=list(kappa1=0.99, kappa2=0.99, prior="BayesPrior", gamma=0.01)),
                                 DCC_config=list(standardize=FALSE),
                                 Connectedness_config = list(
                                   TimeConnectedness=list(generalized=TRUE),
                                   FrequencyConnectedness=list(partition=c(pi,pi/2,0), generalized=TRUE, scenario="ABS"),
                                   R2Connectedness=list(method="pearson", decomposition=TRUE, relative=FALSE, tau=NULL)
                                 )) {
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  model = match.arg(model)
  if (length(connectedness)>1) {
    connectedness = "Time"
  } else {
    connectedness = match.arg(connectedness)
  }

  if (nlag<=0 & connectedness!="R2") {
    stop("nlag needs to be a positive integer")
  }
  
  NAMES = colnames(x)
  k = ncol(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  if (connectedness=="R2") {
    if (!Connectedness_config$R2Connectedness$decomposition) {
      nlag = 0
    }
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
    configuration = list(nlag=nlag, method=VAR_config$VAR$method)
  } else if (model=="QVAR") {
    var_model = QVAR
    if (is.null(VAR_config$QVAR$method)) {
      VAR_config$QVAR$method = "fn"
    }
    configuration = list(nlag=nlag, tau=VAR_config$QVAR$tau, method=VAR_config$QVAR$method)
  } else if (model=="LAD") {
    var_model = LADVAR
    configuration = list(nlag=nlag)
  } else if (model=="LASSO") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=1, nfolds=VAR_config$ElasticNet$nfolds, loss=VAR_config$ElasticNet$loss)
  } else if (model=="Ridge") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=0, nfolds=VAR_config$ElasticNet$nfolds, loss=VAR_config$ElasticNet$loss)
  } else if (model=="Elastic") {
    var_model = ElasticNetVAR
    configuration = list(nlag=nlag, alpha=VAR_config$ElasticNet$alpha, nfolds=VAR_config$ElasticNet$nfolds,
                         loss=VAR_config$ElasticNet$loss, n_alpha=VAR_config$ElasticNet$n_alpha)
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
    if (connectedness!="R2") {
      pb = progress_bar$new(total=t0)
      for (i in 1:t0) {
        fit = var_model(x[i:(i+window.size-1),], configuration=configuration)
        B_t[,,i] = fit$B
        Q_t[,,i] = fit$Q
        pb$tick()
      }
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
      message("The DCC-GARCH connectedness approach is implemented according to:\n Gabauer, D. (2020). Volatility impulse response analysis for DCC-GARCH models: The role of volatility transmission mechanisms. Journal of Forecasting.")
    } else {
      dca = TimeConnectedness(Phi=B_t, Sigma=Q_t, nfore=nfore,
                              generalized=generalized,
                              corrected=corrected)
      if (model=="VAR" && !generalized) {
        message("The (orthogonalized) VAR connectedness approach is implemented according to:\n Diebold, F. X., & Yilmaz, K. (2009). Measuring financial asset return and volatility spillovers, with application to global equity markets. The Economic Journal.")
      } else if (model=="VAR" && generalized) {
        message("The (generalized) VAR connectedness approach is implemented according to:\n Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive directional measurement of volatility spillovers. International Journal of Forecasting.")
      } else if (model=="TVP-VAR") {
        message("The TVP-VAR connectedness approach is implemented according to:\n Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020). Refined measures of dynamic connectedness based on time-varying parameter vector autoregressions. Journal of Risk and Financial Management.")
      } else if (model=="QVAR") {
        message("The QVAR connectedness approach is implemented according to:\n Chatziantoniou, I., Gabauer, D., & Stenfors, A. (2021). Interest rate swaps and the transmission mechanism of monetary policy: A quantile connectedness approach. Economics Letters.")
      } else if (model=="LASSO" || model=="Ridge" || model=="Elastic") {
        message("The Elastic Net and its restricted models, namely, the LASSO and Ridge VAR connectedness approach are implemented according to:\n Gabauer, D., Gupta, R., Marfatia, H. A., & Miller, S. M. (2024). Estimating US housing price network connectedness: Evidence from dynamic Elastic Net, Lasso, and ridge vector autoregressive models. International Review of Economics & Finance.")
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
  } else if (connectedness=="R2") {
    if (Connectedness_config$R2Connectedness$decomposition) {
      if (is.null(Connectedness_config$R2Connectedness$tau)) {
        if (nlag>0) {
          message("The contemporaneous R2 connectedness approach is implemented according to:\n Naeem, M. A., Chatziantoniou, I., Gabauer, D., & Karim, S. (2023). Measuring the G20 Stock Market Return Transmission Mechanism: Evidence From the R2 Connectedness Approach. International Review of Financial Analysis.\n")
          message("The generalized R2 connectedness approach is implemented according to:\n Balli, F., Balli, H. O., Dang, T. H. N., & Gabauer, D. (2023). Contemporaneous and lagged R2 decomposed connectedness approach: New evidence from the energy futures market. Finance Research Letters, 57, 104168.")
        } else {
          message("The contemporaneous R2 connectedness approach is implemented according to:\n Naeem, M. A., Chatziantoniou, I., Gabauer, D., & Karim, S. (2023). Measuring the G20 Stock Market Return Transmission Mechanism: Evidence From the R2 Connectedness Approach. International Review of Financial Analysis.")
        }
        dca = R2Connectedness(x, nlag=nlag, window.size=window.size, method=Connectedness_config$R2Connectedness$method,
                              relative=Connectedness_config$R2Connectedness$relative, tau=NULL)
      } else {
        dca = R2Connectedness(x, nlag=nlag, window.size=window.size, method=Connectedness_config$R2Connectedness$method,
                              relative=Connectedness_config$R2Connectedness$relative, tau=Connectedness_config$R2Connectedness$tau)
      }
    } else {
      fevd = Q_t
      for (i in 1:t0) {
        ct = cor(x[i:(i+window.size-1),], method=Connectedness_config$R2Connectedness$method)^2
        if (Connectedness_config$R2Connectedness$method=="kendall") {
          ct = sin(0.5*pi*ct)
        }
        fevd[,,i] = ct/rowSums(ct)
      }
      dca = TimeConnectedness(FEVD=fevd, corrected=corrected)
      message("The unconditional connectedness approach is implemented according to:\n Gabauer, D, Chatziantoniou, I., & Stenfors, A. (2023). Model-free connectedness measures.")
    }
  } else {
    stop("Connectedness approach does not exist")
  }
  dca
}
