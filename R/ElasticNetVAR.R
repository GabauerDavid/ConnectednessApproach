
#' @title Elastic Net vector autoregression
#' @description Estimation of a VAR using equation-by-equation LASSO, Ridge or Elastic Net regressions.
#' @param x zoo data matrix
#' @param nlag Lag length
#' @param nfolds N-fold cross validation
#' @param loss Loss function
#' @param alpha LASSO is alpha equal 1 and Ridge if alpha equal 0
#' @param n_alpha Creates n-equidistant alpha values
#' @param configuration Model configuration
#' @return Estimate VAR model
#' @examples
#' \donttest{
#' data("dy2012")
#' fit = ElasticNetVAR(dy2012, configuration=list(nlag=1, alpha=1, nfolds=10, loss="mae"))
#' }
#' @import glmnet
#' @importFrom stats predict
#' @references
#' Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor, J., & Tibshirani, R. J. (2012). Strong rules for discarding predictors in lasso‐type problems. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 74(2), 245-266.
#' 
#' Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1), 55-67.
#' 
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
#' 
#' Demirer, M., Diebold, F. X., Liu, L., & Yilmaz, K. (2018). Estimating global bank network connectedness. Journal of Applied Econometrics, 33(1), 1-15.
#' 
#' Gabauer, D., Gupta, R., Marfatia, H., & Miller, S. M. (2020). Estimating US Housing Price Network Connectedness: Evidence from Dynamic Elastic Net, Lasso, and Ridge Vector Autoregressive Models. Lasso, and Ridge Vector Autoregressive Models (July 26, 2020).
#' @author David Gabauer
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet predict.glmnet
#' @export
ElasticNetVAR = function(x, configuration=list(nlag=1, nfolds=10, loss="mae", alpha=NULL, n_alpha=10)) {
  nlag = configuration$nlag
  alpha_ = configuration$alpha
  nfolds = configuration$nfolds
  n_alpha = configuration$n_alpha
  loss = configuration$loss
  if (nlag<=0) {
    stop("nlag needs to be a positive integer")
  }
  if (nfolds<=0) {
    stop("nfolds needs to be a positive integer")
  }
  if (!is(x, "zoo")) {
    stop("Data needs to be of type 'zoo'")
  }
  NAMES = colnames(x)
  k = ncol(x)
  if (is.null(NAMES)) {
    NAMES = 1:k
  }
  if (is.null(alpha_) ) {
    alpha_ = seq(0, 1, 1/n_alpha)
  }
  alpha = Res = B = NULL
  for (i in 1:k) {
    MAE = NULL
    B_ = array(NA, c(nlag*k, length(alpha_)))
    for (j in 1:length(alpha_)) {
      z = embed(x, nlag+1)
      X = z[,-c(1:k)]
      y = z[,i]
      fit = glmnet::cv.glmnet(X, y, alpha=alpha_[j], type.measure=loss, nfolds=nfolds)
      B_[,j] = coef(fit)[-1]
      y_pred = cbind(1, X) %*% coef(fit)
      MAE[j] = mean(abs(y - y_pred[,1]))
      if (MAE[j] <= min(MAE)) {
        Y_pred = y_pred[,1]
      }
    }
    ind = which(MAE==min(MAE))[1]
    B = rbind(B, B_[,ind])
    Res = cbind(Res, y-Y_pred)
    alpha[i] = alpha_[ind]
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(as.character(zoo::index(x)),1)))
  results = list(B=B, Q=Q, alpha=alpha)
}
