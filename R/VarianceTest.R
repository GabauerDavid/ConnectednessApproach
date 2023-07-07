
#' @title Variance Test
#' @description VarianceTest performs variance homogeneity tests including Ftest, Bartlett, Brown-Forsythe and Fligner-Killeen tests.
#' @param formula a formula of the form lhs ~ rhs where lhs gives the sample values and rhs the corresponding groups.
#' @param data a tibble or data frame containing the variables in the formula formula
#' @param alpha the level of significance to assess variance homogeneity. Default is set to alpha = 0.05.
#' @param method a character string to select one of the variance homogeneity tests: "Bartlett", "Brown-Forsythe", "Fisher" and "Fligner-Killeen".
#' @param na.rm Ha logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return Get bivariate portfolio weights
#' @importFrom stats var.test
#' @importFrom stats fligner.test
#' @importFrom stats bartlett.test
#' @importFrom stats complete.cases
#' @importFrom stats model.frame
#' @importFrom stats pf
#' @importFrom stats na.omit
#' @importFrom car leveneTest
#' @references
#' Antonakakis, N., Cunado, J., Filis, G., Gabauer, D., & de Gracia, F. P. (2020). Oil and asset classes implied volatilities: Investment strategies and hedging effectiveness. Energy Economics, 91, 104762.
#' @author David Gabauer
#' @export
VarianceTest = function (formula, data, alpha=0.05, method=c('Bartlett', 'Brown-Forsythe', 'Fligner-Killeen', 'Fisher', 'Levene'), na.rm=TRUE) {
  data = model.frame(formula, data)
  dp = as.character(formula)
  DNAME = paste(dp[[2L]], "and", dp[[3L]])
  if (any(colnames(data) == dp[[3L]]) == FALSE) {
    stop("The name of group variable does not match the variable names in the data. The group variable must be one factor.")
  }
  if (any(colnames(data) == dp[[2L]]) == FALSE) {
    stop("The name of response variable does not match the variable names in the data.")
  }
  y = data[[dp[[2L]]]]
  group = data[[dp[[3L]]]]
  
  if (!(is.factor(group) | is.character(group))) {
    stop("The group variable must be a factor or a character.")
  }
  if (is.character(group)) {
    group = as.factor(group)
  }
  if (!is.numeric(y)) {
    stop("The response must be a numeric variable.")
  }
  if (na.rm) {
    completeObs = complete.cases(y, group)
    y = y[completeObs]
    group = group[completeObs]
  }
  
  if (method=="Brown-Forsythe") {
    n = length(y)
    x.levels = levels(factor(group))
    y.vars = y.means = m = y.n = NULL
    y.mean = mean(y)
    for (i in x.levels) {
      y.vars[i] = var(y[group == i])
      y.means[i] = mean(y[group == i])
      y.n[i] = length(y[group == i])
    }
    for (j in x.levels) {
      m[j] = (1 - y.n[j]/n) * (y.vars[j])/sum((1 - y.n/n) *  (y.vars))
    }
    SSb = sum(y.n * ((y.means - y.mean)^2))
    denom = sum((1 - y.n/n) * (y.vars))
    df1 = length(x.levels) - 1
    df2 = 1/(sum(m^2/(y.n - 1)))
    statistic = as.numeric(SSb/denom)
    p.value = pf(statistic, df1, df2, lower.tail=F)
  } else if (method == "Bartlett") {
    out = stats::bartlett.test(y, group)
    statistic = as.numeric(out$statistic)
    p.value = out$p.value
  } else if (method == "Fligner-Killeen") {
    out = stats::fligner.test(y, group)
    statistic = as.numeric(out$statistic)
    p.value = out$p.value
  } else if (method == "Fisher") {
    groups = as.character(unique(group))
    out = stats::var.test(y[which(group==groups[1])], y[which(group!=groups[1])])
    statistic = out$statistic
    p.value = out$p.value
  } else if (method == "Levene") {
    out = leveneTest(y, group)
    statistic = na.omit(out$`F value`)[1]
    p.value = na.omit(out$`Pr(>F)`)[1]
  } else {
    stop("This variance test does not exist.")
  }
  result = list(statistic=statistic, p.value=p.value)
}
