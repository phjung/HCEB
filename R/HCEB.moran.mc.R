#' Permutation test for Global Heteroscedastic Consistent Empirical Bayes Moran's I
#'
#' A permutation test for Global HC-EB Moran's I statistics by using nsim random permutations of HC-EB rates.
#' @param n a numeric vector of counts of cases (numerator)
#' @param x a numeric vector of populations at risk (denominator)
#' @param se a numeric vector of sampling standard errors of corresponding populations at risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param sig default 0.05, pseudo significance level for permutation test
#' @param nsim number of permutations
#' @references Jung, PH, Thill J-C, Issel M 2018 Spatial Autocorrelation Statistics of Areal Prevalence Rates under High Uncertainty in Denominator Data, Geographical Analysis
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.moran.mc <- function(n, x, se, zero.policy = NULL, listw, sig = 0.05, nsim = 1000)
{
  I.global <- HCEB.moran(n, x, se, listw, zero.policy)
  n1 <- length(listw$neighbours)
  I.perm <- numeric(length = nsim)

  for (i in 1:nsim) {
    perm <- sample(n1, replace=F)
    I.global.perm <- HCEB.moran(x[perm], n[perm], k[perm], b, listw, zero.policy)
    I.perm[i] <- I.global.perm
    rm(I.global.perm, perm)
  }
  p.perm <- sum(t(I.perm)>I.global)/nsim

  res <- list("I_HCEB"=I.global, "p_perm"=p.perm)

  res
}
