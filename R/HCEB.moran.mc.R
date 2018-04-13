#' Permutation test for Global Heteroscedastic Consistent Empirical Bayes Moran's I
#'
#' A permutation test for Global HC-EB Moran's I statistics by using nsim random permutations of HC-EB rate.
#' @param x a numeric vector of counts of cases
#' @param n a numeric vector of populations at risk
#' @param k a numeric vector of margin of errors of x
#' @param b default 0, base populations-at-risk
#' @param listw a listw object created for example by nb2listw
#' @param sig default 0.05, pseudo significance level for permutation test
#' @param nsim number of permutations
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @references #####, ##### and #####. 2018. Small area adjustment of spatial autocorrelation estimators for heteroscedastic sampling errors, Geographical Analysis (under review)
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.moran.mc <- function(x, n, k, b=0, listw, sig = 0.05, nsim, zero.policy = NULL)
{
  I.global <- HCEB.moran(x, n, k, b, listw, zero.policy)
  n1 <- length(listw$neighbours)
  I.perm <- numeric(length = nsim)

  for (i in 1:nsim) {
    perm <- sample(n1, replace=F)
    I.global.perm <- HCEB.moran(x[perm], n[perm], k[perm], b, listw, zero.policy)
    I.perm[i] <- I.global.perm
    rm(I.global.perm, perm)
  }
  p.perm <- sum(t(I.perm)>I.global)/nsim

  res <- list("I"=I.global, "p.perm"=p.perm)

  res
}
environment(HCEB.moran.mc) <- asNamespace("spdep")
