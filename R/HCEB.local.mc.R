#' Permutation test for Local Heteroscedastic Consistent Empirical Bayes Moran's I
#'
#' A permutation test for Local HC-EB Moran's I statistics by using nsim random permutations of HC-EB rate.
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

HCEB.local.mc <- function(x, n, k, b=0, listw, sig = 0.05, nsim, zero.policy = NULL)
{
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours) # Number of areas
  N = sum(b+n)  # sum of the population-at-risk
  m = sum(x)/N  # global average rate

  if (n1 != length(x) | n1!=length(n))
    stop("objects of different length")

  I.local <- HCEB.local(x, n, k, b, listw, zero.policy)
  z <- HCEB.est(x, n, k, b)$HCEB-m
  lz <- lag.listw(listw, z, zero.policy = zero.policy)

  I.perm <- matrix(data=NA, nrow=nsim, ncol=length(x))
  for(i in 1:nsim){
    perm <- sample(n1, replace=F)
    I.local.perm <- HCEB.local(x[perm], n[perm], k[perm], b, listw, zero.policy)
    I.perm[i,] <- I.local.perm
    rm(I.local.perm, perm)
  }

  res <- data.frame(I.local=I.local, p.perm=NA, quad_sig=NA)

  res$p.perm <- rowSums(t(I.perm)>I.local)/nsim
  res$quad_sig <- NA

  res[(z >= 0 & lz >= 0) & (res$p.perm <= sig), "quad_sig"] <- "High-High"
  res[(z <= 0 & lz <= 0) & (res$p.perm <= sig), "quad_sig"] <- "Low-Low"
  res[(z >= 0 & lz <= 0) & (res$p.perm <= sig), "quad_sig"] <- "High-Low"
  res[(z <= 0 & lz >= 0) & (res$p.perm <= sig), "quad_sig"] <- "Low-High"

  res
}
environment(HCEB.local.mc) <- asNamespace("spdep")
