#' Local Heteroscedastic Consistent Empirical Bayes Moran's I statistic
#'
#' Compute Local HC-EB Moran's I. Local Moran's I with each case weighted by HC-EB method
#' @param x a numeric vector of counts of cases
#' @param n a numeric vector of populations at risk
#' @param k a numeric vector of margin of errors of x
#' @param b default 0, base populations-at-risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @references #####, ##### and #####. 2018. Small area adjustment of spatial autocorrelation estimators for heteroscedastic sampling errors, Geographical Analysis (under review)
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.local <- function(x, n, k, b=0, listw, zero.policy = NULL)
{
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours) # Number of areas
  if (n1 != length(x) | n1!=length(n))
    stop("objects of different length")


  r <- x/(b+n)  # raw rate
  N = sum(b+n)  # sum of the population-at-risk
  m = sum(x)/N  # global average rate

  z <- HCEB.est(x, n, k, b)$HCEB-m
  lz <- lag.listw(listw, z, zero.policy = zero.policy)
  s2 <- sum(z^2)/n1

  res <- (z/s2) * lz

  res

}
environment(HCEB.local) <- asNamespace("spdep")
