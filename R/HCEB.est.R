#' Heteroscedastic Consistent Empirical Bayes rate
#'
#' Compute HC-EB rate estimates for rates smoothed to the overall mean by their denominator data uncertainty.
#' @param x a numeric vector of counts of cases
#' @param n a numeric vector of populations at risk
#' @param k a numeric vector of margin of errors of x
#' @param b default 0, base populations-at-risk
#' @references #####, ##### and #####. 2018. Small area adjustment of spatial autocorrelation estimators for heteroscedastic sampling errors, Geographical Analysis (under review)
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.est <- function(x, n, k, b=0)
{

  if (!is.numeric(x))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  if (!is.numeric(n))
    stop(paste(deparse(substitute(n)), "is not a numeric vector"))
  if (!is.numeric(k))
    stop(paste(deparse(substitute(k)), "is not a numeric vector"))
  if (any(is.na(x)))
    stop("NA in at risk population")
  if (any(is.na(n)))
    stop("NA in cases")
  if (any(is.na(k)))
    stop("NA in margin of errors of risk population")
  if (any(n+b < 0))
    stop("negative number of cases")
  if (any(n+b == 0))
    stop("zero value in populations at risk")
  if (length(x) != length(n) | length(n) != length(k) | length(x) != length(k))
    stop("vectors of different length")

  r <- x/(b+n) # raw rate
  N = sum(b+n) # sum of the population-at-risk
  m = sum(x)/N # global average rate
  CV = (k/1.645)/(b+n)
  s2 = sum((b+n)*(r-m)^2)/N
  A = s2-m/(N/233)
  B = m/(b+n)
  C = r^2*CV^2
  w_EB = A/(A+B)
  w_HCEB = A/(A+B+C)

  rate_EB = w_EB*(r-m)+m
  rate_HCEB = w_HCEB*(r-m)+m

  res <- data.frame(raw=r, EB=rate_EB, HCEB=rate_HCEB)
  attr(res, "CV") <- CV
  attr(res, "EB_weight") <- w_EB
  attr(res, "HCEB_weight") <- w_HCEB

  res
}
environment(HCEB.est) <- asNamespace("spdep")
