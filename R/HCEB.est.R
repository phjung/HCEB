#' Heteroscedastic Consistent Empirical Bayes rate
#'
#' Compute HC-EB rate estimates for rates smoothed to the overall mean by their denominator data uncertainty.
#' @param n a numeric vector of counts of cases (numerator)
#' @param x a numeric vector of populations at risk (denominator)
#' @param se a numeric vector of sampling standard errors of corresponding populations at risk
#' @return
#' @references Jung, PH, Thill J-C, Issel M 2018 Spatial Autocorrelation Statistics of Areal Prevalence Rates under High Uncertainty in Denominator Data, Geographical Analysis
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.est <- function(n, x, se)
{
  if (!is.numeric(n))
    stop(paste(deparse(substitute(n)), "is not a numeric vector"))
  if (!is.numeric(x))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  if (!is.numeric(se))
    stop(paste(deparse(substitute(se)), "is not a numeric vector"))
  if (any(is.na(n)))
    stop("NA in cases")
  if (any(is.na(x)))
    stop("NA in at risk population")
  if (any(is.na(se)))
    stop("NA in stadard errors of risk population")
  if (any(n < 0))
    stop("negative number of cases")
  if (any(x < 0))
    stop("negative number of risk population")
  if (length(x) != length(n) | length(n) != length(se) | length(x) != length(se))
    stop("vectors of different length")

  N <- sum(n)
  X <- sum(x)
  m <- N/X # global average rate

  r <- n/x # raw rate
  r <- ifelse(is.infinite(r), m, r) # +/0 => global average rate
  r <- ifelse(is.nan(r), 0, r) # 0/0 => 0
  r_std <- r-m

  cv <- se/x
  s2 <- sum(x*(r_std)^2/x, na.rm=T)

  A <- s2-m/(N/length(x))
  B <- m/x
  C <- cv^2*r*(r-2*m)

  w_EB <- A/(A+B)
  w_EB <- ifelse(w_EB==0, 0, w_EB)

  w_HCEB <- A/(A+B+C)
  w_HCEB <- ifelse(x==0, 0, w_HCEB)

  r_EB <- w_EB*(r_std)+m
  r_HCEB <- w_HCEB*(r_std)+m

  res <- data.frame(raw=r, EB=r_EB, HCEB=r_HCEB)
  attr(res, "weight_EB") <- w_EB
  attr(res, "weight_HCEB") <- w_HCEB
  attr(res, "m") <- m
  res
}
