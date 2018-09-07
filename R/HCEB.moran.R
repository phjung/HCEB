#' Global Heteroscedastic Consistent Empirical Bayes Moran's I estimator
#'
#' Compute Global HC-EB Moran's I. Global Moran's I of areal prevalence rates with uncertainty in denominator.
#' @param n a numeric vector of counts of cases (numerator)
#' @param x a numeric vector of populations at risk (denominator)
#' @param se a numeric vector of sampling standard errors of corresponding populations at risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @references Jung, PH, Thill J-C, Issel M 2018 Spatial Autocorrelation Statistics of Areal Prevalence Rates under High Uncertainty in Denominator Data, Geographical Analysis
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.moran <- function (n, x, se, listw, zero.policy = NULL)
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
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy")
  stopifnot(is.logical(zero.policy))

  n1 <- length(listw$neighbours)
  if (n1 != length(x))
    stop("objects of different length")

  d <- HCEB.est(n, x, se)
  m <- attributes(d)$m

  z <- d$HCEB-m
  lz <- lag.listw(listw, z, zero.policy = zero.policy)
  zz <- sum(z^2)

  I <- (sum(z * lz))/zz
  res <- I

  res
}
