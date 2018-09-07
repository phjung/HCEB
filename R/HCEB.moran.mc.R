#' Permutation test for Global Heteroscedasticity-Consistent Empirical Bayes Moran's I
#'
#' A permutation test for Global HC-EB Moran's I statistics by using nsim random permutations of HC-EB rates.
#' @param n a numeric vector of counts of cases (numerator)
#' @param x a numeric vector of populations at risk (denominator)
#' @param se a numeric vector of sampling standard errors of corresponding populations at risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param sig default 0.05, pseudo significance level for permutation test
#' @param nsim number of permutations
#' @return
#' \item{I.HCEB}{Heteroscedasticity-consistent empirical Bayes Global Moran's I}
#' \item{p.perm}{pseudo p-value of Monte-Carlo permutation test}
#' @references Jung, PH, Thill J-C, Issel M 2018 Spatial Autocorrelation Statistics of Areal Prevalence Rates under High Uncertainty in Denominator Data, Geographical Analysis
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.moran.mc <- function(n, x, se, listw, zero.policy = NULL, sig = 0.05, nsim = 1000)
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
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))

  n1 <- length(listw$neighbours)
  if (n1 != length(x))
    stop("objects of different length")

  I.global <- HCEB.moran(n, x, se, listw, zero.policy)
  I.perm <- numeric(length = nsim)

  for (i in 1:nsim) {
    perm <- sample(n1, replace=F)
    I.global.perm <- HCEB.moran(n[perm], x[perm], se[perm], listw, zero.policy)
    I.perm[i] <- I.global.perm
    rm(I.global.perm, perm)
  }
  p.perm <- ifelse(I.global<0,
                   sum(t(I.perm)<I.global)/nsim,
                   sum(t(I.perm)>I.global)/nsim)

  res <- list("I_HCEB"=I.HCEB, "p_perm"=p.perm)
  res
}
environment(HCEB.moran.mc) <- asNamespace("spdep")
