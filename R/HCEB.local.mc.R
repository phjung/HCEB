#' Permutation test for Local Heteroscedasticity-Consistent Empirical Bayes Moran's I
#'
#' A permutation test for Local HC-EB Moran's I statistics by using nsim random permutations of HC-EB rate.
#' @param n a numeric vector of counts of cases (numerator)
#' @param x a numeric vector of populations at risk (denominator)
#' @param se a numeric vector of sampling standard errors of corresponding populations at risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param sig default 0.05, pseudo significance level for permutation test
#' @param nsim number of permutations
#' @references #####, ##### and #####. 2018. Small area adjustment of spatial autocorrelation estimators for heteroscedastic sampling errors, Geographical Analysis (under review)
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.local.mc <- function(n, x, se, listw, zero.policy = NULL, sig = 0.05, nsim = 1000)
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

  I.local <- HCEB.local(n, x, se, listw, zero.policy)

  I.perm <- matrix(data=NA, nrow=nsim, ncol=length(x))
  for(i in 1:nsim){
    perm <- sample(n1, replace=F)
    I.local.perm <- HCEB.local(n[perm], x[perm], se[perm], listw, zero.policy)
    I.perm[i,] <- I.local.perm
    rm(I.local.perm, perm)
  }

  res <- data.frame(I.local.HCEB=I.local, p.perm=NA, quad_sig=NA)

  res$p.perm <- ifelse(I.local<0,
                       rowSums(t(I.perm)<I.local)/nsim,
                       rowSums(t(I.perm)>I.local)/nsim)

  res$quad_sig <- NA

  res[(z >= 0 & lz >= 0) & (res$p.perm <= sig), "quad_sig"] <- "High-High"
  res[(z <= 0 & lz <= 0) & (res$p.perm <= sig), "quad_sig"] <- "Low-Low"
  res[(z >= 0 & lz <= 0) & (res$p.perm <= sig), "quad_sig"] <- "High-Low"
  res[(z <= 0 & lz >= 0) & (res$p.perm <= sig), "quad_sig"] <- "Low-High"

  res
}
environment(HCEB.local.mc) <- asNamespace("spdep")
