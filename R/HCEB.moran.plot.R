#' Moran scatterplot with Heteroscedasticity-Consistent Empirical Bayes rate
#'
#' A plot of HC-EB spatial areal rate data against its spatially lagged values.
#' @param n a numeric vector of counts of cases (numerator)
#' @param x a numeric vector of populations at risk (denominator)
#' @param se a numeric vector of sampling standard errors of corresponding populations at risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param spChk should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use get.spChkOption()
#' @param labels character labels for points with high influence measures, if set to FALSE, no labels are plotted for points with large influence
#' @references Jung, PH, Thill J-C, Issel M 2018 Spatial Autocorrelation Statistics of Areal Prevalence Rates under High Uncertainty in Denominator Data, Geographical Analysis
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.moran.plot <- function (n, x, se, listw,
                             zero.policy = NULL, spChk = NULL, labels = NULL,
                             xlab = NULL, ylab = NULL, ...)
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

  if (is.null(spChk))
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(n, listw))
    stop("Check of data and weights ID integrity failed")

  labs <- TRUE
  if (is.logical(labels) && !labels)
    labs <- FALSE
  if (is.null(labels) || length(labels) != n1)
    labels <- as.character(attr(listw, "region.id"))

  d <- HCEB.est(n, x, se)

  r <- d$raw # raw rate
  N = sum(n) # sum of the population-at-risk
  m = attributes(d)$m # global average rate

  y <- d$HCEB # HCEB rate

  wy <- lag.listw(listw, y, zero.policy = zero.policy)

  if (is.null(xlab))
    xlab <- "HC-EB Rate"
  if (is.null(ylab))
    ylab <- "spatially lagged HC-EB Rate"

  plot(y, wy, xlab = xlab, ylab = ylab, ...) # Moran Plot - Point

  if (zero.policy) {
    n0 <- wy == 0
    if (any(n0)) {
      symbols(y[n0], wy[n0], inches = FALSE, circles = rep(diff(range(y))/50,
                                                           length(which(n0))), bg = "grey", add = TRUE)
    }
  }

  ywy.lm <- lm((wy-unique(m)) ~ 0+I(y-unique(m))) # Estimateing Fitting line
  abline(coef=c("(Intercept)"=unique(m)-coef(ywy.lm)*unique(m), ywy.lm)) # Drawing fitting line
  abline(h = unique(m), lty = 2) # Global average line
  abline(v = unique(m), lty = 2) # Global average line

  infl.ywy <- influence.measures(ywy.lm)
  is.inf <- which(apply(infl.ywy$is.inf, 1, any))
  points(y[is.inf], wy[is.inf], pch = 9, cex = 1.2)

  if (labs)
    text(y[is.inf], wy[is.inf], labels = labels[is.inf],
         pos = 2, cex = 0.7)
  rownames(infl.ywy$infmat) <- labels
}
environment(HCEB.moran.plot) <- asNamespace("spdep")
