#' Moran scatterplot with Heteroscedastic Consistent Empirical Bayes rate
#'
#' A plot of HC-EB spatial areal rate data against its spatially lagged values.
#' @param x a numeric vector of counts of cases
#' @param n a numeric vector of populations at risk
#' @param k a numeric vector of margin of errors of x
#' @param b default 0, base populations-at-risk
#' @param listw a listw object created for example by nb2listw
#' @param zero.policy default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
#' @param spChk should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use get.spChkOption()
#' @param labels character labels for points with high influence measures, if set to FALSE, no labels are plotted for points with large influence
#' @references #####, ##### and #####. 2018. Small area adjustment of spatial autocorrelation estimators for heteroscedastic sampling errors, Geographical Analysis (under review)
#' @keywords HCEB
#' @export
# @examples
# HCEB.est()

HCEB.moran.plot <- function (x, n, k, b=0, listw,
                             zero.policy = NULL, spChk = NULL, labels = NULL,
                             xlab = NULL, ylab = NULL, ...)
{
  r <- x/(b+n) # raw rate
  N = sum(b+n) # sum of the population-at-risk
  m = sum(x)/N # global average rate
  y <- HCEB.est(x, n, k, b)$HCEB # HCEB rate

  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  xname <- deparse(substitute(x))
  nname <- deparse(substitute(n))
  kname <- deparse(substitute(k))

  if (!is.numeric(x))
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(x)))
    stop("NA in case-count vector")
  if (!is.numeric(n))
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(n)))
    stop("NA in population-at-risk vector")
  if (!is.numeric(k))
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(k)))
    stop("NA in Margin of Error vector")

  n <- length(listw$neighbours)

  if (n != length(x))
    stop("objects of different length")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  labs <- TRUE
  if (is.logical(labels) && !labels)
    labs <- FALSE
  if (is.null(labels) || length(labels) != n)
    labels <- as.character(attr(listw, "region.id"))

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
