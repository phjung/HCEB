## Compute HCEB Rate - DONE

HCEB.est <- function(x, n, k, b=0)
{
  # x: a numeric vector of counts of cases
  # n: a numeric vector of populations at risk
  # k: a numeric vector of margin of errors of x
  # b: base number of denominator

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
  if (any(x < .Machine$double.eps))
    stop("non-positive risk population")
  if (any(n < 0))
    stop("negative number of cases")
  if (any(x == 0))
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

## Moran Plot - HCEB Moran Plot - Done

moran.plot.HCEB <- function (x, n, k, b=0, listw,
                             zero.policy = NULL, spChk = NULL, labels = NULL)
{
  r <- x/(b+n) # raw rate
  N = sum(b+n) # sum of the population-at-risk
  m = sum(x)/N # global average rate
  y <- HCEB.est(x, n, k, b)$rate_HCEB # HCEB rate

  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(quiet))
    quiet <- !get("verbose", envir = .spdepOptions)
  stopifnot(is.vector(x))
  stopifnot(is.logical(quiet))
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))

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
    xlab <- "HCEB Rate"
  if (is.null(ylab))
    ylab <- "spatially lagged HCEB Rate"

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
  if (!quiet)
    summary(infl.ywy)
  invisible(infl.ywy)
}
environment(moran.plot.HCEB) <- asNamespace("spdep")

## Global Moran's I - Compute HCEB Global Moran's I  - Done

HCEB.moran <- function (x, n, k, b=0, listw, zero.policy = NULL, NAOK = FALSE)
{
  # x: a numeric vector of counts of cases
  # n: a numeric vector of populations at risk
  # k: a numeric vector of margin of errors of x
  # b: base number of denominator
  # zero.policy   default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
  # NAOK	        if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.

  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  if (n1 != length(x) | n1!=length(n))
    stop("objects of different length")

  r <- x/(b+n)  # raw rate
  N = sum(b+n)  # sum of the population-at-risk
  m = sum(x)/N  # global average rate

  z <- HCEB.est(x, n, k, b)$rate_HCEB-m
  lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  zz <- sum(z^2)

  I <- (sum(z * lz, na.rm = NAOK))/zz

  res <- list(I = I)
  res
}
environment(HCEB.moran) <- asNamespace("spdep")

## Global Moran's I - HCEB Global Moran's I Permuatation Test

moran.mc(census.tract$rate_std_EB2_WT_2010, listw, nsim)

HCEB.moran.mc <- function(x, n, k, b=0, listw, sig = 0.05, nsim,
                          zero.policy = NULL, NAOK = FALSE)
{
  I.global <- HCEB.moran(x, n, k, b, listw, zero.policy, NAOK)
  n1 <- length(listw$neighbours)
  I.perm <- numeric(length = nsim)

  for (i in 1:nsim) {
    perm <- sample(n1, replace=F)
    I.global.perm <- HCEB.moran(x[perm], n[perm], k[perm], b, listw, zero.policy, NAOK)
    I.perm[i,] <- I.local.perm
    rm(I.global.perm, perm)
    }

  res$p.perm <- sum(t(I.perm)>I.global$I)/nsim

  res
}


## Local Moran's I - Compute HCEB Local Moran's I - Done

HCEB.local <- function(x, n, k, b=0, listw,
                       zero.policy = NULL, NAOK = FALSE)
{
  # x: a numeric vector of counts of cases
  # n: a numeric vector of populations at risk
  # k: a numeric vector of margin of errors of x
  # b: base number of denominator
  # zero.policy   default NULL, use global option value; if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA
  # NAOK	        if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed on to the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or 'Inf' values is regarded as an error.

  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours) # Number of areas
  if (n1 != length(x) | n1!=length(n))
    stop("objects of different length")


  r <- x/(b+n)  # raw rate
  N = sum(b+n)  # sum of the population-at-risk
  m = sum(x)/N  # global average rate

  z <- HCEB.est(x, n, k, b)$rate_HCEB-m
  lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  s2 <- sum(z^2, na.rm = NAOK)/n1

  res[, 1] <- (z/s2) * lz

}

## Local Moran's I - LISA Map - DONE

HCEB.LISA.mc <- function(x, n, k, b=0, listw, sig = 0.05, nsim,
                      zero.policy = NULL, NAOK = FALSE)
{
  # x is a vector of the values on which to calculate the MoranI statistic
  # R, lw, ... are all the arguments passed to the localmoran function

  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours) # Number of areas
  if (n1 != length(x) | n1!=length(n))
    stop("objects of different length")

  I.local <- HCEB.local(x, n, k, b, listw, zero.policy, NAOK)
  z <- HCEB.est(x, n, k, b)$rate_HCEB-m
  lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)

  I.perm <- matrix(data=NA, nrow=nsim, ncol=length(x))
  for(i in 1:nsim){
    perm <- sample(n1, replace=F)
    I.local.perm <- HCEB.local(x[perm], n[perm], k[perm], b, listw, zero.policy, NAOK)
    I.perm[i,] <- I.local.perm
    rm(I.local.perm, perm)
  }

  res <- data.frame(p.perm=NA, quad_sig=NA)

  res$p.perm <- rowSums(t(I.perm)>I.local[,1])/nsim
  res$quad_sig <- NA

  res[(z >= 0 & lz >= 0) & (res$p.perm <= p), "quad_sig"] <- "High-High"
  res[(z <= 0 & lz <= 0) & (res$p.perm <= p), "quad_sig"] <- "Low-Low"
  res[(z >= 0 & lz <= 0) & (res$p.perm <= p), "quad_sig"] <- "High-Low"
  res[(z <= 0 & lz >= 0) & (res$p.perm <= p), "quad_sig"] <- "Low-High"

}
