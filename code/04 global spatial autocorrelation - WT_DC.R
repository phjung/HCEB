############################
## 4. SPATIAL ASSOCIATION: GLOBAL
############################

library("spdep")
load("data/RData/data_TB.RData")

## Moran plot function
moran.plot.0 <- function (x, listw, m, zero.policy = NULL, spChk = NULL, labels = NULL, 
                          xlab = NULL, ylab = NULL, quiet = NULL, ...) 
{ environment(spdep)
  
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
  if (!is.numeric(x)) 
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(x))) 
    stop("NA in X")
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
  wx <- lag.listw(listw, x, zero.policy = zero.policy)
  if (is.null(xlab)) 
    xlab <- xname
  if (is.null(ylab)) 
    ylab <- paste("spatially lagged", xname)
  plot(x, wx, xlab = xlab, ylab = ylab, ...)
  if (zero.policy) {
    n0 <- wx == 0
    if (any(n0)) {
      symbols(x[n0], wx[n0], inches = FALSE, circles = rep(diff(range(x))/50, 
                                                           length(which(n0))), bg = "grey", add = TRUE)
    }
  }
  # for standardized rate (deviation)
  #xwx.lm <- lm(wx ~ 0+x)
  #abline(xwx.lm)
  #abline(h = 0, lty = 2)
  #abline(v = 0, lty = 2)
  
  # for rate
  xwx.lm <- lm((wx-unique(m)) ~ 0+I(x-unique(m)))
  abline(coef=c("(Intercept)"=unique(m)-coef(xwx.lm)*unique(m), xwx.lm))
  abline(h = unique(m), lty = 2)
  abline(v = unique(m), lty = 2)
  
  infl.xwx <- influence.measures(xwx.lm)
  is.inf <- which(apply(infl.xwx$is.inf, 1, any))
  points(x[is.inf], wx[is.inf], pch = 9, cex = 1.2)
  if (labs) 
    text(x[is.inf], wx[is.inf], labels = labels[is.inf], 
         pos = 2, cex = 0.7)
  rownames(infl.xwx$infmat) <- labels
  if (!quiet) 
    summary(infl.xwx)
  invisible(infl.xwx)
}
environment(moran.plot.0) <- asNamespace("spdep")

############################
## 4.0 SPATIAL WEIGHT MATRIX
############################

census.tract.nb2 <- poly2nb(census.tract, queen=T) ## Queen's distance
census.tract.lw <- nb2listw(census.tract.nb2, style="W")

############################
## 4.1 MORAN PLOT: Teen Birth Rate by Teen Women
############################
## 4.1.1 MORAN PLOT: Teen Birth raw rate by Teen Women
############################

## MORAN PLOT: Teen Birth by Total Birth ## Plot Size: 650 X 610
## Figure 3-1: Moran Plot (Raw rate)

png("output/04 global spatial autocorrelation/teen women/figure_4_1_moran_raw_2010_DC.png", width=7, height=7, unit="in", res=300)
moran.plot(census.tract$rate_raw_WT_DC_2010, census.tract.lw, xlab="Teen Birth Rate", ylab="Spatial Lag") # Moran Plot: Teen Birth Rate - Raw rate (2010)
dev.off()

############################
## 4.1.2 MORAN PLOT: Teen Birth EB rate by Teen Women
############################
## MORAN PLOT: Teen Birth by Total Birth ## Plot Size: 650 X 610
## Figure 3-2: Moran Plot (EB rate)
png("output/04 global spatial autocorrelation/teen women/figure_4_2_moran_EB_2010_DC.png", width=7, height=7, unit="in", res=300)
moran.plot.0(census.tract$rate_EB_WT_DC_2010, census.tract.lw, m=weight_WT$m_WT_2010, xlab="Teen Birth Rate (EB Rate)", ylab="Spatial Lag") #Moran Plot: Teen Birth Rate - EB rate (2010)
dev.off()

############################
## 4.2 MORAN TEST
############################
## 4.2.1 MORAN TEST - Raw rate by Teen Women
############################
moran.std.raw.WT.DC.2010 <- moran.mc(census.tract$rate_std_raw_WT_DC_2010, census.tract.lw, nsim=10000)

############################
## 4.2.2 MORAN TEST - EB rate by Teen Women
############################
moran.std.EB.WT.DC.2010 <- moran.mc(census.tract$rate_std_EB_WT_DC_2010, census.tract.lw, nsim=10000)

############################
## 4.2.4 MORAN TEST - Summary
############################
## SUMMARY

moran.WT.DC <- data.frame(year=2010,
                       moran_raw = c(moran.std.raw.WT.DC.2010$statistic),
                       p_raw     = c(moran.std.raw.WT.DC.2010$p.val),
                       moran_EB  = c(moran.std.EB.WT.DC.2010$statistic),
                       p_EB      = c(moran.std.EB.WT.DC.2010$p.val))

moran.WT.DC.I <- moran.WT.DC %>% gather("variable", "I", c(2, 4)) %>% dplyr::select(year, variable, I)
moran.WT.DC.I$variable <- gsub("moran_raw", "Raw Rate", moran.WT.DC.I$variable)
moran.WT.DC.I$variable <- gsub("moran_EB", "EB Rate", moran.WT.DC.I$variable)
moran.WT.DC.I$variable <- factor(moran.WT.DC.I$variable, levels=c("Raw Rate", "EB Rate", "HC-EB Rate"))

moran.WT.DC.p <- moran.WT.DC %>% gather("variable", "p", c(3, 5)) %>% dplyr::select(year, variable, p)
moran.WT.DC.p$variable <- gsub("p_raw", "Raw Rate", moran.WT.DC.p$variable)
moran.WT.DC.p$variable <- gsub("p_EB", "EB Rate", moran.WT.DC.p$variable)

moran.WT.DC <- full_join(moran.WT.DC.I, moran.WT.DC.p)

moran.WT.DC <- moran.WT.DC %>% dplyr::mutate(p_group=factor(ifelse(p<0.001, "p < 0.001",
                                                             ifelse(p<0.01, "p < 0.01",
                                                                    ifelse(p<0.03, "p < 0.03",
                                                                           ifelse(p<0.05, "p < 0.05",
                                                                                  ifelse(p<0.1, "p  < 0.1", "p > 0.1"))))),
                                                      levels=c("p < 0.001",
                                                               "p < 0.01",
                                                               "p < 0.03",
                                                               "p < 0.05",
                                                               "p < 0.1", "p > 0.1")))

rm(moran.WT.DC.I, moran.WT.DC.p)

## 4.3 EXPORT 
save(list=ls(pattern="moran."), file="output/04 global spatial autocorrelation/teen women/result_Moran_WT_DC.RData") ## save moran's I result
write.csv(moran.WT.DC, file="output/04 global spatial autocorrelation/teen women/table_moran_WT_DC.csv", row.names=F)
############################
rm(list=ls())

