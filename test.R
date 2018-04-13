library("devtools")
library("roxygen2")

devtools::load_all("..")

document()
setwd("HCEB")
setwd("..")
install("HCEB")
setwd("HCEB")

library("HCEB")



library("HCEB")

load("C:/Users/pjung1/Research/2017-06 Public Health/data/RData/data_TB.RData")

attach(census.tract@data)

HCEB.est(TB_2010, WT_2010, MOE_WT_2010, b=0.1)

HCEB.moran(TB_2010, WT_2010, MOE_WT_2010, b=0.1, listw=census.tract.lw)

HCEB.moran.mc(TB_2010, WT_2010, MOE_WT_2010, b=0.1, listw=census.tract.lw,
              sig=0.05, nsim=1000)

HCEB.moran.plot(TB_2010, WT_2010, MOE_WT_2010, b=0.1, listw=census.tract.lw)

HCEB.local(TB_2010, WT_2010, MOE_WT_2010, b=0.1, listw=census.tract.lw)

HCEB.local.mc(TB_2010, WT_2010, MOE_WT_2010, b=0.1, listw=census.tract.lw,
              sig=0.05, nsim=10000)
