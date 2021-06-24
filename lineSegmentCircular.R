#-----------------------------------------------------------------------------------------------------
#
# Analysis of fallen log data for 12 Ausplot sites 
# This script:
# 1) Converts 15x15 grid of DEM data to 100x100 slope data for each site (loess interpolated)
# 2) Computes relative angle between each log and the downhill gradient at its base
# 3) Fits and compare various circular distributions using combined set of relative angles
# 4) Exports hyperframe with psp, dhdx, dhdy, and relAngle columns
#
# Results: 
# 1) Trees tend to fall downhill (known and confirmed here)
# 2) The von Mises distribution is the AIC best model: mu = 0, kappa = 0.436, psi = 0
#
# Author: Luke Yates
# Created: 27th Mar 2019
# Last Edit: 30th July 2020
# Previously Kfibrefibre_v4_rel_angle.R
#
#-------------------------------------------  ----------------------------------------------------------

library(spatstat)
library(tidyverse)
library(circular)
library(pbmcapply)

rm(list=(ls()))

select <-  dplyr::select

#----------
# Functions
#----------
#----

# initial loess smoothing followed by spline fitting
# spline allows derivatives, i.e. slope
# returns a function
loess.fit <- function(h, deriv = 0, ncell = 15){
  x <- seq(100/(2*ncell), by = 100/ncell, length.out = ncell)
  f.loess <-loess(h~x, control = loess.control(surface = "direct"), degree = 2)
  quartic.fun <- splinefun(x = x, y = predict(f.loess), method = "natural", ties = mean)
  return(quartic.fun(x, deriv=deriv))
}

# 15x15 matrix converted to x,y,x tibble
# original values are centred on each cells w.r.t 100m x 100m dimensions
# Matrix must be ordered such that mat[cell_x,cell_y] are the cell references with 1,1 being the bottom left
matrixToTib <- function(mat, ncell = 15){
  tib <- mat %>% as_tibble()
  names(tib) <- as.character(round(seq(100/(2*ncell), by = 100/ncell, length.out = ncell),2))
  tib %>% bind_cols(y = round(seq(100/(2*ncell), by = 100/ncell, length.out = ncell),2),.) %>% 
    gather("x","z",-1) %>% 
    select(x,y,z) %>% 
    mutate(x = as.double(x)) %>%  
    arrange(x,y)
}

# interpolate fine scale DEM (100 x 100)
dem100asMatrix <- function(dem.tb) lapply(dem.tb, function(dem.tib){
  dem.loess = loess(z ~ x * y, data = dem.tib, control = loess.control(surface = "direct"), degree = 2)
  m1 = bind_cols(expand_grid(x = 1:100-0.5, y = 1:100-0.5), 
                 z = predict(dem.loess, newdata = expand_grid(x = 1:100-0.5, y = 1:100-0.5))) %>%
    pivot_wider(names_from = y, values_from = z) %>% as.matrix()
  m <- m1[,-1]; rownames(m) <- m1[,1]
  m
})


#------
# Main
#------

# load auxillary functions for the circular package
load("data/CircStatsInR.RData") 

# load log data
data.logs.all <- read.csv("data/logs.csv")
data.logs <- data.logs.all %>% split(data.logs.all$SiteName)

sites <- unique(data.logs.all$SiteName) %>% as.character()
names(sites) <- sites

# import and format DEM data
dem.raw <- lapply(sites, function(code) read_csv(paste0("data/",code,"-dem-grid-lattice.csv")))
dem.mat <- lapply(dem.raw, function(mat) mat[15:1,-1]) # remove row names and set mat[i,j] -> mat[x,y] (Cartesian)
dem.tb <- lapply(dem.mat, matrixToTib) # convert to tibble with X,Y,Z columns

# interpolate fine scale DEM (100 x 100)
dem.100 <- dem100asMatrix(dem.tb)

# calculate fine scale slope in x and y directions
dhdx = lapply(dem.100, function(data.dem) apply(data.dem, MARGIN = 2, loess.fit, deriv = 1, ncell = 100))
dhdy = lapply(dem.100, function(data.dem) t(apply(data.dem, MARGIN = 1, loess.fit, deriv = 1, ncell = 100)))

# create spatstat psp objects
owin.site <- owin(c(0,100), c(0,100))
psp.site <- lapply(data.logs, function(df) psp(df$X1, df$Y1, df$X2, df$Y2, owin.site))

# compute relative angles of logs to negative gradient at log base
# slope.thresh determines slope threshold at base for inclusion in returned set
# slope.greater determines whether retuend values should be above (T) or below(F) the threshold value
getRelAngles <- function(slope.thresh = 0, slope.greater = T, returnSlope = F) {lapply(sites, function(site){
    logs.dir.rad <- with(data.logs[[site]], coord2rad(x = X2-X1, y = Y2-Y1)) 
    ## associate log base with the nearest cell in the 100X100 array
    log.cell.x <- pmax(1,ceiling(data.logs[[site]]$X1))
    log.cell.y <- pmax(1,ceiling(data.logs[[site]]$Y1))
    ## extract x-y components of the gradient at log base
    grad_x <- sapply(1:nrow(data.logs[[site]]), function(i) return(dhdx[[site]][log.cell.x[i],log.cell.y[i]]))
    grad_y <- sapply(1:nrow(data.logs[[site]]), function(i) return(dhdy[[site]][log.cell.x[i],log.cell.y[i]]))
    slope <- sqrt(grad_x^2 + grad_y^2)
    if(returnSlope) return(slope)
    ## negative gradient at log bases
    neg.grad.at.base.rad <- coord2rad(x = -grad_x , y = -grad_y)
    ## relative angles
    angles.relative <- logs.dir.rad - neg.grad.at.base.rad
    if(slope.greater) angles.relative <- angles.relative[slope >= slope.thresh]
    else angles.relative <- angles.relative[slope < slope.thresh]
    angles.relative.circ <- minusPiPlusPi(circular(angles.relative, type = "angles", units = "radians", zero = 0))
    as.numeric(angles.relative.circ)
  })
}

Reduce(c,getRelAngles()) %>% density %>% plot

# create hyperframe
logs <- hyperframe(logs.psp = psp.site, row.names = names(psp.site))
logs$dhdx <- dhdx[rownames(logs)]
logs$dhdy <- dhdy[rownames(logs)]
logs$relAngle <- getRelAngles()
#saveRDS(logs,"../logs_20200731.rds")


# combine and plot all relative angles
all.rel.angles <- Reduce(c,getRelAngles())
all.rel.angles.circ <- minusPiPlusPi(circular(all.rel.angles, type = "angles", units = "radians", zero = 0))
plot(all.rel.angles.circ, cex = 0.5, bin = 45, stack = T, shrink = 1.5, main = "Relative angles - All 12 sites", axes = F, col = "red")
axis.circular(at=NULL, labels = c("0","+90","180","-90"), template = "none", zero =0, tick=TRUE, col = "blue")
#rose.diag(all.rel.angles.circ, bins = 32, col = "darkgrey", prop = 3, add = T, axes = F)

# fit and plot von Mises distribution
fit.vonm <- mle.vonmises(all.rel.angles %% (2*pi))
curve.circular(dvonmises(x, mu =fit.vonm$mu, kappa = fit.vonm$kappa), xlim = c(-1,2), add = T, join = T, lwd = 1)
curve.circular(dunif(x, -2*pi, 2*pi), xlim = c(-1,2), add = T, join = T, lwd = 1)
vMPPQQ(all.rel.angles %% (2*pi),fit.vonm$mu, fit.vonm$kappa) ## QQ plots etc.

## Model comparison. Compute AIC and BIC Jones-Pewson relative to submodels:
## Jones-Pewsey (estimate psi), Cauchy (psi = -1), Cardoid (psi = 1) and von Mises (psi = 0)
AIC_BIC_JP <- JPpsi0AICBIC(all.rel.angles %% (2*pi),0)[c(2,4)]
AIC_BIC_vM <- JPpsi0AICBIC(all.rel.angles %% (2*pi),0)[c(1,3)]
AIC_BIC_Card <- JPpsi0AICBIC(all.rel.angles %% (2*pi),1)[c(1,3)]
AIC_BIC_Cauchy <- JPpsi0AICBIC(all.rel.angles %% (2*pi),-1)[c(1,3)]
compare <- as.matrix(cbind(AIC_BIC_JP,AIC_BIC_vM,AIC_BIC_Card,AIC_BIC_Cauchy))
dimnames(compare) <- list(c("AIC","BIC"), c("JP","vM","Cardoid","Cauchy"))
compare
compare[1,] %>% simplify - min(compare[1,] %>% simplify) # von Mises is the best fit



# force mu = 0 for a simpler model that is centred around the downhill direction
fit.vonm1 <- mle.vonmises(all.rel.angles %% (2*pi), mu = 0)
# mu = 0 has a better AIC (# kappa = 0.436)
sum(dvonmises(all.rel.angles %% (2*pi), fit.vonm1$mu, fit.vonm1$kappa, log = T))*-2 + 2*1
sum(dvonmises(all.rel.angles %% (2*pi), fit.vonm$mu, fit.vonm$kappa, log = T))*-2 + 2*2

# check wrapped normal: It has an inferior AIC to von M, ignore
fit.wrap <- mle.wrappednormal(all.rel.angles %% (2*pi))
sum(log(dwrappednormal(all.rel.angles %% (2*pi), fit.wrap$mu, fit.wrap$rho, fit.wrap$sd)))*-2 + 2*3


# determine threshold for local effect by profile log-likelihood
# use uniform distribution below threshold and von Mises distribution (mu = 0) above threshold

# profile ll for both sites together
thresh.vec <- seq(0.01,0.3,0.005)
ll.thresh <- sapply(thresh.vec, function(thresh){
  ll1 = sapply(getRelAngles(thresh,F), function(angles) dunif(angles, min = -1*pi, max = pi, log = T) %>% sum) %>% sum
  angles2 = circular(Reduce(c,getRelAngles(thresh,T)))
  fitvm = mle.vonmises(angles2, mu = circular(0))
  ll2 = sum(dvonmises(angles2, circular(fitvm$mu), circular(fitvm$kappa), log = T))
  ll1 + ll2
})
thresh.vec[which.max(ll.thresh)] # gradient threshold = 0.21
max(ll.thresh) # -252.183
threshValue <- 0.21 # N.B. Same effective threshold is found when computed for each site separately.

# vm fit for slopes > threshold, for both sites
fit.vm <- mle.vonmises(circular(Reduce(c,getRelAngles(threshValue,T))),mu = circular(0))

# AIC: T-SX
ll1 <- dunif(getRelAngles(threshValue,F)$`T-SX`, min = -1*pi, max = pi, log = T) %>% sum
angles2 = circular(getRelAngles(threshValue,T)$`T-SX`)
ll2 = sum(dvonmises(angles2, circular(fit.vm$mu), circular(fit.vm$kappa), log = T))
ll1 + ll2 # ll = -164.2915, AIC = -2*ll + 2*2 = 332.583  (params = kappa, threshold)

# AIC: W-FR (All slopes are greater than threshold, i.e. ll1 = 0 and local dist applies to all logs)
ll1 <- dunif(getRelAngles(threshValue,F)$`W-FR`, min = -1*pi, max = pi, log = T) %>% sum
angles2 = circular(getRelAngles(threshValue,T)$`W-FR`)
ll2 = sum(dvonmises(angles2, circular(fit.vm$mu), circular(fit.vm$kappa), log = T))
ll1 +ll2 # -87.89153, AIC = -2*ll + 2*2 = 179.7831 (params = kappa, threshold)


#-------
# plots
#-------

threshValue <- 0.21   #

setEPS()
postscript("plots/circDensity.eps", height = 4, width = 8)
par(mfrow = c(1, 2),  pty = "s", oma = c(1,1,1,1), mar = c(1,1,1,1)) 

angles1 <- Reduce(c,getRelAngles(threshValue,F)) %>% circular()
plot(angles1, cex = 0.5, bin = 15, stack = T, shrink = 1.5, main = "Small slopes", axes = F, col = "dodgerblue")
axis.circular(at=NULL, labels = c("0","+90","180","-90"), template = "none", zero =0, tick=TRUE, col = "black", cex = 1)
curve.circular(dunif(x, -2*pi, 2*pi)*1, xlim = c(-1,2), add = T, join = T, lwd = 1.5)

angles2 <- Reduce(c,getRelAngles(threshValue,T)) %>% circular()
plot(angles2, cex = 0.5, bin = 45, stack = T, shrink = 1.5, main = "Large slopes", axes = F, col = "dodgerblue")
axis.circular(at=NULL, labels = c("0","+90","180","-90"), template = "none", zero =0, tick=TRUE, col = "black")
fit.vm = mle.vonmises(angles2, mu = circular(0))
curve.circular(dvonmises(x, mu =fit.vm$mu, kappa = fit.vm$kappa), xlim = c(-1,2), add = T, join = T, lwd = 1.5)

dev.off()

warnings()


