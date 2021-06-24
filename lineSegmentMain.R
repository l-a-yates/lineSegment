#---------------------------------------------------------------------------------
# Spatial pattern analysis of line-segment data in ecology
#
# This code computes the summary statistics L(r) and PCF(r) and produces the 
#  main plots in the manuscript. 
#
# Simulates six Boolean line segment models: 
#  (lambda,theta) = (hom,uniform), (hom,inhom), (inhom,uniform), (inhom,inhom).
#
# Computes S(r), K(r), Kpcf(r), and L(r) for all simulations & observation sets. 
# 
# Analysis applied to Ausplot sites: T-SX, W-FR
#
# Computation of the summary functions L(r) and K(r) are performed with a new  
#  function Lfibre(). This function can be applied to any spatstat psp object
#
# This function Lfibre() makes use of parallel computation. The call to
#  pbmclapply(...) can be replaced with ordinary lapply() if parallel
#  processing is unavailable, but it will be slow.
#
# Uses some results pre-calculated in the accompanying script lineSegmentCircular.R
#
# Authors: Luke Yates, Barry Brook, Jessie Buettel
# File created: 02/03/2019
# Last Edited: 22/12/2020
#---------------------------------------------------------------------------------

library(spatstat)
library(fitdistrplus)
library(circular)
library(pbmcapply)
library(tidyverse)
library(pbapply)
library(latex2exp)
library(RColorBrewer)
library(ggpubr)

rm(list=ls())
select <- dplyr::select

# import function definitions
source("lineSegmentFunctions_v1.R")

max.cores <- 20

# fallen-log data
data.logs <- read.csv("data/logs.csv")
data.site <- split(data.logs,data.logs$SiteName)
sites <- unique(data.logs$SiteName) %>% as.character()
names(sites) <- sites

# gradient data 
dhdx <- readRDS("data/dhdx.rds")
dhdy <- readRDS("data/dhdy.rds")

# local distribution of relative angles. To fit theta = "inhom" models
# theta.par is a 3-vector of (MLE) parameter values for Jones-Pewsey distribution (pre-computed)
theta.par = c(mu = 0, kappa = 0.9339, psi = 0) # distribution for 'large' gradients
slope.thresh <- 0.21 # gradient threshold for large vs small (estimated using profile likelihood)

# create spatstat psp objects
owin.site <- owin(c(0,100), c(0,100)) # observation window
psp.site <- lapply(data.site, function(df) psp(df$X1, df$Y1, df$X2, df$Y2, owin.site))

# simulate models (6 models x 2 sites)
NSIMS <- 99  # global constant

# generate model simulation (fast < 1 min)
if(F){
  sims <- list()
  for (site.name in sites){
    print(site.name)
    sims[[site.name]] <- list()
    sims[[site.name]]$hom.unif = sim.lsp(psp.site[[site.name]], NSIMS, lambda = "hom", theta = "uniform")
    sims[[site.name]]$hom.local = sim.lsp(psp.site[[site.name]], NSIMS, lambda = "hom", theta = "local", theta.par = theta.par, 
                             dhdx = dhdx[[site.name]], dhdy = dhdy[[site.name]],slope.thresh = slope.thresh)
    sims[[site.name]]$hom.resample = sim.lsp(psp.site[[site.name]], NSIMS, lambda = "hom", theta = "resample")
    sims[[site.name]]$inhom.unif = sim.lsp(psp.site[[site.name]], NSIMS, lambda = "inhom", theta = "uniform")
    sims[[site.name]]$inhom.local = sim.lsp(psp.site[[site.name]], NSIMS, lambda = "inhom", theta = "local", theta.par = theta.par, 
                               dhdx = dhdx[[site.name]], dhdy = dhdy[[site.name]], slope.thresh = slope.thresh)
    sims[[site.name]]$inhom.resample = sim.lsp(psp.site[[site.name]], NSIMS, lambda = "inhom", theta = "resample")
  }
}

# evaluate summary functions (slow ~ 20hrs using 20 cores)
if(F){
snum <- 0
Lfs <- list()
  for(site in sites){
    snum <- snum + 1
    Lfs[[site]] <- list()
    SIM_COUNT <- 1
    prog = paste0("Site ", snum, "/2, Model 1/6")
    Lfs[[site]]$hom.unif <- lapply(sims[[site]]$hom.unif, Lfibre,  max.cores = max.cores, prog = prog)
    
    SIM_COUNT <- 1
    prog = paste0("Site ", snum, "/2, Model 2/6")
    Lfs[[site]]$hom.resample <- lapply(sims[[site]]$hom.resample, Lfibre,  max.cores = max.cores, prog = prog)
    
    SIM_COUNT <- 1
    prog = paste0("Site ", snum, "/2, Model 3/6")
    Lfs[[site]]$hom.local <- lapply(sims[[site]]$hom.local, Lfibre,  max.cores = max.cores, prog = prog)
    
    SIM_COUNT <- 1
    prog = paste0("Site ", snum, "/2, Model 4/6")
    Lfs[[site]]$inhom.unif <- lapply(sims[[site]]$inhom.unif, Lfibre,  max.cores = max.cores, prog = prog)
    
    SIM_COUNT <- 1
    prog = paste0("Site ", snum, "/2, Model 5/6")
    Lfs[[site]]$inhom.resample <- lapply(sims[[site]]$inhom.resample, Lfibre,  max.cores = max.cores, prog = prog)
    
    SIM_COUNT <- 1
    prog = paste0("Site ", snum, "/2, Model 6/6")    
    Lfs[[site]]$inhom.local <- lapply(sims[[site]]$inhom.local, Lfibre,  max.cores = max.cores, prog = prog)

    #saveRDS(Lfs,"output/sf.raw_20200819.rds")
  }
}


#--------
# plots
#--------

# load all summary function data: both sites, all models, all simulations + observations, S(r), K(r), L(r), Kpcf(r)
sf.raw <- readRDS("output/sf.raw_20200819.rds")

# merge r, obs and sims into a single tibble for a given summary function
getSfun <- function(flist, fun = "L"){
  sapply(flist, function(x) tibble(x[[fun]])) %>% 
    as_tibble(.name_repair = ~ c("obs",paste0("sim",1:NSIMS))) %>% bind_cols(r = flist[[1]]$r, .)
}

# colours for plots
blues9 <- brewer.pal(9,"Blues")  
reds9 <- brewer.pal(9,"Reds") 

#------------
# L(r) plots
#------------

# main tibble of L(r) data for plots
Lf <- sf.raw %>% lapply(lapply, getSfun, fun = "L")

# plots L function with DCLF progress indicated by the colour of the empirical curve
plotLp<- function(dat, title = NULL, subtitle = NULL){
  dat.long <- pivot_longer(dat,-1,"sim", values_to = "L") 
  dat.obs <- dat %>% select(-r) %>% 
    apply(1, function(row) (101 - rank((row - mean(row))^2))/100) %>% 
    t %>% 
    as.tibble %>% 
    mutate(r = dat$r, L = dat$obs, p = ifelse(obs<=0.05, reds9[7],blues9[9])) %>% 
    select(r,L,p)
  dat.mean <- dat %>% select(-r) %>% apply(1,mean) %>% tibble(r = dat$r, L = .)
  dat.long %>% ggplot(aes(r,L-r)) +
    geom_line(aes(group = sim), alpha = 0.2, color = blues9[5], size = 0.5) +
    geom_line(data = dat.mean, col = "grey10", size = 1, lty = "dashed") +
    geom_line(col = dat.obs$p, data = dat.obs, size = 1) +
    labs(title = title, subtitle = subtitle) +
    theme_light()
}


# compute mean of summary functions
LfunMean <- lapply(sites, function(site){
  Lf[[site]] %>% sapply(function(x) x %>% select(-r) %>% apply(1,mean)) %>% as_tibble() %>% 
    bind_cols(r = Lf[[site]][[1]]$r, obs = Lf[[site]][[1]]$obs, .)
})

# Model Comparison
discrepancy_fun <- function(x, y, q = 1, p =2) (abs(x^q - y^q))^p 
Ldiscrepancy <- LfunMean %>% lapply(function(x) x %>% select(-r) %>% mutate_at(vars(-obs), ~ discrepancy_fun(obs,.)) %>% select(-obs) %>% apply(2,sum))

# AIC (pre-computed)
aicByParts <- list(`T-SX` = c(hom = 1029.8955, inhom = 1007.4468, unif = 330.8179, local = 332.583),
     `W-FR` = c(hom = 663.4451, inhom = 61.1489, unif = 194.8150, local = 179.7831))
aic <- lapply(sites, function(site){
  with(as.list(aicByParts[[site]]), c(hom.unif = hom + unif, hom.resample = NA, hom.local = hom + local,
                                     inhom.unif = inhom + unif, inhom.resample = NA, inhom.local = inhom + local))
})
del.aic <- aic
del.aic$`T-SX`[-c(2,5)] <- aic$`T-SX`[-c(2,5)] - min(aic$`T-SX`[-c(2,5)])
del.aic$`W-FR`[-c(2,5)] <- aic$`W-FR`[-c(2,5)] - min(aic$`W-FR`[-c(2,5)])

# plot lists with model comparison data in subtitle
sumPlots <- lapply(sites, function(site){
  lapply(names(Lf[[site]]), function(x) plotLp(Lf[[site]][[x]], title = paste(site,x), 
                                              subtitle = bquote(Delta*AIC==.(round(del.aic[[site]][[x]],1),)*
                                                                  " "~Discrepancy==.(round(Ldiscrepancy[[site]][[x]],0)))))
})

# plot theme
tt <- theme(plot.title = element_text(colour = blues9[9], size = 12), panel.grid = element_blank())

plots.SX <- list(
  sumPlots$`T-SX`[[1]] + labs(title = "homogeneous-uniform", x = NULL) + tt,
  sumPlots$`T-SX`[[2]] + labs(title = "homogeneous-resample", x = NULL) + tt,
  sumPlots$`T-SX`[[3]] + labs(title = "homogeneous-local", x = NULL) + tt,
  sumPlots$`T-SX`[[4]] + labs(title = "inhomogeneous-uniform", x = NULL) + ylim(c(-0.4,6))+ tt,
  sumPlots$`T-SX`[[5]] + labs(title = "inhomogeneous-resample", x = NULL)+ ylim(c(-1,6)) + tt,
  sumPlots$`T-SX`[[6]] + labs(title = "inhomogeneous-local") + ylim(c(-0.5,6))+ tt)

plots.FR <- list(
  sumPlots$`W-FR`[[1]] + labs(title = "homogeneous-uniform", x = NULL, y = " ") + tt,
  sumPlots$`W-FR`[[2]] + labs(title = "homogeneous-resample", x = NULL, y = " ") + tt,
  sumPlots$`W-FR`[[3]] + labs(title = "homogeneous-local", x = NULL, y = " ") + tt,
  sumPlots$`W-FR`[[4]] + labs(title = "inhomogeneous-uniform", x = NULL, y = " ") + ylim(c(-0.4,6))+ tt,
  sumPlots$`W-FR`[[5]] + labs(title = "inhomogeneous-resample", x = NULL, y = " ")+ ylim(c(-1,6)) + tt,
  sumPlots$`W-FR`[[6]] + labs(title = "inhomogeneous-local", y = " ") + ylim(c(-1,6))+ tt)

plotSX <- ggarrange(plotlist = plots.SX, ncol = 1, nrow = 6, labels = c("A","B","C","D","E","F")) %>%  annotate_figure(top = text_grob("T-SX", color = "grey20", face = "bold", size = 14))
plotFR <- ggarrange(plotlist = plots.FR, ncol = 1, nrow = 6, labels = c("G","H","I","J","K","L")) %>%  annotate_figure(top = text_grob("W-FR", color = "grey20", face = "bold", size = 14))

ggarrange(plotSX, plotFR, nrow = 1, ncol = 2)
#ggsave("Lr_results.eps", dpi = 300, height = 11, width = 8, device = cairo_ps)

#-----------
# pcf plots
#-----------

# requires some variables/objects defined for L(r) plots above

# extract K values suitable for spline-fitting and differentiation
Kpcf <- sf.raw %>% lapply(lapply, getSfun, fun = "Kpcf")

# calls on spatstat::pcf.fv to compute pcf by spline differentiation
computePCF <- function(kvals, r, method = "c"){
  data.frame(r=r, f = kvals) %>% 
    fv(valu = "f") %>% 
    pcf.fv(method = method) %>% 
    as.data.frame() %>% 
    pull(pcf)
}

# plots pcf with DCLF progress indicated by the colour of the empirical curve
plotPCF<- function(dat, method = "c", title = NULL, subtitle = NULL){
  dat <- dat %>% mutate(across(-1, computePCF, r = dat$r, method = method))
  dat.long <- pivot_longer(dat,-1,"sim", values_to = "L") 
  dat.obs <- dat %>% select(-r) %>% 
    apply(1, function(row) (101 - rank((row - mean(row))^2))/100) %>% 
    t %>% 
    as.tibble %>% 
    mutate(r = dat$r, L = dat$obs, p = ifelse(obs<=0.05, reds9[7],blues9[9])) %>% 
    select(r,L,p)
  dat.mean <- dat %>% select(-r) %>% apply(1,mean) %>% tibble(r = dat$r, L = .)
  dat.long %>% ggplot(aes(r,L-1)) +
    geom_line(aes(group = sim), alpha = 0.2, color = blues9[5], size = 0.5) +
    geom_line(data = dat.mean, col = "grey10", size = 1, lty = "dashed") +
    geom_line(col = dat.obs$p, data = dat.obs, size = 1) +
    labs(title = title, subtitle = subtitle, y = "g(r)-1") +
    theme_light()
}

# compute simulation mean
PCFfunMean <- lapply(sites, function(site){
  Kpcf[[site]] %>%
    sapply(function(x) x %>% mutate(across(-1, computePCF, r = x$r)) %>% 
             select(-r) %>% apply(1,mean)) %>% as_tibble() %>% 
    bind_cols(r = Kpcf[[site]][[1]]$r, obs = computePCF(Kpcf[[site]][[1]]$obs, Kpcf[[site]][[1]]$r), .)
})

# Model Comparison
discrepancy_fun <- function(x, y, q = 1, p =2) (abs(x^q - y^q))^p 
discrepancy.pcf <- PCFfunMean %>% lapply(function(x) x %>% select(-r) %>% mutate_at(vars(-obs), ~ discrepancy_fun(obs,.)) %>% select(-obs) %>% apply(2,sum))


# plot lists with model comparison data in subtitle
pcfPlots <- lapply(sites, function(site){
  lapply(names(Kpcf[[site]]), function(x) plotPCF(Kpcf[[site]][[x]], title = paste(site,x), 
                                               subtitle = bquote(Delta*AIC==.(round(del.aic[[site]][[x]],1),)*
                                                                   " "~Discrepancy==.(round(discrepancy.pcf[[site]][[x]],0)))))
})

tt <- theme(plot.title = element_text(colour = blues9[9], size = 12), panel.grid = element_blank())

y_lim <- c(-1,1.5)

plots.pcf.SX <- list(
  pcfPlots$`T-SX`[[1]] + labs(title = "homogeneous-uniform", x = NULL) + lims(y = y_lim) + tt,
  pcfPlots$`T-SX`[[2]] + labs(title = "homogeneous-resample", x = NULL) + ylim(y_lim) + tt,
  pcfPlots$`T-SX`[[3]] + labs(title = "homogeneous-local", x = NULL) + ylim(y_lim) + tt,
  pcfPlots$`T-SX`[[4]] + labs(title = "inhomogeneous-uniform", x = NULL) + ylim(y_lim) + tt,
  pcfPlots$`T-SX`[[5]] + labs(title = "inhomogeneous-resample", x = NULL) + ylim(y_lim) + tt,
  pcfPlots$`T-SX`[[6]] + labs(title = "inhomogeneous-local") + ylim(y_lim) + tt)

y_lim <- c(-1.5,2)

plots.pcf.FR <- list(
  pcfPlots$`W-FR`[[1]] + labs(title = "homogeneous-uniform", x = NULL, y = " ") + lims(y = y_lim) + tt,
  pcfPlots$`W-FR`[[2]] + labs(title = "homogeneous-resample", x = NULL, y = " ") + lims(y = y_lim) + tt,
  pcfPlots$`W-FR`[[3]] + labs(title = "homogeneous-local", x = NULL, y = " ") + lims(y = y_lim) + tt,
  pcfPlots$`W-FR`[[4]] + labs(title = "inhomogeneous-uniform", x = NULL, y = " ") + lims(y = y_lim) + tt,
  pcfPlots$`W-FR`[[5]] + labs(title = "inhomogeneous-resample", x = NULL, y = " ") + lims(y = y_lim) + tt,
  pcfPlots$`W-FR`[[6]] + labs(title = "inhomogeneous-local", y = " ") + lims(y = y_lim) + tt)

plot.pcf.SX <- ggarrange(plotlist = plots.pcf.SX, ncol = 1, nrow = 6, labels = c("A","B","C","D","E","F")) %>%  annotate_figure(top = text_grob("T-SX", color = "grey20", face = "bold", size = 14))
plot.pcf.FR <- ggarrange(plotlist = plots.pcf.FR, ncol = 1, nrow = 6, labels = c("G","H","I","J","K","L")) %>%  annotate_figure(top = text_grob("W-FR", color = "grey20", face = "bold", size = 14))


ggarrange(plot.pcf.SX, plot.pcf.FR, nrow = 1, ncol = 2)
#ggsave("pcf_results.eps", dpi = 300, height = 11, width = 8, device = cairo_ps)

#-----------------
# plot fallen-logs
#-----------------

plotLogs <- function(dat, title = NULL, subtitle = NULL){
  dat %>% 
    ggplot(aes(x = X1, y = Y1)) +
    geom_segment(aes(xend = X2, yend = Y2), alpha = 0.8) +
    geom_point(size = 0.6) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    labs(x = "x", y = "y", title = title, subtitle = subtitle) +
    theme(panel.background = element_rect(fill = "grey95"), panel.grid = element_blank())
}

tt1 <- theme(plot.title = element_text(colour = "grey20", size = 12, face = "bold", hjust = 0.5))

ggarrange(plotLogs(data.site$`T-SX`, "T-SX") + tt1,
          plotLogs(data.site$`W-FR`, "W-FR") + tt1)
#ggsave("fallen_log_plots.eps", dpi = 300, height = 4, width = 8, device = cairo_ps)




