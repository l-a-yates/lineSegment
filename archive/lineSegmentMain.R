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
library(ggplotify)

rm(list=ls())
select <- dplyr::select

# import function definitions
source("lineSegmentFunctions.R")

max.cores <- 5

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
if(T){
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

# restricted randomisations
randomise_psp <- function(psp_object, randomise, nreps, window){
  df = psp_object %>% 
    as.data.frame() %>% 
    mutate(angle = angles.psp(psp_object), 
           length = lengths_psp(psp_object))
  lapply(1:nreps, function(i){
    df %>% mutate(across(all_of(randomise), sample)) %>% 
      mutate(x1 = pmin(x0 + length*cos(angle),100), y1 = pmin(y0 + length*sin(angle),100)) %>% 
      mutate(length = NULL, angle = NULL) %>% 
      mutate(x1 = pmax(x1,0), y1 = pmax(y1,0)) %>% 
      as.psp(window= window)
  })
}

v <- vector("list",NSIMS)
rrand <- list(`T-SX` = list(ang = v, loc = v, len = v, ))
rrand

randomise_psp(psp.site[[1]], randomise = c("x0","y0"), nreps = 99, window = owin.site)


psp.site[[1]] %>% lengths_psp()

sample(1:10)

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


#-----------------------------------
# Extra plots for model explanations
#-----------------------------------

set.seed(2826) # 2854
#base_inhom <- rpoispp(function(x,y){x/100^2}, win = owin.site)
base_inhom <- rpoispp(function(x,y){exp(-(x)/50)/100}, win = owin.site)
base_hom <- rpoispp(intensity(base_inhom), win = owin.site)
intensity(base_inhom)*10000
# resample empirical lengths
lengths <- psp.site %>% map(lengths_psp) %>% Reduce(c,.)
lengths_inhom <- sample(lengths, npoints(base_inhom))
lengths_hom <- sample(lengths, npoints(base_hom))

# simulate homogeneous and inhomogeneous angular distributions
set.seed(7757)
angles_inhom_vonm <- rvonmises(npoints(base_inhom), circular(0), kappa = circular(2))
angles_hom_vonm<- rvonmises(npoints(base_hom), circular(0), kappa = circular(2))
angles_hom_unif <- rvonmises(npoints(base_hom), circular(0), kappa = circular(0))
angles_inhom_unif <- rvonmises(npoints(base_inhom), circular(0), kappa = circular(0))
#rAngles_hom <- runif(npoints(base_hom))*2*pi

#par(mfrow = c(1, 1),  pty = "s", oma = c(1,1,1,1), mar = c(1,1,1,1)) 

# homogeneous rose plot + pdf
rose_theme <- theme(plot.margin = margin(t= 0, l= -40, r = -10, b = -30))

rose_unif <- as.ggplot(expression({
rose.diag(angles_hom_unif, bins = 24, prop =3.5, ticks = F, axes = F, shrink = 1.5, col = blues9[6],
          control.circle = circle.control(lty = "dashed", col = "grey20"))
curve.circular(dvonmises(x, circular(0), kappa = circular(0))*1.3, 
               xlim = c(0,2*pi), add = T, join = T, lwd = 2, lty = "solid", shrink = 0.8) 
})) + coord_equal() + rose_theme

# inhomogeneous rose plot + pdf
rose_vonm <- as.ggplot(expression({
rose.diag(angles_hom_vonm, bins = 25, prop =3.5, ticks = F, axes = F, shrink = 1.5, col = blues9[6],
          control.circle = circle.control(lty = "dashed", col = "grey20"))
curve.circular(dvonmises(x, circular(0), kappa = circular(1.5))*1.5, 
               xlim = c(0,2*pi), add = T, join = T, lwd = 2, lty = "solid", shrink = 0.8)}
)) + coord_equal() + rose_theme


plot_demo_logs <- function(base_ppp, length, angle, title = NULL, subtitle = NULL){
  base_ppp %>% as.data.frame() %>% 
    rename(X1 =x, Y1 = y) %>% 
    mutate(X2 = pmin(X1 + length*cos(angle),100), Y2 = pmin(Y1 + length*sin(angle),100)) %>% 
    mutate(X2 = pmax(X2,0), Y2 = pmax(Y2,0)) %>% 
    ggplot(aes(x = X1, y = Y1)) +
    geom_segment(aes(xend = X2, yend = Y2), alpha = 1, col = blues9[8], size = 0.5) + #,
    #             arrow = arrow(length = unit(1,"mm"), type = "open")) +
    geom_point(size = 1, col = "grey10") +
    coord_fixed(ratio = 1) +
    theme_bw() +
    labs(x = "" , y = "", title = "", subtitle = subtitle) +
    theme(panel.background = element_rect(fill = "grey99"), 
          panel.grid = element_blank(), axis.text = element_blank(), axis.text.x.top = element_text(),
          axis.ticks = element_blank(), axis.title.x = element_text(margin = margin(t = 10)),
          plot.margin = margin(r = -30, l = -30)) +
    lims(x = c(0,100), y = c(0,100))
}

plots_demo_2 <- list(rose_unif %>% annotate_figure(left = "Uniform"),
                     plot_demo_logs(base_hom, lengths_hom, angles_hom_unif) ,
                     plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_unif),
                     rose_vonm %>% annotate_figure(left = "Non-uniform"),
                     plot_demo_logs(base_hom, lengths_hom, angles_hom_vonm),
                     plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_vonm))

final_demo_plot <- ggarrange(plotlist = plots_demo_2, nrow = 2, ncol = 3, widths = c(0.5,1,1),
                             labels = c("","Homogeneous intensity","Inhomogeneous intensity","","",""),
                             hjust = c(0,-0.32,-0.29,0,0,0), vjust = 1.1,
                             #hjust = c(0,-0.7,-0.6,0,0,0), vjust = 1.1,
                             font.label = list(face = "plain", size = 11)) 

ggsave("demo_logs_v1.pdf", width = 6, height = 5)
#ggsave("plots/demo_logs_v1d.pdf", width = 6, height = 5)


rstan::Rhat()





ddsdsdsdsds


library(scales) # for muted


ggplot(df22, aes(x = x, y = val)) +
  geom_ribbon(aes(ymax = val, ymin = 0, group = type)) +
  geom_col(aes(fill = val)) +
  scale_fill_gradient2(position="bottom" , low = "blue", mid = muted("blue"), high = "red", 
                       midpoint = median(df22$val)) 




plots_demo <- list(plot_demo_logs(base_hom, lengths_hom, angles_hom_unif, title ="hom-unif"),
                   plot_demo_logs(base_hom, lengths_hom, angles_hom_vonm, title ="hom-vonm"),
                   plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_unif, title ="inhom-unif"),
                   plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_vonm, title ="unhom-vonm"),
                   rose_hom, rose_inhom)

ggarrange(plotlist = plots_demo, nrow = 3, ncol = 2, heights = c(1,1,0.5), widths = c(1,1)) %>% 
  annotate_figure(bottom = text_grob("Angles"), left = )

ggsave("demo_plots_.png", width = 6, height = 8)



plots_demo <- list(plot_demo_logs(base_hom, lengths_hom, angles_hom_unif, title ="hom-unif"),
                   plot_demo_logs(base_hom, lengths_hom, angles_hom_vonm, title ="hom-vonm"),
                   plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_unif, title ="inhom-unif"),
                   plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_vonm, title ="unhom-vonm"),
                   rose_hom, rose_inhom)

ggarrange(plotlist = plots_demo, nrow = 3, ncol = 2, heights = c(1,1,0.5), widths = c(1,1)) %>% 
  annotate_figure(bottom = text_grob("Angles"), left = text_grob("Angles"))

ggsave("demo_plots_.png", width = 7, height = 8)


plots_demo_2 <- list(rose_unif %>% annotate_figure(left = "Uniform"),
                   plot_demo_logs(base_hom, lengths_hom, angles_hom_unif),
                   plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_unif),
                   rose_vonm %>% annotate_figure(left = "Non-Uniform"),
                   plot_demo_logs(base_hom, lengths_hom, angles_hom_vonm) +
                                    labs(x = "Homogeneous"),
                   plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_vonm) +
                                    labs(x = "Inhomogeneous"))

final_demo_plot <- ggarrange(plotlist = plots_demo_2, nrow = 2, ncol = 3, widths = c(0.5,1,1),
          labels = c("","A","B","","C","D")) 
#  annotate_figure(bottom = text_grob("Density", rot = 0, face = "bold", size = 14), 
#                  left = text_grob("Angle", rot = 90, size = 14, face = "bold"))

ggsave("demo_plots_2.png", width = 6, height = 5)

annotate()


plots_unif <- cowplot::plot_grid(plot_demo_logs(base_hom, lengths_hom, angles_hom_unif, title ="hom-unif"),
                     plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_unif, title ="inhom-unif"), 
                     rose_hom + theme(plot.margin = margin(l=-15, t = -25, b= -15)) ,
                     rel_heights = c(1,1,0.7),
                     nrow = 3, ncol = 1, greedy = F)

plots_vonm <- cowplot::plot_grid(plot_demo_logs(base_hom, lengths_hom, angles_hom_vonm, title ="hom-unif"),
                                 plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_vonm, title ="inhom-unif"), 
                                 rose_inhom + theme(plot.margin = margin(l=-20, t = -25, b = -25)) ,
                                 rel_heights = c(1,1,0.7),
                                 nrow = 3, ncol = 1, greedy = F)
  
cowplot::plot_grid(plots_unif, plots_vonm) %>% annotate_figure(top = "a")
  
  #cowplot::plot_grid(plotlist = plots_demo, nrow = 3, ncol = 2, greedy = T, align = "h")
  

df <- expand.grid(X1 = 1:10, X2 = 1:10)
df$value <- df$X1 * df$X2

p1 <- ggplot(df, aes(X1, X2)) + geom_tile(aes(fill = value))
p2 <- p1 + geom_point(aes(size = value))

# Basic form
p1 + scale_fill_continuous(guide = "colourbar")
p1 + scale_fill_continuous(guide = guide_colourbar())
p1 + guides(fill = guide_colourbar())

# Control styles

# bar size
leg <- get_legend(p1 + guides(fill = guide_colourbar(barwidth = 0.5, barheight = 22)))


# no label
p1 + guides(fill = guide_colourbar(label = FALSE))

# no tick marks
p1 + guides(fill = guide_colourbar(ticks = FALSE))

# label position
p1 + guides(fill = guide_colourbar(label.position = "left"))

# label theme
p1 + guides(fill = guide_colourbar(label.theme = element_text(colour = "blue", angle = 0)))

# small number of bins
p1 + guides(fill = guide_colourbar(nbin = 3))

# large number of bins
p1 + guides(fill = guide_colourbar(nbin = 100))

# make top- and bottom-most ticks invisible
p1 + scale_fill_continuous(limits = c(0,20), breaks = c(0, 5, 10, 15, 20),
                           guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE))

# guides can be controlled independently
p2 +
  scale_fill_continuous(guide = "colourbar") +
  scale_size(guide = "legend")
p2 + guides(fill = "colourbar", size = "legend")

p2 +
  scale_fill_continuous(guide = guide_colourbar(direction = "horizontal")) +
  scale_size(guide = guide_legend(direction = "vertical"))




tibble(x = seq(0,10,0.01)) %>% 
  ggplot() +
  geom_col(aes(y = 1, x = 1, fill = -x),show.legend = T) +
  scale_fill_gradient(low = blues9[3], high = blues9[8]) +
  theme(legend.position = "bottom")
theme_void() +
  guides()
coord_flip()

col_bar_hom <- tibble() %>% 
  ggplot() +
  geom_col(aes(y = 1, x = 1), fill = blues9[6], show.legend = F) +
  theme_void() +
  coord_flip()

col_bar_inhom
col_bar_hom

plot_demo_logs(base_inhom, lengths_inhom, angles_inhom_vonm) +
  labs(x = "Inhomogeneous")


