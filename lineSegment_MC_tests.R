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
max.cores <- 30

# import function definitions
source("lineSegmentFunctions.R")

# fallen-log data
logs <- readRDS("output/logs_2021_07_06.rds")
sites <- logs %>% rownames
names(sites) <- sites

# spatstat psp objects
owin.site <- logs$logs.psp[[1]]$window
psp.site <- logs$logs.psp

# number of model simulations 
NSIMS <- 199 

# Simulate restricted resampling and uniform randomisation
randomise_psp <- function(psp_object, randomise, nsims, window = NULL){
  df = psp_object %>% 
    as.data.frame() %>% 
    mutate(angle = angles.psp(psp_object), 
           length = lengths_psp(psp_object))
  c(list(psp_object), lapply(1:nsims, function(i){
    if("angle-rs" %in% randomise) df = df %>% mutate(angle = sample(angle))
    if("angle-unif" %in% randomise) df = df %>% mutate(angle = runif(n())*2*pi)
    if("length-rs" %in% randomise) df = df %>% mutate(length = sample(length))
    if("location-unif" %in% randomise) df = df %>% mutate(x0 = runif(n())*100, y0 = runif(n())*100)
    df %>%
      mutate(x1 = pmin(x0 + length*cos(angle),max(Window(psp_object)$xrange)), 
             y1 = pmin(y0 + length*sin(angle),max(Window(psp_object)$xrange))) %>% 
      mutate(length = NULL, angle = NULL) %>% 
      mutate(x1 = pmax(x1,0), y1 = pmax(y1,0)) %>% 
      as.psp(window= Window(psp_object))
  }))
}

# list of hypotheses by name
r_names <- c("r-angle-rs","r-angle-unif","r-location-unif",
             "r-length-rs", "r-length-rs-location-unif-angle-unif")

# list of restricted quantities for each hypothesis
r_type <- list(c("angle-rs"), c("angle-unif"), c("location-unif"), c("length-rs"), 
                c("angle-rs","length-rs","location-unif"))

names(r_type) <- r_names

# generate Monte Carlo simulations
set.seed(390403) #sample(1e6,1)
sims <- with(logs, r_type %>% map(~ randomise_psp(logs.psp, .x, NSIMS)))

# compute L(r) ----- slow ~ 6 hours  start: 17:56
if(F){
  Lfs <- list()
  for(site in sites){
    Lfs[[site]] <- list()
    for(model in names(sims[[site]])){
      message(paste("Site:", site, "  Model:", model))
      Lfs[[site]][[model]] <- imap(sims[[site]][[model]], ~ {message(.y); Lfibre(.x, max.cores = max.cores)})
    }
    saveRDS(Lfs,paste0(paste0("output/sf.raw_nsim",NSIMS,"_2021_07_08.rds")))
  }
}


# plots
if(F){
# load all summary function data: both sites, all models, all simulations + observations, S(r), K(r), L(r), Kpcf(r)
sf.raw <- readRDS("output/sf.raw_nsim199_2021_07_08.rds")
sf.raw <- readRDS("output/sf.raw_rr2_20210629.rds")

NSIMS <- length(sf.raw[[1]][[1]])-1
r_limits <- c(0,25) # c(r_min, r_max)

# merge r, obs and sims into a single tibble for a given summary function
getSfun <- function(flist, fun = "L"){
  sapply(flist, function(x) tibble(x[[fun]])) %>% 
    as_tibble(.name_repair = ~ c("obs",paste0("sim",1:NSIMS))) %>% bind_cols(r = flist[[1]]$r, .)
}

# colours for plots
blues9 <- brewer.pal(9,"Blues")  
reds9 <- brewer.pal(9,"Reds") 

# main tibble of L(r) data for plots
Lf <- sf.raw %>% lapply(lapply, getSfun, fun = "L")

# plots L function with DCLF progress indicated by the colour of the empirical curve
plotLp<- function(dat, title = NULL, subtitle = NULL){
  dat.long <- pivot_longer(dat,-1,"sim", values_to = "L") 
  dat.obs <- dat %>% select(-r) %>% 
    apply(1, function(row) ((NSIMS + 2) - rank((row - mean(row))^2))/(NSIMS + 1)) %>% 
    t %>% 
    as_tibble %>% 
    mutate(r = dat$r, L = dat$obs, p = ifelse(obs<=0.05, reds9[7],blues9[9])) %>% 
    select(r,L,p, obs)
  dat.mean <- dat %>% select(-r) %>% apply(1,mean) %>% tibble(r = dat$r, L = .)
  dat.long %>% ggplot(aes(r,L-r)) +
    geom_line(aes(group = sim), alpha = 0.2, color = blues9[5], size = 0.5) +
    geom_line(data = dat.mean, col = "grey10", size = 1, lty = "dashed") +
    geom_line(col = dat.obs$p, data = dat.obs, size = 1) +
    #geom_point(data = dat.obs %>% filter(obs <= 0.05), size = 1) +
    labs(title = title, subtitle = subtitle) +
    theme_light()
}

plotLp<- function(dat, r_limits, p_sig = 0.05, title = NULL, subtitle = NULL){
  dat.long <- pivot_longer(dat,-1,"sim", values_to = "L") 

  dat_limited <- dat %>% filter(r >= r_limits[1], r <= r_limits[2])
  dat_sig_p <- dat_limited %>% 
    select(-r) %>% 
    apply(1, function(row) (row - mean(row))^2) %>% # squared area
    apply(1, cumsum) %>% # cumulative area, i.e., integral
    apply(1, function(row) ((NSIMS + 2) - rank(row))/(NSIMS + 1)) %>% # progressive p-value
    {tibble(r = dat_limited$r, p =.["obs",], L = dat_limited$obs)} %>% 
    mutate(sig = p<=p_sig, grp = c(0,diff(sig))) %>% 
    filter(sig) %>% 
    mutate(sig = NULL, grp = cumsum(grp))
  
  dat.mean <- dat %>% select(-r) %>% apply(1,mean) %>% tibble(r = dat$r, L = .)
  dat.long %>% ggplot(aes(r,L-r)) +
    geom_line(aes(group = sim), alpha = 0.2, color = blues9[5], size = 0.5) +
    geom_line(data = dat.mean, col = "grey10", size = 1, lty = "dashed") +
    geom_line(aes(y = obs -r), col = blues9[9], data = dat, size = 1) +
    geom_line(aes(group = grp), col = reds9[7], data = dat_sig_p, size = 1) +
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
Ldiscrepancy <- LfunMean %>% 
  lapply(function(x) x %>% 
           select(-r) %>%
           mutate_at(vars(-obs), ~ discrepancy_fun(obs,.)) %>% 
           filter(obs >= r_limits[1], obs <= r_limits[2]) %>% select(-obs) %>% apply(2,sum))

r_titles_plot <- c("resampled angles", "random uniform angles", "random uniform locations",
                 "resampled lengths", "CSR")
names(r_titles_plot) <- r_names

# plot lists with model discrepancy values in subtitle
sumPlots <- lapply(sites, function(site){
  lapply(r_names, function(x) plotLp(Lf[[site]][[x]], r_limits = r_limits, title = NULL, 
                                              subtitle = paste0(r_titles_plot[x],"   (D = ", round(Ldiscrepancy[[site]][[x]],0),")")))
})

# plot theme
tt <- theme(plot.subtitle = element_text(colour = blues9[9], size = 10), panel.grid = element_blank())
y_lim <- ylim(c(-0.4,3.5))
np <- length(r_names) # number of plots per site

# adjust labels
plots.SX <- sumPlots$`T-SX` %>% map(~ .x + labs(x = NULL) + y_lim + tt)
plots.FR <- sumPlots$`W-FR` %>% map(~ .x + labs(x = NULL) + y_lim + tt)

# collate
plotSX <- ggarrange(plotlist = plots.SX, ncol = 1, nrow = np, labels = LETTERS[1:np]) %>%  
  annotate_figure(top = text_grob("T-SX", color = "grey20", face = "bold", size = 14),
                  bottom = text_grob("r", color = "black", face = "plain", size = 12))
plotFR <- ggarrange(plotlist = plots.FR, ncol = 1, nrow = np, labels = LETTERS[1:np + np]) %>%  
  annotate_figure(top = text_grob("W-FR", color = "grey20", face = "bold", size = 14),
                  bottom = text_grob("r", color = "black", face = "plain", size = 12))

# final plot
Lr_plot <- ggarrange(plotSX, plotFR, nrow = 1, ncol = 2)

#ggsave("plots/results_MC_tests_nsim99_2021_07_08.pdf", plot = Lr_plot, dpi = 300, height = 9, width = 8, device = cairo_pdf)




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
    as_tibble %>% 
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

}
