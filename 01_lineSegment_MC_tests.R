#---------------------------------------------------------------------------------
# Spatial pattern analysis of line-segment data in ecology
#
# This script computes the summary statistics L(r) and g(r) and produces the 
#  main results figure for the manuscript. 
#
# Computation of the summary functions L(r) and g(r) are performed with a new  
#  function Lfibre(). This function can be applied to any spatstat psp object
#
# The function Lfibre() makes use of parallel computation. The call to
#  pbmclapply(...) can be replaced with ordinary lapply() if parallel
#  processing is unavailable, but it will be slow.
#
# Authors: Luke Yates, Barry Brook, Jessie Buettel
# File created: 02/03/2019
# Last Edited: 12/07/2021
#---------------------------------------------------------------------------------

library(spatstat)
library(pbmcapply)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggplotify)

rm(list=ls())
select <- dplyr::select
max.cores <- 30

# import function definitions
source("00_lineSegmentFunctions.R")

# load fallen-log data
logs <- readRDS("data/logs_2021_07_06.rds")
sites <- logs %>% rownames
names(sites) <- sites

#------------------------
# randomised simulations
#------------------------

# set number of simulations 
NSIMS <- 199 

# list of hypotheses by name
r_names <- c("r-angle-rs","r-angle-unif","r-location-unif",
             "r-length-rs", "r-length-rs-location-unif-angle-unif")

# named list of randomised quantities for each hypothesis
r_type <- list(c("angle-rs"), c("angle-unif"), c("location-unif"), c("length-rs"), 
                c("angle-rs","length-rs","location-unif"))
names(r_type) <- r_names

# generate Monte Carlo simulations
set.seed(390403) #sample(1e6,1)
sims <- with(logs, r_type %>% map(~ randomise_psp(logs.psp, .x, NSIMS)))


#---------------------------------------------------------
# compute summary function estimates ----- slow ~ 6 hours
#---------------------------------------------------------

# compute summary function estimates ----- slow ~ 6 hours
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

#-------
# plots
#-------

if(T){
# load all summary function data: S(r), K(r), L(r), Kpcf(r)
sf.raw <- readRDS("output/sf.raw_nsim199_2021_07_08.rds")

# colours for plots
blues9 <- brewer.pal(9,"Blues")  
reds9 <- brewer.pal(9,"Reds") 

NSIMS <- length(sf.raw[[1]][[1]])-1
r_limits <- c(0,25) # c(r_min, r_max)

# merge r, obs and sims into a single tibble for a given summary function
getSfun <- function(flist, fun = "L"){
  sapply(flist, function(x) tibble(x[[fun]])) %>% 
    as_tibble(.name_repair = ~ c("obs",paste0("sim",1:NSIMS))) %>% bind_cols(r = flist[[1]]$r, .)
}

# main tibble of L(r) data for plots
Lf <- sf.raw %>% lapply(lapply, getSfun, fun = "L")

# plots L function with DCLF progress indicated by the colour of the empirical curve
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
    geom_line(aes(y = obs -r), col = blues9[9], data = dat, size = 0.8) +
    geom_line(aes(group = grp), col = reds9[7], data = dat_sig_p, size = 0.9) +
    labs(title = title, subtitle = subtitle) +
    theme_light()
}

# compute mean of summary functions
LfunMean <- lapply(sites, function(site){
  Lf[[site]] %>% sapply(function(x) x %>% select(-r) %>% apply(1,mean)) %>% as_tibble() %>% 
    bind_cols(r = Lf[[site]][[1]]$r, obs = Lf[[site]][[1]]$obs, .)
})

# total discrepancy on interval r_limits
discrepancy_fun <- function(x, y, q = 1, p =2) (abs(x^q - y^q))^p 
Ldiscrepancy <- LfunMean %>% map(~.x %>% filter(r >= r_limits[1], r <= r_limits[2]) %>% 
                                        select(-r) %>% mutate_at(vars(-obs), ~ discrepancy_fun(obs,.)) %>% 
                                        select(-obs) %>% map(sum))

# titles for plots
r_titles_plot <- c("resampled angles", "random uniform angles", "random uniform locations",
                 "resampled lengths", "CSR"); names(r_titles_plot) <- r_names

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

#ggsave("plots/results_MC_tests_Lfun_nsim199_2021_07_09.pdf", 
#       plot = Lr_plot, dpi = 300, height = 9, width = 8, device = cairo_pdf)


#-----------
# pcf plots
#-----------

# requires some variables/objects defined for L(r) plots above

# extract K values suitable for spline-fitting and differentiation
Kpcf <- sf.raw %>% lapply(lapply, getSfun, fun = "Kpcf")
r_limits_pcf <- c(2,25)

# calls on spatstat::pcf.fv to compute pcf by spline differentiation
computePCF <- function(kvals, r, method = "c"){
  data.frame(r=r, f = kvals) %>% 
    fv(valu = "f") %>% 
    pcf.fv(method = method) %>% 
    as.data.frame() %>% 
    pull(pcf)
}

# plots pcf with DCLF progress indicated by the colour of the empirical curve
plotPCF<- function(dat, method = "c", title = NULL, subtitle = NULL, p_sig = 0.05, r_limits){
  dat <- dat %>% mutate(across(-1, computePCF, r = dat$r, method = method))
  dat.long <- pivot_longer(dat,-1,"sim", values_to = "L") 

  dat_limited <- dat %>% filter(r >= r_limits[1], r <= r_limits[2])
  dat_sig_p <- dat_limited %>% 
    select(-r) %>% 
    apply(1, function(row) (row - mean(row))^2) %>% # squared area
    apply(1, cumsum) %>% # cumulative area, i.e., integral
    apply(1, function(row) ((NSIMS + 2) - rank(row))/(NSIMS + 1)) %>% # progressive p-value
    {tibble(r = dat_limited$r, p =.["obs",], L = dat_limited$obs)} %>% # extract p-values for obs
    mutate(sig = p<=p_sig, grp = c(0,diff(sig))) %>% # boolean accept/reject + indicator for change
    filter(sig) %>% # rejected values only
    mutate(sig = NULL, grp = cumsum(grp)) # grp collates sequential (rejected) r-values into groups
  
  dat.mean <- dat %>% select(-r) %>% apply(1,mean) %>% tibble(r = dat$r, L = .)
  dat.long %>% ggplot(aes(r,L-1)) +
    geom_line(aes(group = sim), alpha = 0.15, color = blues9[5], size = 0.5) +
    geom_line(data = dat.mean, col = "grey10", size = 1, lty = "dashed") +
    #geom_line(col = dat.obs$p, data = dat.obs, size = 1) +
    geom_line(aes(y = obs - 1), col = blues9[9], data = dat, size = 0.8) +
    geom_line(aes(group = grp), col = reds9[7], data = dat_sig_p, size = 0.9) +
    labs(title = title, subtitle = subtitle, y = "g(r)-1") +
    theme_light
}

# compute simulation mean
PCFfunMean <- lapply(sites, function(site){
  Kpcf[[site]] %>%
    sapply(function(x) x %>% mutate(across(-1, computePCF, r = x$r)) %>% 
             select(-r) %>% apply(1,mean)) %>% as_tibble() %>% 
    bind_cols(r = Kpcf[[site]][[1]]$r, obs = computePCF(Kpcf[[site]][[1]]$obs, Kpcf[[site]][[1]]$r), .)
})

# total discrepancy on interval r_limits
discrepancy.pcf <- PCFfunMean %>% map(~.x %>% filter(r >= r_limits_pcf[1], r <= r_limits_pcf[2]) %>% 
  select(-r) %>% mutate_at(vars(-obs), ~ discrepancy_fun(obs,.)) %>% 
  select(-obs) %>% map(sum))

# plot lists with model comparison data in subtitle
pcfPlots <- lapply(sites, function(site){
  lapply(names(Kpcf[[site]]), function(x) plotPCF(Kpcf[[site]][[x]], r_limits = r_limits_pcf,
                                               subtitle = paste0(r_titles_plot[x],"   (D = ", round(discrepancy.pcf[[site]][[x]],0),")")))
})

# adjust labels
plots.pcf.SX <- pcfPlots$`T-SX` %>% map(~ .x + labs(x = NULL) + ylim(c(-1,1.5)) + tt)
plots.pcf.FR <- pcfPlots$`W-FR` %>% map(~ .x + labs(x = NULL) + ylim(c(-1.2,2)) + tt)

# collate
plot.pcf.SX <- ggarrange(plotlist = plots.pcf.SX, ncol = 1, nrow = 5, labels = c("A","B","C","D","E")) %>%  
  annotate_figure(top = text_grob("T-SX", color = "grey20", face = "bold", size = 14),
                  bottom = text_grob("r", color = "black", face = "plain", size = 12))
plot.pcf.FR <- ggarrange(plotlist = plots.pcf.FR, ncol = 1, nrow = 5, labels = c("F","G","H","I","J")) %>%  
  annotate_figure(top = text_grob("W-FR", color = "grey20", face = "bold", size = 14),
                  bottom = text_grob("r", color = "black", face = "plain", size = 12))


pcf_plot_main <- ggarrange(plot.pcf.SX, plot.pcf.FR, nrow = 1, ncol = 2)
#ggsave("pcf_results.eps", dpi = 300, height = 11, width = 8, device = cairo_ps)
#ggsave(plot = pcf_plot_main, "plots/results_MC_tests_pcf_nsim199_2021_07_09.pdf", 
#       dpi = 300, height = 9, width = 8, device = cairo_pdf)
}
