#-------------------------------------------------------
# Spatial statistics for line segment data 
#
# Main function to compute summary statistics for a line-segment pattern
#
# The function Lfibre() can be applied to any spatstat psp object
#
# Authors: Luke Yates
# File created: 02/03/2019
# Last Edited: 10/07/2021
#-------------------------------------------------------

library(spatstat)
library(pbmcapply)

# computes S, K and L function for line segment pattern. 
# input: a psp object (spatstat)
# output: data frame of summary function values (Kpcf is suitable for differentiation to g(r))
Lfibre <- function(psp_object, r_min = 0, r_max = 25, r_increment = 0.1, ..., max.cores = 4){
  on.exit(rm(list = ls()))
  # Calculates individual summands for a single line segment at a given radius r
  SFibre.inner <- function(x0,y0,x1,y1,r,max.weight = 4){
    on.exit(rm(list = ls()))
    psp.i = psp(x0,y0,x1,y1,W)
    length.i = lengths_psp(psp.i)
    owin.net = intersect.owin(dilation(psp.i, r),W)
    edge.corr = (pi*r^2 + 2*r*length.i)/area.owin(owin.net)
    length.total.i = sum(lengths_psp(clip.psp(psp_object, owin.net)))-length.i
    edge.corr = min(edge.corr, max.weight)  #set maximum weight
    edge.corr*length.total.i
  }
  # Computes S(r) for all segments at fixed 'radius' r
  SFibre.r <- function(r){
    on.exit(rm(list = ls()))
    sum(as.numeric(with(psp_object$ends, mapply(SFibre.inner, x0, y0, x1, y1, list(r)))))
  }
  r = seq(r_min + r_increment,r_max,r_increment) # vector of r-values
  W <- psp_object$window
  A = area.owin(W)
  L = sum(lengths_psp(psp_object))
  n = psp_object$n
  Lm = mean(lengths_psp(psp_object)) # mean length
  LSm = mean(lengths_psp(psp_object)^2) # mean of squared lengths
  Sf = as.numeric(pbmclapply(r, SFibre.r, mc.cores = max.cores))/n # eq(1)
  transformToK = function(){(A/(L - Lm))*(Sf + ((2*r)/A)*(LSm - Lm^2)) + Lm^2/pi} # eq(5)
  transformToL = function(){sqrt(transformToK()/pi)-Lm/pi} # eq(6)
  transformToKpcf = function(){(A/(L - Lm))*(Sf + ((2*r)/A)*(LSm - Lm^2)) - 2*r*Lm} # eq(4)
  data.frame(r = r, 
             S = Sf,
             K = transformToK(),
             Kpcf = transformToKpcf(),
             L = transformToL())
} # end Lfibre(...)




# simulates restricted resampling and uniform randomisation
randomise_psp <- function(psp_object, randomise, nsims){
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
} # end randomise_psp()
