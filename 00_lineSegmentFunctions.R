#-------------------------------------------------------
# Spatial statistics for line segment data 
# Functions for S(r), K(r) and L(r) for fallen-log data 
# Luke Yates
# Last Edited: 10/06/2021
#-------------------------------------------------------

library(spatstat)
library(pbmcapply)

# computes S, K and L function for line segment pattern. 
# input: a psp object (spatstat)
# output: data frame of summary function values (Kpcf(r) is suitable for differentiation to g(r))
Lfibre <- function(psp_object, rmax = 25, ..., max.cores = 4){
  on.exit(rm(list = ls()))
  
  # Calculates individual summands for a single line segment at a given radius 
  SFibre.inner <- function(x0,y0,x1,y1,owin.t,r, psp_object, max.weight = 4){
    on.exit(rm(list = ls()))
    psp.i = psp(x0,y0,x1,y1,owin.t)
    length.this = lengths_psp(psp.i)
    owin.net = intersect.owin(dilation(psp.i, r),owin.site)
    edge.corr = (pi*r^2 + 2*r*length.this)/area.owin(owin.net)
    psp.clipped = clip.psp(psp_object, owin.net)
    length.total.this = sum(lengths_psp(psp.clipped))-length.this
    edge.corr = min(edge.corr, max.weight)  #set maximum weight
    edge.corr*length.total.this
  }
  
  # Computes S(r) for all segments at fixed 'radius' r
  SFibre.r <- function(r, psp_object){
    on.exit(rm(list = ls()))
    svalues = with(psp_object$ends, mapply(SFibre.inner, x0, y0, x1, y1, list(psp_object$window), list(r), list(psp_object)))
    sum(as.numeric(svalues))
  }
  
  r.values = seq(0.1,rmax,0.1)  
  A = area.owin(psp_object$window)
  L = sum(lengths_psp(psp_object))
  l.mean = mean(lengths_psp(psp_object))
  l.squared.mean = mean(lengths_psp(psp_object)*lengths_psp(psp_object))
  s.values = as.numeric(pbmclapply(r.values, SFibre.r, psp_object, mc.cores = max.cores))

  transformToK = function(Sf,L, n, Lm, LSm){(A/(L - Lm))*(Sf/n + ((2*r.values)/A)*(LSm - Lm^2)) + Lm^2/pi} # eq(5)
  transformToL = function(Sf,L, n, Lm, LSm){sqrt(transformToK(Sf,L, n, Lm, LSm)/pi)-Lm/pi} # eq(6)
  transformToKpcf = function(Sf,L, n, Lm, LSm){(A/(L - Lm))*(Sf/n + ((2*r.values)/A)*(LSm - Lm^2)) - 2*r.values*Lm} # eq(4)
  
  data.frame(r = r.values, S = s.values/psp_object$n,
             K = transformToK(s.values, L, psp_object$n, l.mean, l.squared.mean),
             Kpcf = transformToKpcf(s.values, L, psp_object$n, l.mean, l.squared.mean),
             L = transformToL(s.values, L, psp_object$n, l.mean, l.squared.mean))
} # end Lfibre(...)
