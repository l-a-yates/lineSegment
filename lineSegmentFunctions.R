#-------------------------------------------------------
# Spatial statistics for line segment data 
# Functions for S(r), K(r) and L(r) for fallen-log data 
# Luke Yates
# Last Edited: 14/01/2020
#-------------------------------------------------------


library(spatstat)
library(fitdistrplus)
library(circular)
library(pbmcapply)

# computes S, K and L function for line segment pattern. 
# input: a psp object (spatstat)
# output: data frame of summary function values (Kpcf(r) is suitable for differentiation to g(r))
Lfibre <- function(psp.site.logs, rmax = 25, ..., max.cores = 4){
  on.exit(rm(list = ls()))
  
  # Calculates individual summands for a single line segment at a given radius 
  SFibre.inner <- function(x0,y0,x1,y1,owin.t,r, psp.logs, max.weight = 4){
    on.exit(rm(list = ls()))
    psp.this = psp(x0,y0,x1,y1,owin.t)
    length.this = lengths_psp(psp.this)
    owin.net = intersect.owin(dilation(psp.this, r),owin.site)
    edge.corr = (pi*r^2 + 2*r*length.this)/area.owin(owin.net)
    psp.clipped = clip.psp(psp.logs, owin.net)
    length.total.this = sum(lengths_psp(psp.clipped))-length.this
    edge.corr = min(edge.corr, max.weight)  #set maximum weight
    edge.corr*length.total.this
  }
  
  # Computes S(r) for all segments at fixed 'radius' r
  SFibre.r <- function(r, psp.this){
    on.exit(rm(list = ls()))
    svalues = with(psp.this$ends, mapply(SFibre.inner, x0, y0, x1, y1, list(psp.this$window), list(r), list(psp.this)))
    sum(as.numeric(svalues))
  }
  
  r.values = seq(0.1,rmax,0.1)  
  A = area.owin(psp.site.logs$window)
  L = sum(lengths_psp(psp.site.logs))
  l.mean = mean(lengths_psp(psp.site.logs))
  l.squared.mean = mean(lengths_psp(psp.site.logs)*lengths_psp(psp.site.logs))
  s.values = as.numeric(pbmclapply(r.values, SFibre.r, psp.site.logs, mc.cores = max.cores))

  transformToK = function(Sf,L, n, Lm, LSm){(A/(L - Lm))*(Sf/n + ((2*r.values)/A)*(LSm - Lm^2)) + Lm^2/pi}
  transformToL = function(Sf,L, n, Lm, LSm){sqrt(transformToK(Sf,L, n, Lm, LSm)/pi)-Lm/pi}
  transformToKpcf = function(Sf,L, n, Lm, LSm){(A/(L - Lm))*(Sf/n + ((2*r.values)/A)*(LSm - Lm^2)) - 2*r.values*Lm}
  
  data.frame(r = r.values, S = s.values/psp.site.logs$n,
             K = transformToK(s.values, L, psp.site.logs$n, l.mean, l.squared.mean),
             Kpcf = transformToKpcf(s.values, L, psp.site.logs$n, l.mean, l.squared.mean),
             L = transformToL(s.values, L, psp.site.logs$n, l.mean, l.squared.mean))
} # end Lfibre(...)



#--------------------
# Generates simulations of the Boolean model with parameters fitted to a psp object (psp.site.logs)
# Length distribution is estimated by fitting a log-normal distribution to the edf of all non-truncated logs
# lambda specifies the intensity model of the log bases: "hom" or "inhom"
# theta specifies the angular distribution model for log angles: "uniform", "resample" or "inhom" 
# For theta = "inhom" the three parameter angular distribution is described by the vector theta.par, the slopes
#  dhdx and dhdy must be provided as well as the threshold value (see code and supplementary materials).
# First list element of the returned list contains the original psp object (the observations)
#--------------------
sim.lsp <- function(psp.site.logs, num_reps = 19, lambda = "hom", theta = "uniform", 
                    dhdx = NULL, dhdy = NULL, theta.par = NULL, slope.thresh = NULL){
  #on.exit(rm(list = ls()))
  print(paste("Generating ", num_reps, " simulations of Boolean model"))
  print(paste("lambda:", lambda, "  theta:", theta))
  num_trees = psp.site.logs$n #; print(paste("Num_trees =", num_trees))

  ## Fit log-normal distribution to log lengths
  thres = 100 # side length of square plot
  intersect.boundary = (psp.site.logs$ends$y0 == thres | psp.site.logs$ends$x0 == thres | psp.site.logs$ends$x1 == thres | psp.site.logs$ends$y1 == thres)
  psp.site.logs.whole = subset(psp.site.logs, !intersect.boundary)
  fit.lengths.lnorm = fitdist(lengths_psp(psp.site.logs.whole), "lnorm")
  
  if(lambda == "inhom") logs.ppm = ppm(ppp(psp.site.logs$ends$x0, psp.site.logs$ends$y0, as.owin(psp.site.logs)), ~ polynom(x,y,2))
  
  sims = lapply(1:num_reps, function(i){
    if(lambda == "inhom"){
      sim = simulate(logs.ppm, drop = T)
      rX1 = sim$x
      rY1 = sim$y
    } else {
      rX1 = runif(num_trees, 0, 100)
      rY1 = runif(num_trees, 0, 100)
    }
    
    rAngles = (sample(1:360, length(rX1), replace = T))*pi/180
    
    # gradient-fitted angles
    if(theta == "inhom"){
      log.cell.x = pmax(1,ceiling(rX1/(100/15)))  # determine cell
      log.cell.y = pmax(1,ceiling(rY1/(100/15)))  
      grad_x = sapply(1:length(rX1), function(i) return(dhdx[log.cell.x[i],log.cell.y[i]])) # get grad at cell
      grad_y = sapply(1:length(rX1), function(i) return(dhdy[log.cell.x[i],log.cell.y[i]]))
      slope = sqrt(grad_x^2 + grad_y^2)
      neg.grad.angle = coord2rad(x = -grad_x , y = -grad_y) # compute angle of 'downhill' (negative) gradient
      rel.angles.sim = rvonmises(length(rX1), mu = theta.par['mu'], kappa = theta.par['kappa'])
      rAngles.sim = as.numeric((neg.grad.angle + rel.angles.sim))  # add local gradient to sampled relative angle
      # replace uniform angles with simulated angles when the slope at log base is above the threshold (if provided)
      if(is.null(slope.thresh)) rAngles <- rAngles.sim 
      else rAngles[slope >= slope.thresh] <- rAngles.sim[slope >= slope.thresh]
    }
    
    if(theta == "resample") rAngles <- sample(angles.psp(psp.site.logs, directed=T), length(rX1), replace = T)
    
    rlengths = rlnorm(length(rX1),fit.lengths.lnorm$estimate["meanlog"], fit.lengths.lnorm$estimate["sdlog"])
    
    rX2 = rX1 + rlengths*cos(rAngles); rY2 <- rY1 + rlengths*sin(rAngles) # endpoints
    rX2 <- pmin(rX2,100); rX2 <- pmax(rX2, 0); rY2 <- pmin(rY2,100); rY2 <- pmax(rY2, 0) # truncate logs at owin boundary
    
    psp(round(rX1,1),round(rY1,1), round(rX2,1), round(rY2,1), owin.site)
    
  }) # end lapply
  
  c(list(psp.site.logs), sims)
} # end sim.lsp(...)

