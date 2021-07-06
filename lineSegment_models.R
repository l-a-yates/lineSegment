
library(spatstat)
library(tidyverse)
library(circular)
library(pbmcapply)

rm(list=(ls()))

select <-  dplyr::select

# load extra functions from use with circular package
#load("data/CircStatsInR.RData", envir = circ <- new.env()) 

# load spatstat::hyperframe object
logs <- readRDS("output/logs_2021_07_06.rds")
logs$angle <- with(logs, as.data.frame(logs.psp) %>% {coord2rad(.$x1-.$x0, .$y1-.$y0)})
logs$slope <- with(logs, sqrt(dhdx^2 + dhdy^2))
logs$slope.im <- with(logs, as.im(slope, W = Window(logs.psp)))

#--------------------
# 1: intensity models
#--------------------

# define point pattern of log bases
logs$bases.ppp <- with(logs, as.data.frame(logs.psp) %>% ppp(x = .$x0, y = .$y0, window = Window(logs.psp)))

# fit all models
m_intensity <- list()
m_intensity$hom <- with(logs, ppm(bases.ppp ~ 1))
m_intensity$slope <- with(logs, ppm(bases.ppp ~ 1 + slope.im))
m_intensity$xy <- with(logs, ppm(bases.ppp ~ 1 + polynom(x,y,2)))
m_intensity$xy_slope <- with(logs, ppm(bases.ppp ~ 1 + polynom(x,y,2) + slope.im))

# model comparison
m_intensity %>% map("T-SX") %>% map_dbl(AIC) %>% {. - min(.)} %>% round(3)
m_intensity %>% map("W-FR") %>% map_dbl(AIC) %>% {. - min(.)} %>% round(3)


#--------------------
# 2: angle models
#--------------------

# fits the von Mises circular distribution to supplied set of angles
# set mu and/or kappa equal to zero for simpler submodels
# returns the fitted object, the log likelihood, and the AIC estimate
fit_vonmises <- function(angles, mu = NULL, kappa = NULL){
  fit <- mle.vonmises(angles, mu = mu, kappa = kappa)
  loglik <- dvonmises(angles, fit$mu, fit$kappa, log = T) %>% sum
  K <- as.numeric(is.null(mu)) + as.numeric(is.null(kappa))
  aic <- -2*loglik + 2*K
  list(fit = fit, loglik = loglik, aic = aic)
}

# calculate slope at each log base and add to hyperframe
logs$slope_i <- with(logs, {
  bases.ppp %>% as.data.frame() %>% 
    mutate(across(everything(), ~ pmax(1,ceiling(.x)))) %>% 
    mutate(dhdx = logs$dhdx$`T-SX`[x], dhdy = logs$dhdy$`T-SX`[y]) %>% 
    mutate(slope = sqrt(dhdx^2 + dhdy^2)) %>% 
    pull(slope)
})

# fits the local slope-dependent angle model
# location (mu_i) = negative gradient at base i (-\nabla_i)
# log concentration (log kappa_i) = beta0 + beta1*slope_i
fit_vm_local <- function(angles, slopes){
  # negative log-likelihood for local non-uniform angle model
  nll_vm_local <- function(params, angles, slopes){
    b0 <- params[1]
    b1 <- params[2]
    kappa <- exp(b0 + b1*slopes) # log link
    sapply(1:length(angles), function(i){
      dvonmises(circular(angles[i]),circular(0),circular(kappa[i]), log = T)
    }) %>% sum %>% {.*-1}
  }
  opt <- optim(c(0,0), nll_vm_local, angles = circular(angles), slopes = slopes)
  list(par = opt$par, loglik = -1*opt$value, aic = 2*opt$value + 4)
}

# fit all models
m_angle <- list()
m_angle$uniform <- with(logs, fit_vonmises(angle, mu = circular(0), kappa = circular(0)))
m_angle$vM <- with(logs, fit_vonmises(angle))
m_angle$vM_rel <- with(logs, fit_vonmises(circular(relAngle), mu = circular(0)))
m_angle$vM_rel_slope <- with(logs, fit_vm_local(circular(relAngle), slope_i))

# model comparison
m_angle %>% map("T-SX") %>% map_dbl("aic") %>% {. - min(.)} %>% round(3)
m_angle %>% map("W-FR") %>% map_dbl("aic") %>% {. - min(.)} %>% round(3)



