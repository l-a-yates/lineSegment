
library(spatstat)
library(fitdistrplus)
library(circular)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggplotify)


rm(list=(ls()))
select <-  dplyr::select

# fallen-log data
logs <- readRDS("output/logs_2021_07_06.rds")
sites <- logs %>% rownames
names(sites) <- sites

# spatstat psp objects
owin.site <- logs$logs.psp[[1]]$window
psp.site <- logs$logs.psp

#-----------------
# plot data sets
#-----------------

plotLogs <- function(site_psp, title = NULL, subtitle = NULL){
  site_psp %>% 
    as.data.frame %>% 
    ggplot(aes(x = x0, y = y0)) +
    geom_segment(aes(xend = x1, yend = y1), alpha = 1, col = blues9[8], size = 0.5) +
    geom_point(size = 1, col = "grey10") +
    coord_fixed(ratio = 1) +
    theme_bw() +
    labs(x = " ", y = " ", title = title, subtitle = subtitle) +
    theme(panel.background = element_rect(fill = "grey99"), 
          panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()) 
    #theme(panel.background = element_rect(fill = "grey95"), panel.grid = element_blank())
}

tt1 <- theme(plot.title = element_text(colour = "grey10", size = 12, face = "bold", hjust = 0.5))

ggarrange(plotLogs(psp.site$`T-SX`, "T-SX") + tt1,
          plotLogs(psp.site$`W-FR`, "W-FR") + tt1)

#ggsave("plots/fallen_log_plots_2021_07_06.pdf", dpi = 300, height = 4, width = 8, device = cairo_pdf)
#ggsave("fallen_log_plots.eps", dpi = 300, height = 4, width = 8, device = cairo_ps)


#-----------------------------------
# Extra plots for model explanations
#-----------------------------------

set.seed(2826) # 2854
#base_inhom <- rpoispp(function(x,y){x/100^2}, win = owin.site)
base_inhom <- rpoispp(function(x,y){exp(-(x)/50)/100}, win = owin.site)
base_hom <- rpoispp(intensity(base_inhom), win = owin.site)

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


#ggsave("demo_logs_v1.pdf", width = 6, height = 5)
#ggsave("plots/demo_logs_v1d.pdf", width = 6, height = 5)