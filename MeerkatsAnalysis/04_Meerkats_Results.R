##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al. 2022, Ecology).
# 
# Eva Conquet, Arpat Ozgul, Daniel T. Blumstein, Kenneth B. Armitage,
# Madan K. Oli, Julien G. A. Martin, Tim H. Clutton-Brock, and Maria Paniw
#
# This script contains the code to obtain the result plots for the meerkat analysis. 
#
##########################################################################################################


###########################################################################
#
# 1. House keeping and loading libraries and data ----
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading libraries ----
# -----------------------

load.librairies = function(){
  library(ggplot2)
  library(boot)
}

load.librairies()


## 1.3. Loading data ----
# ------------------

load("MeerkatsProjectionsResults.RData")




###########################################################################
#
# 2. Analysing the results ----
#
###########################################################################

## 2.1. Building a data frame containing the metrics for each treatment ----
# ---------------------------------------------------------------------

sim.results = data.frame(treatment = c(rep("Control.Density.Dep", 500), rep("Control.Optimum.Dens", 500), rep("High.Seas.S.Surv.Density.Dep", 500), rep("High.Seas.S.Surv.Optimum.Dens", 500), rep("High.Seas.H.Surv.Density.Dep", 500), rep("High.Seas.H.Surv.Optimum.Dens", 500), rep("High.Seas.D.Surv.Density.Dep", 500), rep("High.Seas.D.Surv.Optimum.Dens", 500), rep("High.Seas.Emig.Density.Dep", 500), rep("High.Seas.Emig.Optimum.Dens", 500), rep("High.Seas.Trans.Density.Dep", 500), rep("High.Seas.Trans.Optimum.Dens", 500), rep("High.Seas.H.Recruit.Density.Dep", 500), rep("High.Seas.H.Recruit.Optimum.Dens", 500), rep("High.Seas.D.Recruit.Density.Dep", 500), rep("High.Seas.D.Recruit.Optimum.Dens", 500), rep("Low.Seas.S.Surv.Density.Dep", 500), rep("Low.Seas.S.Surv.Optimum.Dens", 500), rep("Low.Seas.H.Surv.Density.Dep", 500), rep("Low.Seas.H.Surv.Optimum.Dens", 500), rep("Low.Seas.D.Surv.Density.Dep", 500), rep("Low.Seas.D.Surv.Optimum.Dens", 500), rep("Low.Seas.Emig.Density.Dep", 500), rep("Low.Seas.Emig.Optimum.Dens", 500), rep("Low.Seas.Trans.Density.Dep", 500), rep("Low.Seas.Trans.Optimum.Dens", 500), rep("Low.Seas.H.Recruit.Density.Dep", 500), rep("Low.Seas.H.Recruit.Optimum.Dens", 500), rep("Low.Seas.D.Recruit.Density.Dep", 500), rep("Low.Seas.D.Recruit.Optimum.Dens", 500)),
                         scenario = c(rep("Control", 1000), rep("High seasonality", 7000), rep("Low seasonality", 7000)), 
                         density = c(rep(c(rep("Dens.Dep", 500), rep("Opt.Dens", 500)), 15)), 
                         vital.rate = c(rep("Control", 1000), rep(c(rep("S survival", 1000), rep("H survival", 1000), rep("D survival", 1000), rep("H emigration", 1000), rep("H-D transition", 1000), rep("H recruitment", 1000), rep("D recruitment", 1000)), 2)), 
                         stoch.lambda = c(control.dens.dep.sim.lambda, control.opt.dens.sim.lambda, high.seas.S.surv.dens.dep.sim.lambda, high.seas.S.surv.opt.dens.sim.lambda, high.seas.H.surv.dens.dep.sim.lambda, high.seas.H.surv.opt.dens.sim.lambda, high.seas.D.surv.dens.dep.sim.lambda, high.seas.D.surv.opt.dens.sim.lambda, high.seas.emig.dens.dep.sim.lambda, high.seas.emig.opt.dens.sim.lambda, high.seas.trans.dens.dep.sim.lambda, high.seas.trans.opt.dens.sim.lambda, high.seas.h.recruit.dens.dep.sim.lambda, high.seas.h.recruit.opt.dens.sim.lambda, high.seas.d.recruit.dens.dep.sim.lambda, high.seas.d.recruit.opt.dens.sim.lambda, low.seas.S.surv.dens.dep.sim.lambda, low.seas.S.surv.opt.dens.sim.lambda, low.seas.H.surv.dens.dep.sim.lambda, low.seas.H.surv.opt.dens.sim.lambda, low.seas.D.surv.dens.dep.sim.lambda, low.seas.D.surv.opt.dens.sim.lambda, low.seas.emig.dens.dep.sim.lambda, low.seas.emig.opt.dens.sim.lambda, low.seas.trans.dens.dep.sim.lambda, low.seas.trans.opt.dens.sim.lambda, low.seas.h.recruit.dens.dep.sim.lambda, low.seas.h.recruit.opt.dens.sim.lambda, low.seas.d.recruit.dens.dep.sim.lambda, low.seas.d.recruit.opt.dens.sim.lambda), 
                         var.lambda = c(control.dens.dep.var.lambda, control.opt.dens.var.lambda, high.seas.S.surv.dens.dep.var.lambda, high.seas.S.surv.opt.dens.var.lambda, high.seas.H.surv.dens.dep.var.lambda, high.seas.H.surv.opt.dens.var.lambda, high.seas.D.surv.dens.dep.var.lambda, high.seas.D.surv.opt.dens.var.lambda, high.seas.emig.dens.dep.var.lambda, high.seas.emig.opt.dens.var.lambda, high.seas.trans.dens.dep.var.lambda, high.seas.trans.opt.dens.var.lambda, high.seas.h.recruit.dens.dep.var.lambda, high.seas.h.recruit.opt.dens.var.lambda, high.seas.d.recruit.dens.dep.var.lambda, high.seas.d.recruit.opt.dens.var.lambda, low.seas.S.surv.dens.dep.var.lambda, low.seas.S.surv.opt.dens.var.lambda, low.seas.H.surv.dens.dep.var.lambda, low.seas.H.surv.opt.dens.var.lambda, low.seas.D.surv.dens.dep.var.lambda, low.seas.D.surv.opt.dens.var.lambda, low.seas.emig.dens.dep.var.lambda, low.seas.emig.opt.dens.var.lambda, low.seas.trans.dens.dep.var.lambda, low.seas.trans.opt.dens.var.lambda, low.seas.h.recruit.dens.dep.var.lambda, low.seas.h.recruit.opt.dens.var.lambda, low.seas.d.recruit.dens.dep.var.lambda, low.seas.d.recruit.opt.dens.var.lambda), 
                         ext.prob = c(control.sim.dens.dep[[2]], control.sim.opt.dens[[2]], high.seas.S.surv.dens.dep[[2]], high.seas.S.surv.opt.dens[[2]], high.seas.H.surv.dens.dep[[2]], high.seas.H.surv.opt.dens[[2]], high.seas.D.surv.dens.dep[[2]], high.seas.D.surv.opt.dens[[2]], high.seas.emig.dens.dep[[2]], high.seas.emig.opt.dens[[2]], high.seas.trans.dens.dep[[2]], high.seas.trans.opt.dens[[2]], high.seas.h.recruit.dens.dep[[2]], high.seas.h.recruit.opt.dens[[2]], high.seas.d.recruit.dens.dep[[2]], high.seas.d.recruit.opt.dens[[2]], low.seas.S.surv.dens.dep[[2]], low.seas.S.surv.opt.dens[[2]], low.seas.H.surv.dens.dep[[2]], low.seas.H.surv.opt.dens[[2]], low.seas.D.surv.dens.dep[[2]], low.seas.D.surv.opt.dens[[2]], low.seas.emig.dens.dep[[2]], low.seas.emig.opt.dens[[2]], low.seas.trans.dens.dep[[2]], low.seas.trans.opt.dens[[2]], low.seas.h.recruit.dens.dep[[2]], low.seas.h.recruit.opt.dens[[2]], low.seas.d.recruit.dens.dep[[2]], low.seas.d.recruit.opt.dens[[2]]))


## 2.2. Mean and quantiles of metrics for each scenario ----
# -----------------------------------------------------

# Stochastic lambda

aggregate(stoch.lambda ~ vital.rate + scenario + density, data = sim.results, FUN = "mean")
aggregate(stoch.lambda ~ vital.rate + scenario + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


# Variance

aggregate(var.lambda ~ vital.rate + scenario + density, data = sim.results, FUN = "mean")
aggregate(var.lambda ~ vital.rate + scenario + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


# Extinction probability

aggregate(ext.prob ~ vital.rate + scenario + density, data = sim.results, FUN = "mean")


###########################################################################
#
# 3. Plotting the results ----
#
###########################################################################

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim.results$vital.rate = factor(sim.results$vital.rate, 
                                 levels = c("Control", "S survival", "H survival", "D survival", "H emigration", "H-D transition", "H recruitment", "D recruitment")) # Changing the order of the factors for plotting

sim.results$treatment = factor(sim.results$treatment, 
                                levels = c("Control.Density.Dep", "High.Seas.S.Surv.Density.Dep", "Low.Seas.S.Surv.Density.Dep", "High.Seas.H.Surv.Density.Dep", "Low.Seas.H.Surv.Density.Dep", "High.Seas.D.Surv.Density.Dep", "Low.Seas.D.Surv.Density.Dep", "High.Seas.Emig.Density.Dep", "Low.Seas.Emig.Density.Dep", "High.Seas.Trans.Density.Dep", "Low.Seas.Trans.Density.Dep", "High.Seas.H.Recruit.Density.Dep", "Low.Seas.H.Recruit.Density.Dep", "High.Seas.D.Recruit.Density.Dep", "Low.Seas.D.Recruit.Density.Dep", "Control.Optimum.Dens", "High.Seas.S.Surv.Optimum.Dens", "Low.Seas.S.Surv.Optimum.Dens", "High.Seas.H.Surv.Optimum.Dens", "Low.Seas.H.Surv.Optimum.Dens", "High.Seas.D.Surv.Optimum.Dens", "Low.Seas.D.Surv.Optimum.Dens", "High.Seas.Emig.Optimum.Dens", "Low.Seas.Emig.Optimum.Dens", "High.Seas.Trans.Optimum.Dens", "Low.Seas.Trans.Optimum.Dens", "High.Seas.H.Recruit.Optimum.Dens", "Low.Seas.H.Recruit.Optimum.Dens", "High.Seas.D.Recruit.Optimum.Dens", "Low.Seas.D.Recruit.Optimum.Dens")) # Changing the order of the factors for plotting

facet_labels = c(Dens.Dep = "Density dependence", Opt.Dens = "Average density")

agg.meerkats.var = aggregate(var.lambda ~ treatment, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))
agg.meerkats.lambda = aggregate(stoch.lambda ~ treatment, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


# Stochastic lambda

png(filename = "StochLambda.png",
    width = 5000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.stoch.lambda = ggplot(sim.results,  aes(x = vital.rate,  y = stoch.lambda,  fill = scenario)) +
 facet_wrap(~ density, scale = "free", labeller = labeller(density = facet_labels)) +
 geom_boxplot(ymin = agg.meerkats.lambda$stoch.lambda[, 1], ymax = agg.meerkats.lambda$stoch.lambda[, 3]) +
 stat_summary(fun = mean,  colour = "darkred",  geom = "point",  
        shape = 17,  size = 3, show.legend = FALSE, position = position_dodge(width = 0.75)) +
 labs(title = "Stochastic log λ of the meerkat population\ndepending on the simulation scenario",  
      x = "Perturbed vital rate",  
      y = "Stochastic log growth rate") +
 scale_fill_manual(values = cbbPalette, 
          name = "Scenario") +
 theme_bw() +
 theme(axis.title.x = element_text(size = 40, colour = "black",  margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"),  
       axis.title.y = element_text(size = 40, colour = "black",  margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"),  
       axis.text.x = element_text(size = 25, colour = "black",  margin = margin(t = 10, r = 0, b = 0, l = 0),  angle = 45,  hjust = 1),  
       axis.text.y = element_text(size = 25, colour = "black",  margin = margin(t = 0, r = 10, b = 0, l = 0)),  
       plot.title = element_text(size = 50,  hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"),  
       legend.text = element_text(size = 20),  
       legend.title = element_text(face = "bold", size = 30),  
       legend.position = "right",  
       legend.key.size = unit(6, "lines"), 
       strip.text.x = element_text(size = 20, face = "bold"))

plot.stoch.lambda

dev.off()


# Var(lambda)

tiff(filename = "VarLambda.tiff", 
  width = 8, 
  height = 3.5, 
  units = "in", 
  bg = "white", 
  res = 600, 
  compression = "lzw")

plot.var.lambda = ggplot(sim.results, aes(x = vital.rate, y = var.lambda, fill = scenario)) +
 facet_wrap(~ density, labeller = labeller(density = facet_labels)) +
 geom_boxplot(ymin = agg.meerkats.var$var.lambda[, 1], ymax = agg.meerkats.var$var.lambda[, 3], outlier.size = 0.1, size = 0.1) +
 stat_summary(fun = mean, colour = "darkred", geom = "point", 
              shape = 17, size = 1, show.legend = TRUE, position = position_dodge(width = 0.75)) +
 labs(title = "", x = "Perturbed vital rate", y = "Variance of annual \u03BB") +
 scale_fill_manual(values = cbbPalette, 
          name = "Scenario") +
 theme_bw() +
 theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
       axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
       axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
       axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
       legend.text = element_text(size = 7), 
       legend.title = element_text(size = 9), 
       legend.position = "right", 
       legend.key.size = unit(2, "lines"), 
       strip.text.x = element_text(size = 7),
       strip.background = element_blank())

plot.var.lambda

dev.off()




###########################################################################
#
# 4. Saving files and data ----
#
###########################################################################

write.csv(sim.results, "SimResults_Meerkats.csv", row.names = F)
