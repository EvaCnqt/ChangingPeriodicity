##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script contains the code to obtain the result plots for the meerkat analysis. 
#
# Author: Eva Conquet
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

load("MeerkatsProjectionsResults_SeasonOnly.RData")




###########################################################################
#
# 2. Analysing the results ----
#
###########################################################################

## 2.1. Building a data frame containing the metrics for each treatment ----
# ---------------------------------------------------------------------

sim.results = data.frame(treatment = c(rep("Control", 500), rep("High.Seas.S.Surv", 500), rep("High.Seas.H.Surv", 500), rep("High.Seas.D.Surv", 500), rep("High.Seas.Emig", 500), rep("High.Seas.Trans", 500), rep("High.Seas.H.Recruit", 500), rep("High.Seas.D.Recruit", 500), rep("Low.Seas.S.Surv", 500), rep("Low.Seas.H.Surv", 500), rep("Low.Seas.D.Surv", 500), rep("Low.Seas.Emig", 500), rep("Low.Seas.Trans", 500), rep("Low.Seas.H.Recruit", 500), rep("Low.Seas.D.Recruit", 500)),
                         scenario = c(rep("Control", 500), rep("High seasonality", 3500), rep("Low seasonality", 3500)), 
                         vital.rate = c(rep("Control", 500), rep(c(rep("S survival", 500), rep("H survival", 500), rep("D survival", 500), rep("H emigration", 500), rep("H-D transition", 500), rep("H recruitment", 500), rep("D recruitment", 500)), 2)), 
                         stoch.lambda = c(control.sim.lambda, high.seas.S.surv.sim.lambda, high.seas.H.surv.sim.lambda, high.seas.D.surv.sim.lambda, high.seas.emig.sim.lambda, high.seas.trans.sim.lambda, high.seas.h.recruit.sim.lambda, high.seas.d.recruit.sim.lambda, low.seas.S.surv.sim.lambda, low.seas.H.surv.sim.lambda, low.seas.D.surv.sim.lambda, low.seas.emig.sim.lambda, low.seas.trans.sim.lambda, low.seas.h.recruit.sim.lambda, low.seas.d.recruit.sim.lambda), 
                         var.lambda = c(control.var.lambda, high.seas.S.surv.var.lambda, high.seas.H.surv.var.lambda, high.seas.D.surv.var.lambda, high.seas.emig.var.lambda, high.seas.trans.var.lambda, high.seas.h.recruit.var.lambda, high.seas.d.recruit.var.lambda, low.seas.S.surv.var.lambda, low.seas.H.surv.var.lambda, low.seas.D.surv.var.lambda, low.seas.emig.var.lambda, low.seas.trans.var.lambda, low.seas.h.recruit.var.lambda, low.seas.d.recruit.var.lambda), 
                         ext.prob = c(control.sim[[2]], high.seas.S.surv[[2]], high.seas.H.surv[[2]], high.seas.D.surv[[2]], high.seas.emig[[2]], high.seas.trans[[2]], high.seas.h.recruit[[2]], high.seas.d.recruit[[2]], low.seas.S.surv[[2]], low.seas.H.surv[[2]], low.seas.D.surv[[2]], low.seas.emig[[2]], low.seas.trans[[2]], low.seas.h.recruit[[2]], low.seas.d.recruit[[2]]))


## 2.2. Mean and quantiles of metrics for each scenario ----
# -----------------------------------------------------

# Stochastic lambda

aggregate(stoch.lambda ~ vital.rate + scenario, data = sim.results, FUN = "mean")
aggregate(stoch.lambda ~ vital.rate + scenario, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


# Variance

aggregate(var.lambda ~ vital.rate + scenario, data = sim.results, FUN = "mean")
aggregate(var.lambda ~ vital.rate + scenario, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


# Extinction probability

aggregate(ext.prob ~ vital.rate + scenario, data = sim.results, FUN = "mean")


###########################################################################
#
# 3. Plotting the results ----
#
###########################################################################

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim.results$vital.rate = factor(sim.results$vital.rate, 
                                 levels = c("Control", "S survival", "H survival", "D survival", "H emigration", "H-D transition", "H recruitment", "D recruitment")) # Changing the order of the factors for plotting

sim.results$treatment = factor(sim.results$treatment, 
                                levels = c("Control", "High.Seas.S.Surv", "Low.Seas.S.Surv", "High.Seas.H.Surv", "Low.Seas.H.Surv", "High.Seas.D.Surv", "Low.Seas.D.Surv", "High.Seas.Emig", "Low.Seas.Emig", "High.Seas.Trans", "Low.Seas.Trans", "High.Seas.H.Recruit", "Low.Seas.H.Recruit", "High.Seas.D.Recruit", "Low.Seas.D.Recruit")) # Changing the order of the factors for plotting

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
       strip.text.x = element_text(size = 7))

plot.var.lambda

dev.off()




###########################################################################
#
# 4. Saving files and data
#
###########################################################################

write.csv(sim.results, "SimResults_Meerkats_SeasonOnly.csv", row.names = F)
