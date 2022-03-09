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
  library(cowplot)
}

load.librairies()


## 1.3. Loading data ----
# ------------------

load("MeerkatsProjectionsResults_DensityOnly.RData")

res.meerkats.seasonOnly = read.csv("SimResults_Meerkats_SeasonOnly.csv")




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
                         stoch.lambda = c(control.dens.dep.sim.lambda, high.seas.S.surv.dens.dep.sim.lambda, high.seas.H.surv.dens.dep.sim.lambda, high.seas.D.surv.dens.dep.sim.lambda, high.seas.emig.dens.dep.sim.lambda, high.seas.trans.dens.dep.sim.lambda, high.seas.h.recruit.dens.dep.sim.lambda, high.seas.d.recruit.dens.dep.sim.lambda, low.seas.S.surv.dens.dep.sim.lambda, low.seas.H.surv.dens.dep.sim.lambda, low.seas.D.surv.dens.dep.sim.lambda, low.seas.emig.dens.dep.sim.lambda, low.seas.trans.dens.dep.sim.lambda, low.seas.h.recruit.dens.dep.sim.lambda, low.seas.d.recruit.dens.dep.sim.lambda), 
                         var.lambda = c(control.dens.dep.var.lambda, high.seas.S.surv.dens.dep.var.lambda, high.seas.H.surv.dens.dep.var.lambda, high.seas.D.surv.dens.dep.var.lambda, high.seas.emig.dens.dep.var.lambda, high.seas.trans.dens.dep.var.lambda, high.seas.h.recruit.dens.dep.var.lambda, high.seas.d.recruit.dens.dep.var.lambda, low.seas.S.surv.dens.dep.var.lambda, low.seas.H.surv.dens.dep.var.lambda, low.seas.D.surv.dens.dep.var.lambda, low.seas.emig.dens.dep.var.lambda, low.seas.trans.dens.dep.var.lambda, low.seas.h.recruit.dens.dep.var.lambda, low.seas.d.recruit.dens.dep.var.lambda), 
                         ext.prob = c(control.sim.dens.dep[[2]], high.seas.S.surv.sim.dens.dep[[2]], high.seas.H.surv.sim.dens.dep[[2]], high.seas.D.surv.sim.dens.dep[[2]], high.seas.emig.sim.dens.dep[[2]], high.seas.trans.sim.dens.dep[[2]], high.seas.h.recruit.sim.dens.dep[[2]], high.seas.d.recruit.sim.dens.dep[[2]], low.seas.S.surv.sim.dens.dep[[2]], low.seas.H.surv.sim.dens.dep[[2]], low.seas.D.surv.sim.dens.dep[[2]], low.seas.emig.sim.dens.dep[[2]], low.seas.trans.sim.dens.dep[[2]], low.seas.h.recruit.sim.dens.dep[[2]], low.seas.d.recruit.sim.dens.dep[[2]]))


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

plot.stoch.lambda.densityOnly = ggplot(sim.results,  aes(x = vital.rate,  y = stoch.lambda,  fill = scenario)) +
  geom_boxplot(ymin = agg.meerkats.lambda$stoch.lambda[, 1], ymax = agg.meerkats.lambda$stoch.lambda[, 3], outlier.size = 0.1, size = 0.1) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", 
               shape = 17, size = 1, show.legend = TRUE, position = position_dodge(width = 0.75)) +
  labs(title = "No season",  
       x = "Perturbed vital rate",  
       y = "Stochastic log \u03BB") +
  scale_fill_manual(values = cbbPalette, 
                    name = "Scenario") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", 
        legend.key.size = unit(2, "lines"), 
        strip.text.x = element_text(size = 7),
        plot.title = element_text(size = 9, hjust = 0.5))

plot.stoch.lambda.densityOnly

dev.off()


# Var(lambda)

tiff(filename = "VarLambda.tiff", 
     width = 8, 
     height = 3.5, 
     units = "in", 
     bg = "white", 
     res = 600, 
     compression = "lzw")

plot.var.lambda.densityOnly = ggplot(sim.results, aes(x = vital.rate, y = var.lambda, fill = scenario)) +
  geom_boxplot(ymin = agg.meerkats.var$var.lambda[, 1], ymax = agg.meerkats.var$var.lambda[, 3], outlier.size = 0.1, size = 0.1) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", 
               shape = 17, size = 1, show.legend = TRUE, position = position_dodge(width = 0.75)) +
  labs(title = "No season", x = "Perturbed vital rate", y = "Variance of annual \u03BB") +
  scale_fill_manual(values = cbbPalette, 
                    name = "Scenario") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", 
        legend.key.size = unit(2, "lines"), 
        strip.text.x = element_text(size = 7),
        plot.title = element_text(size = 9, hjust = 0.5))

plot.var.lambda.densityOnly

dev.off()



res.meerkats.seasonOnly$vital.rate = factor(res.meerkats.seasonOnly$vital.rate, 
                                levels = c("Control", "S survival", "H survival", "D survival", "H emigration", "H-D transition", "H recruitment", "D recruitment")) # Changing the order of the factors for plotting

res.meerkats.seasonOnly$treatment = factor(res.meerkats.seasonOnly$treatment, 
                               levels = c("Control", "High.Seas.S.Surv", "Low.Seas.S.Surv", "High.Seas.H.Surv", "Low.Seas.H.Surv", "High.Seas.D.Surv", "Low.Seas.D.Surv", "High.Seas.Emig", "Low.Seas.Emig", "High.Seas.Trans", "Low.Seas.Trans", "High.Seas.H.Recruit", "Low.Seas.H.Recruit", "High.Seas.D.Recruit", "Low.Seas.D.Recruit")) # Changing the order of the factors for plotting

agg.meerkats.var.seasonOnly = aggregate(var.lambda ~ treatment, data = res.meerkats.seasonOnly, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))
agg.meerkats.lambda.seasonOnly = aggregate(stoch.lambda ~ treatment, data = res.meerkats.seasonOnly, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


# Stochastic lambda

png(filename = "StochLambda.png",
    width = 5000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.stoch.lambda.seasonOnly = ggplot(res.meerkats.seasonOnly,  aes(x = vital.rate,  y = stoch.lambda,  fill = scenario)) +
  geom_boxplot(ymin = agg.meerkats.lambda.seasonOnly$stoch.lambda[, 1], ymax = agg.meerkats.lambda.seasonOnly$stoch.lambda[, 3], outlier.size = 0.1, size = 0.1) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", 
               shape = 17, size = 1, show.legend = TRUE, position = position_dodge(width = 0.75)) +
  labs(title = "No density dependence",  
       x = "Perturbed vital rate",  
       y = "Stochastic log \u03BB") +
  scale_fill_manual(values = cbbPalette, 
                    name = "Scenario") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", 
        legend.key.size = unit(2, "lines"), 
        strip.text.x = element_text(size = 7),
        plot.title = element_text(size = 9, hjust = 0.5))

plot.stoch.lambda.seasonOnly

dev.off()


# Var(lambda)

tiff(filename = "VarLambda.tiff", 
     width = 8, 
     height = 3.5, 
     units = "in", 
     bg = "white", 
     res = 600, 
     compression = "lzw")

plot.var.lambda.seasonOnly = ggplot(res.meerkats.seasonOnly, aes(x = vital.rate, y = var.lambda, fill = scenario)) +
  geom_boxplot(ymin = agg.meerkats.var.seasonOnly$var.lambda[, 1], ymax = agg.meerkats.var.seasonOnly$var.lambda[, 3], outlier.size = 0.1, size = 0.1) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", 
               shape = 17, size = 1, show.legend = TRUE, position = position_dodge(width = 0.75)) +
  labs(title = "No density dependence", x = "Perturbed vital rate", y = "Variance of annual \u03BB") +
  scale_fill_manual(values = cbbPalette, 
                    name = "Scenario") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", 
        legend.key.size = unit(2, "lines"), 
        strip.text.x = element_text(size = 7),
        plot.title = element_text(size = 9, hjust = 0.5))

plot.var.lambda.seasonOnly

dev.off()



plot.legend = ggplot(res.meerkats.seasonOnly, aes(x = vital.rate, y = var.lambda, fill = scenario)) +
  geom_boxplot(ymin = agg.meerkats.var.seasonOnly$var.lambda[, 1], ymax = agg.meerkats.var.seasonOnly$var.lambda[, 3], outlier.size = 0.1, size = 0.1) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", 
               shape = 17, size = 1, show.legend = TRUE, position = position_dodge(width = 0.75)) +
  labs(title = "", x = "Perturbed vital rate", y = "Variance of annual \u03BB") +
  scale_fill_manual(values = cbbPalette, 
                    name = "") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "bottom", 
        legend.key.size = unit(2, "lines"), 
        strip.text.x = element_text(size = 7),
        legend.background = element_rect(colour = "transparent", fill = "transparent"))


results.fig = cowplot::ggdraw() +
  draw_plot(plot.stoch.lambda.densityOnly, x = 0, y = 0.525, width = 0.5, height = 0.48) +
  draw_plot(plot.var.lambda.densityOnly, x = .5, y = 0.525, width = .5, height = .48) +
  draw_plot(plot.stoch.lambda.seasonOnly, x = 0, y = 0, width=.5, height=0.48) +
  draw_plot(plot.var.lambda.seasonOnly, x = 0.5, y = 0, width = .5, height = 0.48) +
  draw_plot(get_legend(plot.legend), x = 0.45, y = 0.4, width = 0.1, height = 0.2) +
  draw_plot_label(label = c("(a)", "(b)", "(c)", "(d)"), 
                  size = 10, 
                  fontface = "plain",
                  x = c(0.02, 0.52, 0.02, 0.52), y = c(0.995, 0.995, 0.47, 0.47))


tiff(filename = "Meerkats_SeasonAndDensityOnly.tiff",
     width=6,
     height=6,
     units="in",
     bg="white",
     res=800,
     compression = "lzw")

results.fig

dev.off()

###########################################################################
#
# 4. Saving files and data
#
###########################################################################

write.csv(sim.results, "SimResults_Meerkats_SeasonOnly.csv", row.names = F)
write.csv(agg.meerkats.lambda, "Meerkats_DensityOnly_StochLambda_MedianQuantiles.csv", row.names = F)
write.csv(agg.meerkats.var, "Meerkats_DensityOnly_VarLambda_MedianQuantiles.csv", row.names = F)
write.csv(agg.meerkats.lambda.seasonOnly, "Meerkats_SeasonOnly_StochLambda_MedianQuantiles.csv", row.names = F)
write.csv(agg.meerkats.var.seasonOnly, "Meerkats_SeasonOnly_VarLambda_MedianQuantiles.csv", row.names = F)
