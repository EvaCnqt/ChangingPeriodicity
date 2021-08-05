##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script contains the code to obtain the result plots for the dewy pine analysis. 
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
  library(boot)
  library(cowplot)
  library(ggplot2)
}

load.librairies()


## 1.3. Loading data ----
# ------------------

load("DewyPineProjectionsResults.RData")




###########################################################################
#
# 2. Analysing the results ----
#
###########################################################################

## 2.1. Building a data frame for the extinction probability, stochastic lambda, and var(lambda) ----
# ----------------------------------------------------------------------------------------------

sim.results = data.frame(treatment = c(rep(c(rep("Control", 500), rep("TSF0", 500), rep("TSF1", 500), rep("TSF2", 500), rep("TSF3", 500), rep("TSF>3 - 0.1", 500), rep("TSF>3 - 0.3", 500), rep("TSF>3 - 0.5", 500), rep("TSF>3 - 0.7", 500), rep("TSF>3 - 0.9", 500), rep("TSF>3 - 1", 500)), 2)), 
                         density = c(rep("DensDep", 5500), rep("AvgAb", 5500)), 
                         stoch.lambda = c(control.sim.lambda, perturb.state1.sim.lambda, perturb.state2.sim.lambda, perturb.state3.sim.lambda, perturb.state4.sim.lambda, perturb.state5.proportion.10.sim.lambda, perturb.state5.proportion.30.sim.lambda, perturb.state5.proportion.50.sim.lambda, perturb.state5.proportion.70.sim.lambda, perturb.state5.proportion.90.sim.lambda, perturb.state5.sim.lambda, control.avgab.sim.lambda, perturb.state1.avgab.sim.lambda, perturb.state2.avgab.sim.lambda, perturb.state3.avgab.sim.lambda, perturb.state4.avgab.sim.lambda, perturb.state5.avgab.proportion.10.sim.lambda, perturb.state5.avgab.proportion.30.sim.lambda, perturb.state5.avgab.proportion.50.sim.lambda, perturb.state5.avgab.proportion.70.sim.lambda, perturb.state5.avgab.proportion.90.sim.lambda, perturb.state5.avgab.sim.lambda), 
                         var.lambda = c(control.var.lambda, perturb.state1.var.lambda, perturb.state2.var.lambda, perturb.state3.var.lambda, perturb.state4.var.lambda, perturb.state5.proportion.10.var.lambda, perturb.state5.proportion.30.var.lambda, perturb.state5.proportion.50.var.lambda, perturb.state5.proportion.70.var.lambda, perturb.state5.proportion.90.var.lambda, perturb.state5.var.lambda, control.avgab.var.lambda, perturb.state1.avgab.var.lambda, perturb.state2.avgab.var.lambda, perturb.state3.avgab.var.lambda, perturb.state4.avgab.var.lambda, perturb.state5.avgab.proportion.10.var.lambda, perturb.state5.avgab.proportion.30.var.lambda, perturb.state5.avgab.proportion.50.var.lambda, perturb.state5.avgab.proportion.70.var.lambda, perturb.state5.avgab.proportion.90.var.lambda, perturb.state5.avgab.var.lambda), 
                         ext.prob = c(control.sim[[3]], perturb.state1[[3]], perturb.state2[[3]], perturb.state3[[3]], perturb.state4[[3]], perturb.state5.proportion.10[[3]], perturb.state5.proportion.30[[3]], perturb.state5.proportion.50[[3]], perturb.state5.proportion.70[[3]], perturb.state5.proportion.90[[3]], perturb.state5[[3]], control.sim.avgab[[3]], perturb.state1.avgab[[3]], perturb.state2.avgab[[3]], perturb.state3.avgab[[3]], perturb.state4.avgab[[3]], perturb.state5.avgab.proportion.10[[3]], perturb.state5.avgab.proportion.30[[3]], perturb.state5.avgab.proportion.50[[3]], perturb.state5.avgab.proportion.70[[3]], perturb.state5.avgab.proportion.90[[3]], perturb.state5.avgab[[3]]))


# Reordering factor levels for plotting
sim.results$treatment = factor(sim.results$treatment, levels = c("Control", "TSF0", "TSF1", "TSF2", "TSF3", "TSF>3 - 0.1", "TSF>3 - 0.3", "TSF>3 - 0.5", "TSF>3 - 0.7", "TSF>3 - 0.9", "TSF>3 - 1"))


## 2.2. Magnitude of changes between scenarios (mean and quantiles) ----
# -----------------------------------------------------------------

# Stochastic lambda
aggregate(stoch.lambda ~ treatment + density, data = sim.results, FUN = "mean")
aggregate(stoch.lambda ~ treatment + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))

# Variance
aggregate(var.lambda ~ treatment + density, data = sim.results, FUN = "mean")
aggregate(var.lambda ~ treatment + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))

aggregate(ext.prob ~ treatment + density, data = sim.results, FUN = "mean")




###########################################################################
#
# 3. Plotting the results ----
#
###########################################################################

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim.results$density <- factor(sim.results$density, levels = c("DensDep", "AvgAb"))

facet_labels <- c(DensDep = "Density dependence", AvgAb = "Average abundance")

agg.stoch.lambda = aggregate(stoch.lambda~treatment+density,data=sim.results,FUN=function(x) quantile(x,c(0.025,0.5,0.975)))
agg.var.lambda = aggregate(var.lambda~treatment+density,data=sim.results,FUN=function(x) quantile(x,c(0.025,0.5,0.975)))

agg.ext.prob = aggregate(ext.prob~treatment+density,data=sim.results,FUN=mean)
agg.ext.prob.quantiles = aggregate(ext.prob~treatment+density,data=sim.results,FUN=function(x) quantile(x,c(0.025,0.5,0.975)))


## 3.1. Stochastic lambda ----
# -----------------------

plot.stoch.lambda = ggplot(sim.results, aes(x = treatment, y = stoch.lambda)) +
  facet_wrap( ~ density, labeller = labeller(density = facet_labels)) +
  geom_boxplot(ymin = agg.stoch.lambda$stoch.lambda[, 1], ymax = agg.stoch.lambda$stoch.lambda[, 3]) +
  stat_summary(fun = mean, colour = "darkred", aes(shape = "Mean variance"), geom = "point", size = 6) +
  labs(x = "Perturbed state", y = "Variance of annual \u03BB") +
  scale_shape_manual("", values = c("Mean variance" = 17)) +
  scale_x_discrete(labels = c("Control", expression("TSF"["0"]), expression("TSF"["1"]), expression("TSF"["2"]), expression("TSF"["3"]), expression("TSF"[">3"]~"(0.1)"), expression("TSF"[">3"]~"(0.3)"), expression("TSF"[">3"]~"(0.5)"), expression("TSF"[">3"]~"(0.7)"), expression("TSF"[">3"]~"(0.9)"), expression("TSF"[">3"]~"(all)"))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 6, r =0 , b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 6, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(2, "lines"), 
        strip.text.x = element_text(size = 20, face = "bold"))


## 3.2. Combined plot (Fig. S10.3) ----
# --------------------------------

plot.legend.appendix = ggplot(agg.ext.prob, aes(x = treatment, y = ext.prob, colour = density)) +
  geom_point(size = 6, position = position_dodge(1)) +
  labs(x = "Perturbed state", y = "Extinction probability") +
  scale_colour_manual(values = c("#000000", "#009E73"),
                      labels = c("Average\ndensity", "Density\ndependence"),
                      name = "Scenario") +
  scale_x_discrete(labels = c("Control", expression("TSF"["0"]), expression("TSF"["1"]), expression("TSF"["2"]), expression("TSF"["3"]), expression("TSF"[">3"]~"(0.1)"), expression("TSF"[">3"]~"(0.3)"), expression("TSF"[">3"]~"(0.5)"), expression("TSF"[">3"]~"(0.7)"), expression("TSF"[">3"]~"(0.9)"), expression("TSF"[">3"]~"(all)"))) +
  ylim(0, 1) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 6, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 6, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "right", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(2, "lines"))


plot.var.lambda.appendix = ggplot(sim.results, aes(x = treatment, y = var.lambda)) +
  facet_wrap(~ density,labeller = labeller(density = facet_labels)) +
  geom_boxplot(ymin = agg.var.lambda$var.lambda[, 1], ymax = agg.var.lambda$var.lambda[, 3]) +
  stat_summary(fun = mean, colour = "darkred", aes(shape = "Mean variance"), geom = "point", size = 6) +
  labs(x = "Perturbed state", y = "Variance of annual \u03BB") +
  scale_shape_manual("", values = c("Mean variance" = 17)) +
  scale_x_discrete(labels = c("Control", expression("TSF"["0"]), expression("TSF"["1"]), expression("TSF"["2"]), expression("TSF"["3"]), expression("TSF"[">3"]~"(0.1)"), expression("TSF"[">3"]~"(0.3)"), expression("TSF"[">3"]~"(0.5)"), expression("TSF"[">3"]~"(0.7)"), expression("TSF"[">3"]~"(0.9)"), expression("TSF"[">3"]~"(all)"))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 6, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 6, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(2, "lines"), 
        strip.text.x = element_text(size = 20, face = "bold"))


plot.ext.prob.appendix = ggplot(agg.ext.prob, aes(x = treatment, y = ext.prob, colour = density)) +
  geom_point(size = 6, position = position_dodge(1)) +
  labs(x = "Perturbed state", y = "Extinction probability") +
  scale_colour_manual(values = c("#000000", "#009E73"),
                      labels = c("Average density", "Density dependence"),
                      name = "Scenario") +
  scale_x_discrete(labels = c("Control", expression("TSF"["0"]), expression("TSF"["1"]), expression("TSF"["2"]), expression("TSF"["3"]), expression("TSF"[">3"]~"(0.1)"), expression("TSF"[">3"]~"(0.3)"), expression("TSF"[">3"]~"(0.5)"), expression("TSF"[">3"]~"(0.7)"), expression("TSF"[">3"]~"(0.9)"), expression("TSF"[">3"]~"(all)"))) +
  ylim(0, 1) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 6, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 6, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none", legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(2, "lines"))


dp.appendix = ggdraw() +
  draw_plot(plot.var.lambda.appendix, x = 0, y = -0.01, width = .55, height = 1) +
  draw_plot(plot.ext.prob.appendix, x = .55, y = -0.01, width = .33, height = 1) +
  draw_plot(get_legend(plot.legend.appendix), x = 0.663, y = -0.01, width = .55, height = 1.2) +
  draw_plot_label(label = c("(a)", "(b)"), size = 26,
                  x = c(0, 0.537), y = c(1, 1))


png(filename = "DPAppendix10.png",
    width=6000,
    height=2000,
    units="px",
    bg="white",
    res=300,
    type = "cairo")

dp.appendix

dev.off()




###########################################################################
#
# 4. Saving files ----
#
###########################################################################

write.csv(sim.results,"SimResults_DP.csv",row.names = F)
write.csv(add.sim.results,"SimResults_Additional_DP.csv",row.names = F)
