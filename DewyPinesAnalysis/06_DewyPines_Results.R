##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
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

load("DewyPine_Projections_Results.RData")




###########################################################################
#
# 2. Analysing the results ----
#
###########################################################################

## 2.1. Building a data frame for the extinction probability, stochastic lambda, and var(lambda) ----
# ----------------------------------------------------------------------------------------------

sim.results = data.frame(treatment = c(rep(c(rep("Control", 500), rep("TSF>3", 500), rep("TSF3 and >3", 500), rep("TSF2, 3, and >3", 500), rep("TSF1, 2, 3, and >3", 500), rep("All TSFs", 500)), 8)), 
                         disturbance = c(rep("Stochastic fire (p = 0.066)", 6000), rep("Stochastic fire (p = 0.033)", 6000), rep("Periodic fire (f = 15 years)", 6000), rep("Periodic fire (f = 30 years)", 6000)),
                         density = c(rep(c(rep("DensDep", 3000), rep("AvgDens", 3000)), 4)), 
                         stoch.lambda = c(stochasticFire15.control.sim.lambda, stochasticFire15.perturb.state5.sim.lambda, stochasticFire15.perturb.state4_5.sim.lambda, stochasticFire15.perturb.state3_4_5.sim.lambda, stochasticFire15.perturb.state2_3_4_5.sim.lambda, stochasticFire15.perturb.stateAll.sim.lambda,
                                          stochasticFire15.control.avgdens.sim.lambda, stochasticFire15.perturb.state5.avgdens.sim.lambda, stochasticFire15.perturb.state4_5.avgdens.sim.lambda, stochasticFire15.perturb.state3_4_5.avgdens.sim.lambda, stochasticFire15.perturb.state2_3_4_5.avgdens.sim.lambda, stochasticFire15.perturb.stateAll.avgdens.sim.lambda,
                                          stochasticFire30.control.sim.lambda, stochasticFire30.perturb.state5.sim.lambda, stochasticFire30.perturb.state4_5.sim.lambda, stochasticFire30.perturb.state3_4_5.sim.lambda, stochasticFire30.perturb.state2_3_4_5.sim.lambda, stochasticFire30.perturb.stateAll.sim.lambda,
                                          stochasticFire30.control.avgdens.sim.lambda, stochasticFire30.perturb.state5.avgdens.sim.lambda, stochasticFire30.perturb.state4_5.avgdens.sim.lambda, stochasticFire30.perturb.state3_4_5.avgdens.sim.lambda, stochasticFire30.perturb.state2_3_4_5.avgdens.sim.lambda, stochasticFire30.perturb.stateAll.avgdens.sim.lambda,
                                          periodicFire15.control.sim.lambda, periodicFire15.perturb.state5.sim.lambda, periodicFire15.perturb.state4_5.sim.lambda, periodicFire15.perturb.state3_4_5.sim.lambda, periodicFire15.perturb.state2_3_4_5.sim.lambda, periodicFire15.perturb.stateAll.sim.lambda,
                                          periodicFire15.control.avgdens.sim.lambda, periodicFire15.perturb.state5.avgdens.sim.lambda, periodicFire15.perturb.state4_5.avgdens.sim.lambda, periodicFire15.perturb.state3_4_5.avgdens.sim.lambda, periodicFire15.perturb.state2_3_4_5.avgdens.sim.lambda, periodicFire15.perturb.stateAll.avgdens.sim.lambda,
                                          periodicFire30.control.sim.lambda, periodicFire30.perturb.state5.sim.lambda, periodicFire30.perturb.state4_5.sim.lambda, periodicFire30.perturb.state3_4_5.sim.lambda, periodicFire30.perturb.state2_3_4_5.sim.lambda, periodicFire30.perturb.stateAll.sim.lambda,
                                          periodicFire30.control.avgdens.sim.lambda, periodicFire30.perturb.state5.avgdens.sim.lambda, periodicFire30.perturb.state4_5.avgdens.sim.lambda, periodicFire30.perturb.state3_4_5.avgdens.sim.lambda, periodicFire30.perturb.state2_3_4_5.avgdens.sim.lambda, periodicFire30.perturb.stateAll.avgdens.sim.lambda), 
                         var.lambda = c(stochasticFire15.control.var.lambda, stochasticFire15.perturb.state5.var.lambda, stochasticFire15.perturb.state4_5.var.lambda, stochasticFire15.perturb.state3_4_5.var.lambda, stochasticFire15.perturb.state2_3_4_5.var.lambda, stochasticFire15.perturb.stateAll.var.lambda,
                                        stochasticFire15.control.avgdens.var.lambda, stochasticFire15.perturb.state5.avgdens.var.lambda, stochasticFire15.perturb.state4_5.avgdens.var.lambda, stochasticFire15.perturb.state3_4_5.avgdens.var.lambda, stochasticFire15.perturb.state2_3_4_5.avgdens.var.lambda, stochasticFire15.perturb.stateAll.avgdens.var.lambda,
                                        stochasticFire30.control.var.lambda, stochasticFire30.perturb.state5.var.lambda, stochasticFire30.perturb.state4_5.var.lambda, stochasticFire30.perturb.state3_4_5.var.lambda, stochasticFire30.perturb.state2_3_4_5.var.lambda, stochasticFire30.perturb.stateAll.var.lambda,
                                        stochasticFire30.control.avgdens.var.lambda, stochasticFire30.perturb.state5.avgdens.var.lambda, stochasticFire30.perturb.state4_5.avgdens.var.lambda, stochasticFire30.perturb.state3_4_5.avgdens.var.lambda, stochasticFire30.perturb.state2_3_4_5.avgdens.var.lambda, stochasticFire30.perturb.stateAll.avgdens.var.lambda,
                                        periodicFire15.control.var.lambda, periodicFire15.perturb.state5.var.lambda, periodicFire15.perturb.state4_5.var.lambda, periodicFire15.perturb.state3_4_5.var.lambda, periodicFire15.perturb.state2_3_4_5.var.lambda, periodicFire15.perturb.stateAll.var.lambda,
                                        periodicFire15.control.avgdens.var.lambda, periodicFire15.perturb.state5.avgdens.var.lambda, periodicFire15.perturb.state4_5.avgdens.var.lambda, periodicFire15.perturb.state3_4_5.avgdens.var.lambda, periodicFire15.perturb.state2_3_4_5.avgdens.var.lambda, periodicFire15.perturb.stateAll.avgdens.var.lambda,
                                        periodicFire30.control.var.lambda, periodicFire30.perturb.state5.var.lambda, periodicFire30.perturb.state4_5.var.lambda, periodicFire30.perturb.state3_4_5.var.lambda, periodicFire30.perturb.state2_3_4_5.var.lambda, periodicFire30.perturb.stateAll.var.lambda,
                                        periodicFire30.control.avgdens.var.lambda, periodicFire30.perturb.state5.avgdens.var.lambda, periodicFire30.perturb.state4_5.avgdens.var.lambda, periodicFire30.perturb.state3_4_5.avgdens.var.lambda, periodicFire30.perturb.state2_3_4_5.avgdens.var.lambda, periodicFire30.perturb.stateAll.avgdens.var.lambda), 
                         ext.prob = c(stochasticFire15.control.sim[[3]], stochasticFire15.perturb.state5[[3]], stochasticFire15.perturb.state4_5[[3]], stochasticFire15.perturb.state3_4_5[[3]], stochasticFire15.perturb.state2_3_4_5[[3]], stochasticFire15.perturb.stateAll[[3]],
                                      stochasticFire15.control.sim.avgdens[[3]], stochasticFire15.perturb.state5.avgdens[[3]], stochasticFire15.perturb.state4_5.avgdens[[3]], stochasticFire15.perturb.state3_4_5.avgdens[[3]], stochasticFire15.perturb.state2_3_4_5.avgdens[[3]], stochasticFire15.perturb.stateAll.avgdens[[3]],
                                      stochasticFire30.control.sim[[3]], stochasticFire30.perturb.state5[[3]], stochasticFire30.perturb.state4_5[[3]], stochasticFire30.perturb.state3_4_5[[3]], stochasticFire30.perturb.state2_3_4_5[[3]], stochasticFire30.perturb.stateAll[[3]],
                                      stochasticFire30.control.sim.avgdens[[3]], stochasticFire30.perturb.state5.avgdens[[3]], stochasticFire30.perturb.state4_5.avgdens[[3]], stochasticFire30.perturb.state3_4_5.avgdens[[3]], stochasticFire30.perturb.state2_3_4_5.avgdens[[3]], stochasticFire30.perturb.stateAll.avgdens[[3]],
                                      periodicFire15.control.sim[[3]], periodicFire15.perturb.state5[[3]], periodicFire15.perturb.state4_5[[3]], periodicFire15.perturb.state3_4_5[[3]], periodicFire15.perturb.state2_3_4_5[[3]], periodicFire15.perturb.stateAll[[3]],
                                      periodicFire15.control.sim.avgdens[[3]], periodicFire15.perturb.state5.avgdens[[3]], periodicFire15.perturb.state4_5.avgdens[[3]], periodicFire15.perturb.state3_4_5.avgdens[[3]], periodicFire15.perturb.state2_3_4_5.avgdens[[3]], periodicFire15.perturb.stateAll.avgdens[[3]],
                                      periodicFire30.control.sim[[3]], periodicFire30.perturb.state5[[3]], periodicFire30.perturb.state4_5[[3]], periodicFire30.perturb.state3_4_5[[3]], periodicFire30.perturb.state2_3_4_5[[3]], periodicFire30.perturb.stateAll[[3]],
                                      periodicFire30.control.sim.avgdens[[3]], periodicFire30.perturb.state5.avgdens[[3]], periodicFire30.perturb.state4_5.avgdens[[3]], periodicFire30.perturb.state3_4_5.avgdens[[3]], periodicFire30.perturb.state2_3_4_5.avgdens[[3]], periodicFire30.perturb.stateAll.avgdens[[3]]))


# Reordering factor levels for plotting
sim.results$treatment = factor(sim.results$treatment, levels = c("Control", "TSF>3", "TSF3 and >3", "TSF2, 3, and >3", "TSF1, 2, 3, and >3", "All TSFs"))


## 2.2. Magnitude of changes between scenarios (mean and quantiles) ----
# -----------------------------------------------------------------

# Stochastic lambda
aggregate(stoch.lambda ~ treatment + disturbance + density, data = sim.results, FUN = "mean")
aggregate(stoch.lambda ~ treatment + disturbance + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))

cbind(aggregate(stoch.lambda ~ treatment + disturbance + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975))), mean = aggregate(stoch.lambda ~ treatment + disturbance + density, data = sim.results, FUN = "mean")$stoch.lambda)

# Variance
aggregate(var.lambda ~ treatment + disturbance + density, data = sim.results, FUN = "mean")
aggregate(var.lambda ~ treatment + disturbance + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))

aggregate(ext.prob ~ treatment + disturbance + density, data = sim.results, FUN = "mean")




###########################################################################
#
# 3. Plotting the results ----
#
###########################################################################

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim.results$density <- factor(sim.results$density, levels = c("DensDep", "AvgDens"))
sim.results$disturbance <- factor(sim.results$disturbance, levels = c("Stochastic fire (p = 0.066)", "Stochastic fire (p = 0.033)", "Periodic fire (f = 15 years)", "Periodic fire (f = 30 years)"))

facet_labels <- c(DensDep = "Density dependence", AvgDens = "Average density")

agg.stoch.lambda = aggregate(stoch.lambda ~ treatment + disturbance + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))
agg.var.lambda = aggregate(var.lambda ~ treatment + disturbance + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))

agg.stoch.lambda$treatment = factor(agg.stoch.lambda$treatment, levels = c("Control", "TSF>3", "TSF3 and >3", "TSF2, 3, and >3", "TSF1, 2, 3, and >3", "All TSFs"))
agg.stoch.lambda$density <- factor(agg.stoch.lambda$density, levels = c("DensDep", "AvgDens"))
agg.stoch.lambda$disturbance <- factor(agg.stoch.lambda$disturbance, levels = c("Stochastic fire (p = 0.066)", "Stochastic fire (p = 0.033)", "Periodic fire (f = 15 years)", "Periodic fire (f = 30 years)"))

agg.stoch.lambda = rbind(agg.stoch.lambda[which(agg.stoch.lambda$treatment == "Control" & agg.stoch.lambda$density == "DensDep"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF>3" & agg.stoch.lambda$density == "DensDep"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF3 and >3" & agg.stoch.lambda$density == "DensDep"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF2, 3, and >3" & agg.stoch.lambda$density == "DensDep"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF1, 2, 3, and >3" & agg.stoch.lambda$density == "DensDep"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "All TSFs" & agg.stoch.lambda$density == "DensDep"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "Control" & agg.stoch.lambda$density == "AvgDens"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF>3" & agg.stoch.lambda$density == "AvgDens"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF3 and >3" & agg.stoch.lambda$density == "AvgDens"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF2, 3, and >3" & agg.stoch.lambda$density == "AvgDens"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "TSF1, 2, 3, and >3" & agg.stoch.lambda$density == "AvgDens"), ],
                         agg.stoch.lambda[which(agg.stoch.lambda$treatment == "All TSFs" & agg.stoch.lambda$density == "AvgDens"), ])

agg.var.lambda$treatment = factor(agg.var.lambda$treatment, levels = c("Control", "TSF>3", "TSF3 and >3", "TSF2, 3, and >3", "TSF1, 2, 3, and >3", "All TSFs"))
agg.var.lambda$density <- factor(agg.var.lambda$density, levels = c("DensDep", "AvgDens"))
agg.var.lambda$disturbance <- factor(agg.var.lambda$disturbance, levels = c("Stochastic fire (p = 0.066)", "Stochastic fire (p = 0.033)", "Periodic fire (f = 15 years)", "Periodic fire (f = 30 years)"))

agg.var.lambda = rbind(agg.var.lambda[which(agg.var.lambda$treatment == "Control" & agg.var.lambda$density == "DensDep"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF>3" & agg.var.lambda$density == "DensDep"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF3 and >3" & agg.var.lambda$density == "DensDep"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF2, 3, and >3" & agg.var.lambda$density == "DensDep"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF1, 2, 3, and >3" & agg.var.lambda$density == "DensDep"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "All TSFs" & agg.var.lambda$density == "DensDep"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "Control" & agg.var.lambda$density == "AvgDens"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF>3" & agg.var.lambda$density == "AvgDens"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF3 and >3" & agg.var.lambda$density == "AvgDens"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF2, 3, and >3" & agg.var.lambda$density == "AvgDens"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "TSF1, 2, 3, and >3" & agg.var.lambda$density == "AvgDens"), ],
                         agg.var.lambda[which(agg.var.lambda$treatment == "All TSFs" & agg.var.lambda$density == "AvgDens"), ])

agg.ext.prob = aggregate(ext.prob ~ treatment + disturbance + density, data = sim.results, FUN = mean)
agg.ext.prob.quantiles = aggregate(ext.prob ~ treatment + disturbance + density, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


## 3.1. Stochastic lambda ----
# -----------------------

plot.stoch.lambda = ggplot(sim.results, aes(x = treatment, y = stoch.lambda, fill = disturbance)) +
  facet_wrap( ~ density, labeller = labeller(density = facet_labels)) +
  geom_boxplot(ymin = agg.stoch.lambda$stoch.lambda[, 1], ymax = agg.stoch.lambda$stoch.lambda[, 3], size = 0.1, outlier.size = 0.1) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean variance", group = disturbance), geom = "point", shape = 17, size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Perturbed state", y = "Stochastic log \u03BB") +
  scale_shape_manual("", values = c("Mean variance" = 17)) +
  scale_x_discrete(labels = c("Control", expression("TSF"[">3"]), expression("TSF"["3 and >3"]), expression("TSF"["2, 3, and >3"]), expression("TSF"["1, 2, 3, and >3"]), "All TSF")) +
  scale_fill_manual(name = "Disturbance type",
                    breaks = unique(sim.results$disturbance),
                    labels = unique(sim.results$disturbance),
                    values = cbbPalette[c(4, 6:8)]) +
  theme_bw() +
  theme(axis.title.x = element_text(size=9,colour="black", margin = margin(t=0,r=0,b=0,l=0)), 
        axis.title.y = element_text(size=9, colour = "black"), 
        axis.text.x = element_text(size=7,colour="black", margin = margin(t=2,r=0,b=0,l=0),angle=45,hjust=1), 
        axis.text.y = element_text(size=7,colour="black", margin = margin(t=0,r=2,b=0,l=0)), 
        legend.text = element_text(size=7), 
        legend.title = element_text(face="bold",size=9), 
        legend.position = "right", 
        legend.key.size = unit(1,"lines"),
        strip.text.x = element_text(size=7),
        panel.grid.major = element_line(size = 0.1, colour = "lightgrey"),
        panel.grid.minor = element_blank())

tiff(filename = "DP_StochasticLambda.tiff",
     width=6,
     height=3,
     units="in",
     bg="white",
     res=600,
     compression = "lzw")

plot.stoch.lambda

dev.off()


## 3.2. Combined plot (Fig. S10.3) ----
# --------------------------------

plot.legend.appendix = ggplot(agg.ext.prob, aes(x = treatment, y = ext.prob, colour = disturbance)) +
  facet_wrap( ~ density, labeller = labeller(density = facet_labels)) +
  geom_point(size = 3, position = position_dodge(1)) +
  labs(x = "Perturbed state", y = "Extinction probability") +
  scale_x_discrete(labels = c("Control", expression("TSF"[">3"]), expression("TSF"["3 and >3"]), expression("TSF"["2, 3, and >3"]), expression("TSF"["1, 2, 3, and >3"]), "All TSFs")) +
  scale_colour_manual(name = "Disturbance regime",
                    breaks = unique(sim.results$disturbance),
                    labels = c("Stochastic fire\n(p = 1/15)", "Stochastic fire\n(p = 1/30)", "Periodic fire\n(every 15 years)", "Periodic fire\n(every 30 years)"),
                    values = cbbPalette[c(4, 6:8)]) +
  ylim(0, 1) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7, margin = margin(t = 2)), 
        legend.title = element_text(size = 9), 
        legend.position = "right", 
        legend.key.size = unit(1.5, "lines"))


plot.var.lambda.appendix = ggplot(sim.results, aes(x = treatment, y = var.lambda, fill = disturbance)) +
  facet_wrap(~ density,labeller = labeller(density = facet_labels)) +
  geom_boxplot(ymin = agg.var.lambda$var.lambda[, 1], ymax = agg.var.lambda$var.lambda[, 3], outlier.size = 0.1, size = 0.1) +
  stat_summary(fun = mean, colour = cbbPalette[5], aes(shape = "Mean variance", group = disturbance), geom = "point", size = 1, position = position_dodge(width = 0.75)) +
  labs(x = "Perturbed state", y = "Variance of annual \u03BB") +
  scale_shape_manual("", values = c("Mean variance" = 17)) +
  scale_x_discrete(labels = c("Control", expression("TSF"[">3"]), expression("TSF"["3 and >3"]), expression("TSF"["2, 3, and >3"]), expression("TSF"["1, 2, 3, and >3"]), "All TSFs")) +
  scale_fill_manual(name = "Disturbance regime",
                    breaks = unique(sim.results$disturbance),
                    labels = unique(sim.results$disturbance),
                    values = cbbPalette[c(4, 6:8)]) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(2, "lines"), 
        strip.text.x = element_text(size = 7, face = "plain"))


plot.ext.prob.appendix = ggplot(agg.ext.prob, aes(x = treatment, y = ext.prob, colour = disturbance)) +
  facet_wrap( ~ density, labeller = labeller(density = facet_labels)) +
  geom_point(size = 2, position = position_dodge(1)) +
  labs(x = "Perturbed state", y = "Extinction probability") +
  scale_x_discrete(labels = c("Control", expression("TSF"[">3"]), expression("TSF"["3 and >3"]), expression("TSF"["2, 3, and >3"]), expression("TSF"["1, 2, 3, and >3"]), "All TSFs")) +
  scale_colour_manual(name = "Disturbance regime",
                    breaks = unique(sim.results$disturbance),
                    labels = unique(sim.results$disturbance),
                    values = cbbPalette[c(4, 6:8)]) +
  ylim(0, 1) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(2, "lines"), 
        strip.text.x = element_text(size = 7, face = "plain"))


dp.appendix = ggdraw() +
  draw_plot(plot.var.lambda.appendix, x = 0, y = 0.5, width = 0.7, height = 0.5) +
  draw_plot(plot.ext.prob.appendix, x = 0, y = 0, width = .7, height = 0.5) +
  draw_plot(get_legend(plot.legend.appendix), x = 0.75, y = 0.3, width = .2, height = 0.6) +
  draw_plot_label(label = c("(a)", "(b)"), size = 10, fontface = "plain",
                  x = c(0, 0), y = c(1, 0.505))


tiff(filename = "DPAppendix7.tiff",
    width=5,
    height=5,
    units="in",
    bg="white",
    res=800,
    compression = "lzw")

dp.appendix

dev.off()




###########################################################################
#
# 4. Saving files ----
#
###########################################################################

write.csv(sim.results,"SimResults_DP.csv",row.names = F)
write.csv(agg.stoch.lambda, "DP_StochLambda_MedianQuantiles.csv", row.names = F)
write.csv(agg.var.lambda, "DP_VarLambda_MedianQuantiles.csv", row.names = F)
