##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script contains the code to obtain the result plots for the marmot analysis. 
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

load("MarmotsProjectionsResults.RData")




###########################################################################
#
# 2. Analysing the results ----
#
###########################################################################

## 2.1. Building a data frame containing the metrics for each treatment ----
# ---------------------------------------------------------------------

sim.results = data.frame(treatment = c(rep("Control", 500), rep("High.Seas.Y.Surv", 500), rep("High.Seas.NR.Surv", 500), rep("High.Seas.R.Surv", 500), rep("Low.Seas.Y.Surv", 500), rep("Low.Seas.NR.Surv", 500), rep("Low.Seas.R.Surv", 500)), 
                         stoch.lambda = c(control.sim.lambda, high.seas.Y.surv.sim.lambda, high.seas.NR.surv.sim.lambda, high.seas.R.surv.sim.lambda, low.seas.Y.surv.sim.lambda, low.seas.NR.surv.sim.lambda, low.seas.R.surv.sim.lambda), 
                         var.lambda = c(control.var.lambda, high.seas.Y.surv.var.lambda, high.seas.NR.surv.var.lambda, high.seas.R.surv.var.lambda, low.seas.Y.surv.var.lambda, low.seas.NR.surv.var.lambda, low.seas.R.surv.var.lambda),
                         ext.prob = c(control.sim[[2]], high.seas.Y.surv[[2]], high.seas.NR.surv[[2]], high.seas.R.surv[[2]], low.seas.Y.surv[[2]], low.seas.NR.surv[[2]], low.seas.R.surv[[2]]),
                         vital.rate = c(rep("Control", 500), rep("Y survival", 500), rep("NR survival", 500), rep("R survival", 500), rep("Y survival", 500), rep("NR survival", 500), rep("R survival", 500)), 
                         scenario = c(rep("Control", 500), rep("High seasonality", 1500), rep("Low seasonality", 1500)))



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
# 3. Plotting the results
#
###########################################################################

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sim.results$vital.rate <- factor(sim.results$vital.rate, levels = c("Control", "Y survival", "NR survival", "R survival")) # Changing the order of the factors for plotting
sim.results$treatment <- factor(sim.results$treatment, levels = c("Control","High.Seas.Y.Surv","Low.Seas.Y.Surv","High.Seas.NR.Surv","Low.Seas.NR.Surv","High.Seas.R.Surv","Low.Seas.R.Surv")) # Changing the order of the factors for plotting

agg.marmots.stoch.lambda = aggregate(stoch.lambda ~ treatment, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))
agg.marmots.var.lambda = aggregate(var.lambda ~ treatment, data = sim.results, FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))


## 3.1. Stochastic lambda ----
# -----------------------

plot.stoch.lambda.fig = ggplot(sim.results, aes(x = vital.rate, y = stoch.lambda, fill = scenario)) +
  geom_boxplot(ymin = agg.marmots.stoch.lambda$stoch.lambda[, 1], ymax = agg.marmots.stoch.lambda$stoch.lambda[, 3]) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", shape = 17, size = 3, position = position_dodge(width = 0.75), show.legend = T) +
  labs(x = "", y = "Stochastic log \u03BB") +
  scale_fill_manual(values = cbbPalette,
                    name = "Scenario") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.text = element_text(size = 25), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none")


## 3.2. Var(lambda) ----
# -----------------

plot.var.lambda.fig = ggplot(sim.results, aes(x = vital.rate, y = var.lambda, fill = scenario)) +
  geom_boxplot(ymin = agg.marmots.var.lambda$var.lambda[, 1], ymax = agg.marmots.var.lambda$var.lambda[, 3]) +
  stat_summary(fun = mean, colour = "darkred", geom = "point", shape = 17, size = 3, position = position_dodge(width = 0.75), show.legend = T) +
  labs(x = "Perturbed vital rate", y = "Variance of annual \u03BB") +
  scale_fill_manual(values = cbbPalette,
                    name = "Scenario") +
  ylim(0, 0.07) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), angle = 20, hjust = 1), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "right", 
        legend.key.size = unit(4, "lines"))




###########################################################################
#
# 4. Plotting the absolute changes
# between seasons for NR survival (Figure S5.1)
#
###########################################################################

abs.changes = data.frame(year = as.numeric(rownames(coef(survR)$year)), 
                         abs.change = abs(inv.logit(coef(survR)$year$`(Intercept)` + coef(survR)$year$seasonwinter) - inv.logit(coef(survR)$year$`(Intercept)`)))
threshold.abs.changes = quantile(abs.changes$abs.change)[3]

plot.abs.change = ggplot(abs.changes, aes(x = year, y = abs.change)) +
  geom_point(size = 3) +
  geom_abline(intercept = threshold.abs.changes, slope = 0, colour = "darkred", lwd = 2) +
  ylab(expression(bold(Delta*"R survival (Summer-Winter)"))) +
  xlab("Year") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.text = element_text(size = 25), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none")

png(filename = "FigS5.png",
    width=4000,
    height=3000,
    units="px",
    bg="white",
    res=300,
    type = "cairo")

plot.abs.change

dev.off()




###########################################################################
#
# 5. Saving results file
#
###########################################################################

write.csv(sim.results,"SimResults_Marmots.csv",row.names = F)
