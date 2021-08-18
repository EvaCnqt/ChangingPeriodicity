############################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script uses the data of a dewy pine population in South-Eastern Spain collected between 2012 and 2019. 
# The aim of this script is to model the survival, transitions, and reproductive rates of four life-history stages: seedlings, juveniles, and small and large reproductive adults. We model the rates in different post-fire habitat states in natural populations (i.e., natural fire return) and populations under fire management policies and exposed to human disturbances, i.e., livestock trampling and grazing. In this script, we model the vital rates for time since fire (TSF) 1 to 3.
#
# Author: Eva Conquet
#
###########################################################################

###########################################################################
#
# 1. House keeping and loading libraries and data
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading libraries ----
# -----------------------

load.librairies = function(){
 library(lme4)
 library(bbmle)
 library(MASS)
 library(boot)
 library(ggplot2)
 library(multcomp)
 library(MuMIn)
}

load.librairies()


## 1.3. Loading and preparing data ----
# --------------------------------

data.droso = read.csv("DewyPinesTSF1_to_3_Data.csv")
head(data.droso)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Separate different sites & TSFs

   # Site C

siteC = data.droso[data.droso$site == "siteC", ]
siteC.TSF1 = data.droso[data.droso$site == "siteC" & data.droso$TSF == "one", ]
siteC.TSF2 = data.droso[data.droso$site == "siteC" & data.droso$TSF == "two", ]
siteC.TSF3 = data.droso[data.droso$site == "siteC" & data.droso$TSF == "three", ]

   # Site E

siteE = data.droso[data.droso$site == "siteE", ]
siteE.TSF1 = data.droso[data.droso$site == "siteE" & data.droso$TSF == "one", ]
siteE.TSF2 = data.droso[data.droso$site == "siteE" & data.droso$TSF == "two", ]
siteE.TSF3 = data.droso[data.droso$site == "siteE" & data.droso$TSF == "three", ]

   # Site F

siteF = data.droso[data.droso$site == "siteF", ]
siteF.TSF1 = data.droso[data.droso$site == "siteF" & data.droso$TSF == "one", ]
siteF.TSF2 = data.droso[data.droso$site == "siteF" & data.droso$TSF == "two", ]
siteF.TSF3 = data.droso[data.droso$site == "siteF" & data.droso$TSF == "three", ]


# Compute average density per square (i.e. population density)

   # Site C - TSF1

siteC.TSF1$IDsquare = NA
siteC.TSF1$IDsquare = paste(siteC.TSF1$transect, siteC.TSF1$subQuadrat, sep = "_")
nbsquares.siteC.TSF1 = length(unique(siteC.TSF1$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteC.TSF1, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean)
density.init.siteC.TSF1 = sum(density.per.square$density)/nbsquares.siteC.TSF1

   # Site C - TSF2

siteC.TSF2$IDsquare = NA
siteC.TSF2$IDsquare = paste(siteC.TSF2$transect, siteC.TSF2$subQuadrat, sep = "_")
nbsquares.siteC.TSF2 = length(unique(siteC.TSF2$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteC.TSF2, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean) 
density.init.siteC.TSF2 = sum(density.per.square$density)/nbsquares.siteC.TSF2

   # Site C - TSF3

siteC.TSF3$IDsquare = NA
siteC.TSF3$IDsquare = paste(siteC.TSF3$transect, siteC.TSF3$subQuadrat, sep = "_")
nbsquares.siteC.TSF3 = length(unique(siteC.TSF3$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteC.TSF3, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean) 
density.init.siteC.TSF3 = sum(density.per.square$density)/nbsquares.siteC.TSF3


   # Site E - TSF1

siteE.TSF1$IDsquare = NA
siteE.TSF1$IDsquare = paste(siteE.TSF1$transect, siteE.TSF1$subQuadrat, sep = "_")
nbsquares.siteE.TSF1 = length(unique(siteE.TSF1$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteE.TSF1, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean) 
density.init.siteE.TSF1 = sum(density.per.square$density)/nbsquares.siteE.TSF1

   # Site E - TSF2

siteE.TSF2$IDsquare = NA
siteE.TSF2$IDsquare = paste(siteE.TSF2$transect, siteE.TSF2$subQuadrat, sep = "_")
nbsquares.siteE.TSF2 = length(unique(siteE.TSF2$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteE.TSF2, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean) 
density.init.siteE.TSF2 = sum(density.per.square$density)/nbsquares.siteE.TSF2

   # Site E - TSF3

siteE.TSF3$IDsquare = NA
siteE.TSF3$IDsquare = paste(siteE.TSF3$transect, siteE.TSF3$subQuadrat, sep = "_")
nbsquares.siteE.TSF3 = length(unique(siteE.TSF3$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteE.TSF3, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean) 
density.init.siteE.TSF3 = sum(density.per.square$density)/nbsquares.siteE.TSF3


   # Site F - TSF1

siteF.TSF1$IDsquare = NA
siteF.TSF1$IDsquare = paste(siteF.TSF1$transect, siteF.TSF1$subQuadrat, sep = "_") # unique combination of transect and square
nbsquares.siteF.TSF1 = length(unique(siteF.TSF1$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteF.TSF1, function(x) unique(x)) # density per square
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean)
density.init.siteF.TSF1 = sum(density.per.square$density)/nbsquares.siteF.TSF1

   # Site F - TSF2

siteF.TSF2$IDsquare = NA
siteF.TSF2$IDsquare = paste(siteF.TSF2$transect, siteF.TSF2$subQuadrat, sep = "_")
nbsquares.siteF.TSF2 = length(unique(siteF.TSF2$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteF.TSF2, function(x) unique(x))
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean)
density.init.siteF.TSF2 = sum(density.per.square$density)/nbsquares.siteF.TSF2

   # Site F - TSF3

siteF.TSF3$IDsquare = NA
siteF.TSF3$IDsquare = paste(siteF.TSF3$transect, siteF.TSF3$subQuadrat, sep = "_")
nbsquares.siteF.TSF3 = length(unique(siteF.TSF3$IDsquare)) # Number of squares
density.per.square = aggregate(density ~ IDsquare + time, data = siteF.TSF3, function(x) unique(x))
density.per.square = aggregate(density ~ IDsquare, data = density.per.square, mean)
density.init.siteF.TSF3 = sum(density.per.square$density)/nbsquares.siteF.TSF3


nbsquares.per.site.per.TSF = expand.grid(site = unique(data.droso$site),
                                        TSF = unique(data.droso$TSF))
nbsquares.per.site.per.TSF$nbsquares = c(nbsquares.siteC.TSF1,
                                        nbsquares.siteE.TSF1,
                                        nbsquares.siteF.TSF1,
                                        nbsquares.siteC.TSF2,
                                        nbsquares.siteE.TSF2,
                                        nbsquares.siteF.TSF2,
                                        nbsquares.siteC.TSF3,
                                        nbsquares.siteE.TSF3,
                                        nbsquares.siteF.TSF3) 


# Changing the order of the TSF factor

data.droso$TSF = factor(data.droso$TSF, levels = c("one", "two", "three"))




###########################################################################
#
# 2. Fitting Generalized Linear Models (GLMs) ----
# to estimate the vital rates for TSF0 to TSF3 
#
###########################################################################

## 2.1. Survival (surv) ----
# ---------------------

## 2.1.1. Site C - Seedlings ----
# --------------------------

table(data.droso$surv[which(data.droso$site == "siteC" & data.droso$stage == "SD")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "SD")]) # There are very few or no observations in TSF2 and 3. We only fit a model for TSF1 and use the observed values for TSF2 and 3.

survSD.C1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SD" & data.droso$TSF == "one", ], family = binomial)
survSD.C2 = glm(surv ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SD" & data.droso$TSF == "one", ], family = binomial)
survSD.C3 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SD" & data.droso$TSF == "one", ], family = binomial)

AICctab(survSD.C1, survSD.C2, survSD.C3, base = TRUE) # survSD.C2 is the best model.

summary(survSD.C2)
survSD.siteC.TSF1 = survSD.C2
survSD.siteC.TSF2 = 0
survSD.siteC.TSF3 = 0


# Plotting the observed data vs the predictions of the model

survSD.siteC.obs = aggregate(surv ~ density, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "SD"), ], FUN = mean)

survSD.siteC.pred = data.frame(density = survSD.siteC.obs$density, 
                               pred = predict(survSD.siteC.TSF1, newdata = data.frame(density = survSD.siteC.obs$density), type = "response"))
colnames(survSD.siteC.pred) = c("density", "pred")


# Obs vs. Pred

plot(survSD.siteC.obs$surv, survSD.siteC.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SD) survival - Site C")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survSD.siteC.new.data = expand.grid(density = seq(1, 60))


# 95% confidence intervals

survSD.siteC.predict = predict(survSD.siteC.TSF1, survSD.siteC.new.data, se.fit = T)
survSD.siteC.new.data$pred = inv.logit(survSD.siteC.predict$fit)
survSD.siteC.new.data$lwr = inv.logit(survSD.siteC.predict$fit - (1.96 * survSD.siteC.predict$se.fit))
survSD.siteC.new.data$upr = inv.logit(survSD.siteC.predict$fit + (1.96 * survSD.siteC.predict$se.fit))


# Plotting the predictions

plot.survSD.siteC = ggplot(survSD.siteC.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Seedling survival - Site C, TSF1") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_SDSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSD.siteC

dev.off()


## 2.1.2. Site C - Juveniles ----
# --------------------------

table(data.droso$surv[which(data.droso$site == "siteC" & data.droso$stage == "J")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "J")]) # There are very few or no observations in TSF2 and 3. We only fit a model for TSF1 and use the observed values for TSF2 and 3.

survJ.C1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)
survJ.C2 = glm(surv ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)
survJ.C3 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)

AICctab(survJ.C1, survJ.C2, survJ.C3, base = T) # survJ.C1 and survJ.C2 are in the 2 dAIC is the best model.

summary(survJ.C1)
survJ.siteC.TSF1 = survJ.C1
survJ.siteC.TSF2 = mean(data.droso$surv[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "two"])
survJ.siteC.TSF3 = mean(data.droso$surv[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "three"])
survJ.siteC = aggregate(surv ~ TSF, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "J"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

survJ.siteC.new.data = data.frame(pred = NA, lwr = NA, upr = NA)

# 95% confidence intervals

survJ.siteC.predict = predict(survJ.siteC.TSF1, se.fit = T)
survJ.siteC.new.data$pred = unique(inv.logit(survJ.siteC.predict$fit))
survJ.siteC.new.data$lwr = unique(inv.logit(survJ.siteC.predict$fit - (1.96 * survJ.siteC.predict$se.fit)))
survJ.siteC.new.data$upr = unique(inv.logit(survJ.siteC.predict$fit + (1.96 * survJ.siteC.predict$se.fit)))


# Plotting the predictions

plot.survJ.siteC = ggplot(survJ.siteC.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Juvenile survival - Site C, TSF1") + 
 scale_x_discrete(labels = c("TSF1", "TSF2", "TSF3")) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_JSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survJ.siteC

dev.off()


## 2.1.3. Site C - Small reproductives ----
# ------------------------------------

table(data.droso$surv[which(data.droso$site == "siteC" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "SR")]) 

survSR.C1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C2 = glm(surv ~ TSF, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C3 = glm(surv ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C4 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C5 = glm(surv ~ TSF + density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C6 = glm(surv ~ TSF + density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C7 = glm(surv ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C8 = glm(surv ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
survSR.C9 = glm(surv ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)

AICctab(survSR.C1, survSR.C2, survSR.C3, survSR.C4, survSR.C5, survSR.C6, survSR.C7, survSR.C8, survSR.C9, base = T) # survSR.C5, survSR.C6, survSR.C7, and survSR.C9 are in the 2 dAIC range. survSR.C5 is the simplest model.

summary(survSR.C5)
survSR.siteC = survSR.C5


# Plotting the observed data vs the predictions of the model

survSR.siteC.obs = aggregate(surv ~ TSF + density, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "SR"), ], FUN = mean) 

survSR.siteC.pred = survSR.siteC.obs[, c(1, 2)]
survSR.siteC.pred$pred = predict(survSR.siteC, newdata = survSR.siteC.pred, type = "response")


# Obs vs. Pred

plot(survSR.siteC.obs$surv, survSR.siteC.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of small reproductives (SR) survival - Site C")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survSR.siteC.new.data = expand.grid(TSF = c("two", "three"), density = seq(0, 60))


# 95% confidence intervals

survSR.siteC.predict = predict(survSR.siteC, survSR.siteC.new.data, se.fit = T)
survSR.siteC.new.data$pred = inv.logit(survSR.siteC.predict$fit)
survSR.siteC.new.data$lwr = inv.logit(survSR.siteC.predict$fit - (1.96 * survSR.siteC.predict$se.fit))
survSR.siteC.new.data$upr = inv.logit(survSR.siteC.predict$fit + (1.96 * survSR.siteC.predict$se.fit))

survSR.siteC.new.data$TSF = factor(survSR.siteC.new.data$TSF, levels = c("two", "three")) # Changing the order of the TSF factor


# Plotting the predictions

plot.survSR.siteC = ggplot(survSR.siteC.new.data, aes(x = density, y = pred, group = TSF, colour = TSF)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr, fill = TSF), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Small reproductive survival - Site C") + 
 scale_fill_manual(name = "TSF", 
   breaks = c("two", "three"), 
   labels = c("TSF2", "TSF3"), 
   values = c(cbbPalette[2], cbbPalette[1])) + 
 scale_colour_manual(name = "TSF", 
   breaks = c("two", "three"), 
   labels = c("TSF2", "TSF3"), 
   values = c(cbbPalette[2], cbbPalette[1])) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_SRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSR.siteC

dev.off()


## 2.1.4. Site C - Large reproductives ----
# ------------------------------------

table(data.droso$surv[which(data.droso$site == "siteC" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "LR")]) # There are very few or no observations in TSF2 and 3. We only fit a model for TSF1 and use the observed values for TSF2 and 3.

survLR.C1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
survLR.C2 = glm(surv ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
survLR.C3 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)

AICctab(survLR.C1, survLR.C2, survLR.C3, base = T) # survLR.C1 and survLR.C2 are in the 2 dAIC range. survLR.C1 is the simplest model.

summary(survLR.C1)
survLR.siteC.TSF3 = survLR.C1
survLR.siteC.TSF2 = mean(data.droso$surv[which(data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "two")])
survLR.siteC = aggregate(surv ~ TSF, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "LR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

survLR.siteC.new.data = data.frame(pred = NA, lwr = NA, upr = NA)

# 95% confidence intervals

survLR.siteC.predict = predict(survLR.siteC.TSF3, se.fit = T)
survLR.siteC.new.data$pred = unique(inv.logit(survLR.siteC.predict$fit))
survLR.siteC.new.data$lwr = unique(inv.logit(survLR.siteC.predict$fit - (1.96 * survLR.siteC.predict$se.fit)))
survLR.siteC.new.data$upr = unique(inv.logit(survLR.siteC.predict$fit + (1.96 * survLR.siteC.predict$se.fit)))


# Plotting the predictions

plot.survLR.siteC = ggplot(survLR.siteC.new.data, aes(x = 1, y = pred)) + 
   geom_point(size = 3) + 
   geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
   xlab("") + 
   ylab("Survival probability") + 
   ggtitle("Large reproductive survival - Site C, TSF3") + 
   theme_bw() + 
   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_LRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survLR.siteC

dev.off()


## 2.1.5. Site E - Seedlings ----
# --------------------------

table(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "SD")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "SD")]) # There are very few observations in TSF2. We only fit a model for TSF1 and 3 and use the observed values for TSF2.

survSD.E1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E2 = glm(surv ~ TSF, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E3 = glm(surv ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E4 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E5 = glm(surv ~ TSF + density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E6 = glm(surv ~ TSF + density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E7 = glm(surv ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E8 = glm(surv ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)
survSD.E9 = glm(surv ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two", ], family = binomial)

AICctab(survSD.E1, survSD.E2, survSD.E3, survSD.E4, survSD.E5, survSD.E6, survSD.E7, survSD.E8, survSD.E9, base = TRUE) # survSD.E7 and survSD.E8 are in the 2 dAIC range. survJ.E7 is the simplest model.

summary(survSD.E7)
survSD.siteE.TSF1_3 = survSD.E7
survSD.siteE.TSF2 = mean(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF == "two")])


# Plotting the observed data vs the predictions of the model

survSD.siteE.obs = aggregate(surv ~ TSF + density, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "SD" & data.droso$TSF != "two"), ], FUN = mean)

survSD.siteE.pred = survSD.siteE.obs[, c(1, 2)]
survSD.siteE.pred$pred = predict(survSD.siteE.TSF1_3, newdata = survSD.siteE.pred, type = "response")


# Obs vs. Pred

plot(survSD.siteE.obs$surv, survSD.siteE.pred$pred, pch = 16, ylim = c(0, 0.7), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SD) survival - Site E")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survSD.siteE.new.data = expand.grid(TSF = c("one", "three"), density = seq(1, 60))


# 95% confidence intervals

survSD.siteE.predict = predict(survSD.siteE.TSF1_3, survSD.siteE.new.data, se.fit = T)
survSD.siteE.new.data$pred = inv.logit(survSD.siteE.predict$fit)
survSD.siteE.new.data$lwr = inv.logit(survSD.siteE.predict$fit - (1.96 * survSD.siteE.predict$se.fit))
survSD.siteE.new.data$upr = inv.logit(survSD.siteE.predict$fit + (1.96 * survSD.siteE.predict$se.fit))

survSD.siteE.new.data$TSF = factor(survSD.siteE.new.data$TSF, levels = c("one", "three")) # Changing the order of the TSF factor


# Plotting the predictions

plot.survSD.siteE = ggplot(survSD.siteE.new.data, aes(x = density, y = pred, colour = TSF)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr, fill = TSF), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Seedling survival - Site E") + 
 scale_fill_manual(name = "TSF", 
                   breaks = unique(survSD.siteE.new.data$TSF), 
                   labels = c("TSF1", "TSF3"), 
                   values = c(cbbPalette[2], cbbPalette[7])) + 
 scale_colour_manual(name = "TSF", 
                     breaks = unique(survSD.siteE.new.data$TSF), 
                     labels = c("TSF1", "TSF3"), 
                     values = c(cbbPalette[2], cbbPalette[7])) +  
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_SDSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSD.siteE

dev.off()


## 2.1.6. Site E - Juveniles ----
# --------------------------

table(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "J")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "J")]) 

survJ.E1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E2 = glm(surv ~ TSF, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E3 = glm(surv ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E4 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E5 = glm(surv ~ TSF + density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E6 = glm(surv ~ TSF + density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E7 = glm(surv ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E8 = glm(surv ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
survJ.E9 = glm(surv ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)

AICctab(survJ.E1, survJ.E2, survJ.E3, survJ.E4, survJ.E5, survJ.E6, survJ.E7, survJ.E8, base = T) # survJ.E6, and survJ.E8 are in the 2 dAIC range. survJ.E6 is the simplest model.

summary(survJ.E6) 
survJ.siteE = glm(surv ~ TSF + density + I(density^2), data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)


# Plotting the observed data vs the predictions of the model

survJ.siteE.obs = aggregate(surv ~ TSF + density, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "J"), ], FUN = mean) 

survJ.siteE.pred = survJ.siteE.obs[, c(1, 2)]
survJ.siteE.pred$pred = predict(survJ.siteE, newdata = survJ.siteE.pred, type = "response")


# Obs vs. Pred

plot(survJ.siteE.obs$surv, survJ.siteE.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of juveniles (J) survival - Site E")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survJ.siteE.new.data = expand.grid(TSF = levels(data.droso$TSF), density = seq(0, 60))


# 95% confidence intervals

survJ.siteE.predict = predict(survJ.siteE, survJ.siteE.new.data, se.fit = T)
survJ.siteE.new.data$pred = inv.logit(survJ.siteE.predict$fit)
survJ.siteE.new.data$lwr = inv.logit(survJ.siteE.predict$fit - (1.96 * survJ.siteE.predict$se.fit))
survJ.siteE.new.data$upr = inv.logit(survJ.siteE.predict$fit + (1.96 * survJ.siteE.predict$se.fit))

survJ.siteE.new.data$TSF = factor(survJ.siteE.new.data$TSF, levels = c("one", "two", "three")) # Changing the order of the TSF factor


# Plotting the predictions 

plot.survJ.siteE = ggplot(survJ.siteE.new.data, aes(x = density, y = pred, group = TSF, colour = TSF)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr, fill = TSF), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Juvenile survival - Site E") + 
 scale_fill_manual(name = "TSF", 
   breaks = levels(data.droso$TSF), 
   labels = c("TSF1", "TSF2", "TSF3"), 
   values = c(cbbPalette[2], cbbPalette[1], cbbPalette[7])) + 
 scale_colour_manual(name = "TSF", 
   breaks = levels(data.droso$TSF), 
   labels = c("TSF1", "TSF2", "TSF3"), 
   values = c(cbbPalette[2], cbbPalette[1], cbbPalette[7])) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_JSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survJ.siteE

dev.off()


## 2.1.7. Site E - Small reproductives ----
# ------------------------------------

table(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "SR")]) # There are very few observations in TSF2. We only fit a model for TSF3 and use the observed values for TSF2.

survSR.E1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)
survSR.E2 = glm(surv ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)
survSR.E3 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)

AICctab(survSR.E1, survSR.E2, survSR.E3, base = T) # survSR.E2 and survSR.E3 are in the 2 dAIC range. survSR.E2 is the simplest model.

summary(survSR.E2)
survSR.siteE.TSF3 = survSR.E2
survSR.siteE.TSF2 = mean(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "two")])


# Plotting the observed data vs the predictions of the model

survSR.siteE.obs = aggregate(surv ~ density, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three"), ], FUN = mean) 

survSR.siteE.pred = data.frame(density = survSR.siteE.obs[, 1])
survSR.siteE.pred$pred = predict(survSR.siteE.TSF3, newdata = survSR.siteE.pred, type = "response")


# Obs vs. Pred

plot(survSR.siteE.obs$surv, survSR.siteE.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of small reproductive adult (SR) survival - Site E")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survSR.siteE.new.data = expand.grid(density = seq(0, 60))


# 95% confidence intervals

survSR.siteE.predict = predict(survSR.siteE.TSF3, survSR.siteE.new.data, se.fit = T)
survSR.siteE.new.data$pred = inv.logit(survSR.siteE.predict$fit)
survSR.siteE.new.data$lwr = inv.logit(survSR.siteE.predict$fit - (1.96 * survSR.siteE.predict$se.fit))
survSR.siteE.new.data$upr = inv.logit(survSR.siteE.predict$fit + (1.96 * survSR.siteE.predict$se.fit))


# Plot the predictions

plot.survSR.siteE = ggplot(survSR.siteE.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Density") + 
 ylab("Survival probability") + 
 ggtitle("Small reproductive survival - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_SRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSR.siteE

dev.off()


## 2.1.8. Site E - Large reproductives ----
# ------------------------------------

table(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "LR")]) # There are very few observations in TSF2. We only fit a model for TSF3 and use the observed values for TSF2.

survLR.E1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
survLR.E2 = glm(surv ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
survLR.E3 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)

AICctab(survLR.E1, survLR.E2, survLR.E3, base = T) # survLR.E2 and survLR.E3 are in the 2 dAIC range. survLR.E2 is the simplest model.

summary(survLR.E2)

survLR.siteE.TSF3 = survLR.E2
survLR.siteE.TSF2 = mean(data.droso$surv[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "two")])

# Plotting the observed data vs the predictions of the model
survLR.siteE.obs = aggregate(surv ~ density, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three"), ], FUN = mean)

survLR.siteE.pred = data.frame(density = survLR.siteE.obs[, 1])
survLR.siteE.pred$pred = predict(survLR.siteE.TSF3, newdata = survLR.siteE.pred, type = "response")


# Obs vs. Pred

plot(survLR.siteE.obs$surv, survLR.siteE.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of large reproductives (LR) survival - Site E")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survLR.siteE.new.data = expand.grid(TSF = c("three"), density = seq(0, 60))

# 95% confidence intervals
survLR.siteE.predict = predict(survLR.siteE.TSF3, survLR.siteE.new.data, se.fit = T)
survLR.siteE.new.data$pred = inv.logit(survLR.siteE.predict$fit)
survLR.siteE.new.data$lwr = inv.logit(survLR.siteE.predict$fit - (1.96 * survLR.siteE.predict$se.fit))
survLR.siteE.new.data$upr = inv.logit(survLR.siteE.predict$fit + (1.96 * survLR.siteE.predict$se.fit))


# Plotting the predictions

plot.survLR.siteE = ggplot(survLR.siteE.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Large reproductive survival - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_LRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survLR.siteE

dev.off()



## 2.1.9. Site F - Seedlings ----
# -------------------------- 

table(data.droso$surv[which(data.droso$site == "siteF" & data.droso$stage == "SD")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "SD")]) # There are very few observations in TSF1 and none in TSF2. We only fit a model for TSF3 and use the observed values for TSF1 and 2.

survSD.F1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "SD" & data.droso$TSF == "three", ], family = binomial)
survSD.F2 = glm(surv ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "SD" & data.droso$TSF == "three", ], family = binomial)
survSD.F3 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "SD" & data.droso$TSF == "three", ], family = binomial)

AICctab(survSD.F1, survSD.F2, survSD.F3, base = TRUE) # survSD.F2 and survSD.F3 are in the 2 dAIC range. survJ.F2 is the simplest model.

summary(survSD.F2)
survSD.siteF.TSF3 = survSD.F2
survSD.siteF.TSF1 = mean(data.droso$surv[which(data.droso$site == "siteF" & data.droso$stage == "SD" & data.droso$TSF == "one")])
survSD.siteF.TSF2 = 0


# Plotting the observed data vs the predictions of the model

survSD.siteF.obs = aggregate(surv ~ density, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "SD" & data.droso$TSF == "three"), ], FUN = mean)

survSD.siteF.pred = data.frame(density = survSD.siteF.obs$density, 
                               pred = predict(survSD.siteF.TSF3, newdata = data.frame(density = survSD.siteF.obs$density), type = "response"))


# Obs vs. Pred

plot(survSD.siteF.obs$surv, survSD.siteF.pred$pred, pch = 16, ylim = c(0, 0.7), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SD) survival - Site F")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survSD.siteF.new.data = expand.grid(density = seq(0, 60))

# 95% confidence intervals
survSD.siteF.predict = predict(survSD.siteF.TSF3, survSD.siteF.new.data, se.fit = T)
survSD.siteF.new.data$pred = inv.logit(survSD.siteF.predict$fit)
survSD.siteF.new.data$lwr = inv.logit(survSD.siteF.predict$fit - (1.96 * survSD.siteF.predict$se.fit))
survSD.siteF.new.data$upr = inv.logit(survSD.siteF.predict$fit + (1.96 * survSD.siteF.predict$se.fit))


# Plotting the predictions

plot.survSD.siteF = ggplot(survSD.siteF.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Seedling survival - Site F") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_SDSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSD.siteF

dev.off()


## 2.1.10. Site F - Juveniles ----
# ---------------------------

table(data.droso$surv[which(data.droso$site == "siteF" & data.droso$stage == "J")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "J")]) # There are very few observations in TSF2. We only fit a model for TSF1 and 3 and use the observed values for TSF2.

survJ.F1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F2 = glm(surv ~ TSF, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F3 = glm(surv ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F4 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F5 = glm(surv ~ TSF + density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F6 = glm(surv ~ TSF + density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F7 = glm(surv ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F8 = glm(surv ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)
survJ.F9 = glm(surv ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two", ], family = binomial)

AICctab(survJ.F1, survJ.F2, survJ.F3, survJ.F4, survJ.F5, survJ.F6, survJ.F7, survJ.F8, survJ.F9, base = T) # survJ.F1, survJ.F2, and survJ.F3 are in the 2 dAIC range. survJ.F1 is the simplest model.

summary(survJ.F1) 
survJ.siteF.TSF1_3 = survJ.F1
survJ.siteF.TSF2 = mean(data.droso$surv[which(data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF == "two")])
survJ.siteF = aggregate(surv ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "J"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

survJ.siteF.new.data = expand.grid(TSF = c("one", "three"))
survJ.siteF.predict = predict(survJ.siteF.TSF1_3, newdata = survJ.siteF.new.data, type = "response", se.fit = T)

# 95% confidence intervals

survJ.siteF.new.data$pred = inv.logit(survJ.siteF.predict$fit)
survJ.siteF.new.data$lwr = inv.logit(survJ.siteF.predict$fit - (1.96 * survJ.siteF.predict$se.fit))
survJ.siteF.new.data$upr = inv.logit(survJ.siteF.predict$fit + (1.96 * survJ.siteF.predict$se.fit))


# Plotting the predictions

plot.survJ.siteF = ggplot(survJ.siteF.new.data, aes(x = TSF, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 ylab("Survival probability") + 
 ggtitle("Juvenile survival - Site F") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_JSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survJ.siteF

dev.off()


## 2.1.11. Site F - Small reproductives ----
# -------------------------------------

table(data.droso$surv[which(data.droso$site == "siteF" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "SR")]) # There are very few observations all TSFs. We use the observed values.

survSR.siteF = aggregate(surv ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "SR"), ], FUN = mean)


# Plotting the predictions

plot.survSR.siteF = ggplot(survSR.siteF, aes(x = TSF, y = surv)) + 
 geom_point(size = 3) + 
 ylab("Survival probability") + 
 ggtitle("Small reproductive survival - Site F") + 
 scale_x_discrete(labels = c("TSF2", "TSF3")) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_SRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSR.siteF

dev.off()


## 2.1.12. Site F - Large reproductives ----
# -------------------------------------

table(data.droso$surv[which(data.droso$site == "siteF" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "LR")]) 

survLR.F1 = glm(surv ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F2 = glm(surv ~ TSF, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F3 = glm(surv ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F4 = glm(surv ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F5 = glm(surv ~ TSF + density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F6 = glm(surv ~ TSF + density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F7 = glm(surv ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F8 = glm(surv ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
survLR.F9 = glm(surv ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)

AICctab(survLR.F1, survLR.F2, survLR.F3, survLR.F4, survLR.F5, survLR.F6, survLR.F7, survLR.F8, survLR.F9, base = T) # survLR.F2 and survLR.F5 are in the 2 dAIC range. survLR.F2 is the simplest model.

summary(survLR.F2)
survLR.siteF = survLR.F2


# Plotting the observed data vs the predictions of the model

survLR.siteF.obs = aggregate(surv ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "LR"), ], FUN = mean)

survLR.siteF.pred = data.frame(TSF = survLR.siteF.obs$TSF, 
                               pred = predict(survLR.F2, newdata = data.frame(TSF = survLR.siteF.obs$TSF), type = "response"))


# Obs vs. Pred

plot(survLR.siteF.obs$surv, survLR.siteF.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of large reproductives (LR) survival - Site F")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survLR.siteF.new.data = expand.grid(TSF = c("two", "three"))

# 95% confidence intervals

survLR.siteF.predict = predict(survLR.siteF, survLR.siteF.new.data, se.fit = T)
survLR.siteF.new.data$pred = inv.logit(survLR.siteF.predict$fit)
survLR.siteF.new.data$lwr = inv.logit(survLR.siteF.predict$fit - (1.96 * survLR.siteF.predict$se.fit))
survLR.siteF.new.data$upr = inv.logit(survLR.siteF.predict$fit + (1.96 * survLR.siteF.predict$se.fit))

survLR.siteF.new.data$TSF = factor(survLR.siteF.new.data$TSF, levels = c("two", "three")) # Changing the order of the TSF factor


# Plotting the predictions

plot.survLR.siteF = ggplot(survLR.siteF.new.data, aes(x = TSF, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("Time Since Fire") + 
 ylab("Survival probability") + 
 ggtitle("Large reproductive survival - Site F") + 
 scale_x_discrete(labels = c("TSF2", "TSF3")) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_LRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survLR.siteF

dev.off()



## 2.2. Transition from juvenile (J) to either small reproductive (SR) or large reproductive (LR) (transitionJ) ----
# -------------------------------------------------------------------------------------------------------------

## 2.2.1. Site C ----
# --------------

table(data.droso$transitionJ[which(data.droso$site == "siteC" & data.droso$stage == "J")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "J")]) # There are very few observations in TSF2 and 3. We only fit a model for TSF1 and use the observed values for TSF2 and 3.

transJ.C1 = glm(transitionJ ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)
transJ.C2 = glm(transitionJ ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)
transJ.C3 = glm(transitionJ ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)

AICctab(transJ.C1, transJ.C2, transJ.C3, base = T) # transJ.C1, transJ.C2, and transJ.C3 are in the 2 dAIC range. transJ.C1 is the simplest model, but we can also check the other to see whether it brings something to the predictions.

summary(transJ.C1)
transJ.siteC.TSF1 = transJ.C1
transJ.siteC.TSF2 = mean(data.droso$transitionJ[which(data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "two")], na.rm = T)
transJ.siteC.TSF3 = mean(data.droso$transitionJ[which(data.droso$site == "siteC" & data.droso$stage == "J" & data.droso$TSF == "three")], na.rm = T)
transJ.siteC = aggregate(transitionJ ~ TSF, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "J"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

transJ.siteC.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

transJ.siteC.predict = predict(transJ.siteC.TSF1, se.fit = T)
transJ.siteC.new.data$pred = unique(inv.logit(transJ.siteC.predict$fit))
transJ.siteC.new.data$lwr = unique(inv.logit(transJ.siteC.predict$fit - (1.96 * transJ.siteC.predict$se.fit)))
transJ.siteC.new.data$upr = unique(inv.logit(transJ.siteC.predict$fit + (1.96 * transJ.siteC.predict$se.fit)))


# Plotting the predictions

plot.transJ.siteC = ggplot(transJ.siteC.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Transition probability") + 
 ggtitle("Juvenile to large adult transition - Site C, TSF1") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_JTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transJ.siteC

dev.off()


## 2.2.2. Site E ----
# --------------

table(data.droso$transitionJ[which(data.droso$site == "siteE" & data.droso$stage == "J")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "J")])

transJ.E1 = glm(transitionJ ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E2 = glm(transitionJ ~ TSF, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E3 = glm(transitionJ ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E4 = glm(transitionJ ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E5 = glm(transitionJ ~ TSF + density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E6 = glm(transitionJ ~ TSF + density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E7 = glm(transitionJ ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E8 = glm(transitionJ ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)
transJ.E9 = glm(transitionJ ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)

AICctab(transJ.E1, transJ.E2, transJ.E3, transJ.E4, transJ.E5, transJ.E6, transJ.E7, transJ.E8, base = T) # transJ.E6 and transJ.E8 are in the 2 dAIC range. transJ.E6 is the simplest model.

summary(transJ.E6)

transJ.siteE = glm(transitionJ ~ TSF + density + I(density^2), data = data.droso[data.droso$site == "siteE" & data.droso$stage == "J", ], family = binomial)


# Plotting the observed data vs the predictions of the model

transJ.siteE.obs = aggregate(transitionJ ~ TSF + density, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "J"), ], FUN = mean)

transJ.siteE.pred = transJ.siteE.obs[, c(1, 2)]
transJ.siteE.pred$pred = predict(transJ.siteE, newdata = transJ.siteE.pred, type = "response")


# Obs vs. Pred

plot(transJ.siteE.obs$transitionJ, transJ.siteE.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted transition probability", main = "Observed vs. Predicted values of J to adult transition - Site E")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transJ.siteE.new.data = expand.grid(TSF = levels(data.droso$TSF), 
     density = seq(0, 60))


# 95% confidence intervals

transJ.siteE.predict = predict(transJ.siteE, transJ.siteE.new.data, se.fit = T)
transJ.siteE.new.data$pred = inv.logit(transJ.siteE.predict$fit)
transJ.siteE.new.data$lwr = inv.logit(transJ.siteE.predict$fit-(1.96*transJ.siteE.predict$se.fit))
transJ.siteE.new.data$upr = inv.logit(transJ.siteE.predict$fit + (1.96*transJ.siteE.predict$se.fit))

transJ.siteE.new.data$TSF = factor(transJ.siteE.new.data$TSF, levels = c("one", "two", "three"))


# Plotting the predictions

plot.transJ.siteE = ggplot(transJ.siteE.new.data, aes(x = density, y = pred, group = TSF, colour = TSF)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr, fill = TSF), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Transition probability") + 
 ggtitle("Juvenile to large adult transition - Site E") + 
 scale_fill_manual(name = "TSF", 
   breaks = c("one", "two", "three"), 
   labels = c("TSF1", "TSF2", "TSF3"), 
   values = c(cbbPalette[2], cbbPalette[1], cbbPalette[7])) + 
 scale_colour_manual(name = "TSF", 
   breaks = c("one", "two", "three"), 
   labels = c("TSF1", "TSF2", "TSF3"), 
   values = c(cbbPalette[2], cbbPalette[1], cbbPalette[7])) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_JTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transJ.siteE

dev.off()


## 2.2.3. Site F ----
# --------------

table(data.droso$transitionJ[which(data.droso$site == "siteF" & data.droso$stage == "J")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "J")]) # There are very few observations in TSF2 and 3. We only fit a model for TSF1 and use the observed values for TSF2 and 3.

transJ.F1 = glm(transitionJ ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)
transJ.F2 = glm(transitionJ ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)
transJ.F3 = glm(transitionJ ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF == "one", ], family = binomial)

AICctab(transJ.F1, transJ.F2, base = T) # transJ.F1 and transJ.F2 are in the 2 dAIC range. transJ.F1 is the simplest model.

summary(transJ.F1)
transJ.siteF.TSF1 = transJ.F1
transJ.siteF.TSF2 = 0
transJ.siteF.TSF3 = mean(data.droso$transitionJ[which(data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF == "three")], na.rm = T)
transJ.siteF = aggregate(transitionJ ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "J" & data.droso$TSF != "two"), ], FUN = mean)
transJ.siteF = rbind(transJ.siteF, c("two", 0.0))


# Computing the predictions of the model and the 95 % CI

transJ.siteF.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

transJ.siteF.predict = predict(transJ.siteF.TSF1, se.fit = T)
transJ.siteF.new.data$pred = unique(inv.logit(transJ.siteF.predict$fit))
transJ.siteF.new.data$lwr = unique(inv.logit(transJ.siteF.predict$fit - (1.96 * transJ.siteF.predict$se.fit)))
transJ.siteF.new.data$upr = unique(inv.logit(transJ.siteF.predict$fit + (1.96 * transJ.siteF.predict$se.fit)))


# Plotting the predictions

plot.transJ.siteF = ggplot(transJ.siteF.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Transition probability") + 
 ggtitle("Juvenile to large adult transition - Site F, TSF1") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_JTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transJ.siteF

dev.off()


## 2.3. Transition from small reproductive (SR) to either small reproductive (SR = stasis) or large reproductive (LR = growth) (transitionSR) ----
# ---------------------------------------------------------------------------------------------------------------------------------------------

## 2.3.1. Site C ----
# --------------

table(data.droso$transitionSR[which(data.droso$site == "siteC" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "SR")]) # There are very few observations in all TSFs. We use the observations.

transSR.siteC = aggregate(transitionSR ~ TSF, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "SR"), ], FUN = mean)


# Plotting the predictions

plot.transSR.siteC = ggplot(transSR.siteC, aes(x = TSF, y = transitionSR)) + 
 geom_point(size = 3) + 
 xlab("TSF") + 
 ylab("Transition probability") + 
 ggtitle("Small to large adult transition - Site C") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_SRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transSR.siteC

dev.off()


## 2.3.2. Site E ----
# --------------

table(data.droso$transitionSR[which(data.droso$site == "siteE" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "SR")]) # There are very few observations in TSF2. We can fit a model for TSF3 and use the observations for TSF2.

transSR.E1 = glm(transitionSR ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)
transSR.E2 = glm(transitionSR ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)
transSR.E3 = glm(transitionSR ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)

AICctab(transSR.E1, transSR.E2, transSR.E3, base = T) # transSR.E1, transSR.E2, and transSR.E3 are in the 2 dAIC range. transSR.E1 is the simplest model.

summary(transSR.E1)
transSR.siteE.TSF3 = transSR.E1
transSR.siteE.TSF2 = mean(data.droso$transitionSR[which(data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "two")])
transSR.siteE = aggregate(transitionSR ~ TSF, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "SR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

transSR.siteE.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

transSR.siteE.predict = predict(transSR.siteE.TSF3, se.fit = T)
transSR.siteE.new.data$pred = unique(inv.logit(transSR.siteE.predict$fit))
transSR.siteE.new.data$lwr = unique(inv.logit(transSR.siteE.predict$fit - (1.96 * transSR.siteE.predict$se.fit)))
transSR.siteE.new.data$upr = unique(inv.logit(transSR.siteE.predict$fit + (1.96 * transSR.siteE.predict$se.fit)))


# Plotting the predictions

plot.transSR.siteE = ggplot(transSR.siteE.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Transition probability") + 
 ggtitle("Small to large adult transition - Site E, TSF2") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_SRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transSR.siteE

dev.off()


## 2.3.3. Site F ----
# --------------

table(data.droso$transitionSR[which(data.droso$site == "siteF" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "SR")]) # There are very few observations in all TSFs. We can use the observations.

transSR.siteF = aggregate(transitionSR ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "SR"), ], FUN = mean)


# Plotting the predictions

plot.transSR.siteF = ggplot(transSR.siteF, aes(x = TSF, y = transitionSR)) + 
 geom_point(size = 3) + 
 xlab("Time Since Fire (TSF)") + 
 ylab("Transition probability") + 
 ggtitle("Small to large adult transition - Site F") + 
 scale_x_discrete(labels = c("TSF2", "TSF3")) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_SRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transSR.siteF

dev.off()


## 2.4. Transition from large reproductive (LR) to either small reproductive (SR = shrinkage) or large reproductive (LR = stasis) (transitionSR) ----
# ---------------------------------------------------------------------------------------------------------------------------------------------

## 2.4.1. Site C ----
# --------------

table(data.droso$transitionLR[which(data.droso$site == "siteC" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "LR")]) # There are very few observations in TSF2. We can fit a model for TSF3 and use the observations for TSF2.

transLR.C1 = glm(transitionLR ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
transLR.C2 = glm(transitionLR ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
transLR.C3 = glm(transitionLR ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)

AICctab(transLR.C1, transLR.C2, transLR.C3, base = T) # transLR.C2 is the best model.

summary(transLR.C2)
transLR.siteC.TSF3 = transLR.C2
transLR.siteC.TSF2 = 0


# Plotting the observed data vs the predictions of the model

transLR.siteC.obs = aggregate(transitionLR ~ density, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three"), ], FUN = mean)

transLR.siteC.pred = data.frame(density = transLR.siteC.obs[, 1])
transLR.siteC.pred$pred = predict(transLR.siteC.TSF3, newdata = transLR.siteC.pred, type = "response")


# Obs vs. Pred

plot(transLR.siteC.obs$transitionLR, transLR.siteC.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted transition probability", main = "Observed vs. Predicted values of LR to LR transition - Site C")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transLR.siteC.new.data = expand.grid(density = seq(0, 60))


# 95% confidence intervals

transLR.siteC.predict = predict(transLR.siteC.TSF3, transLR.siteC.new.data, se.fit = T)
transLR.siteC.new.data$pred = inv.logit(transLR.siteC.predict$fit)
transLR.siteC.new.data$lwr = inv.logit(transLR.siteC.predict$fit - (1.96 * transLR.siteC.predict$se.fit))
transLR.siteC.new.data$upr = inv.logit(transLR.siteC.predict$fit + (1.96 * transLR.siteC.predict$se.fit))


# Plotting the predictions

plot.transLR.siteC = ggplot(transLR.siteC.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Transition probability") + 
 ggtitle("Large to small adult transition - Site C, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_LRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transLR.siteC

dev.off()


## 2.4.2. Site E ----
# --------------

table(data.droso$transitionLR[which(data.droso$site == "siteE" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "LR")]) # There are very few observations in TSF2. We can fit a model for TSF3 and use the observations for TSF2.

transLR.E1 = glm(transitionLR ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
transLR.E2 = glm(transitionLR ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
transLR.E3 = glm(transitionLR ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)

AICctab(transLR.E1, transLR.E2, transLR.E3, base = T) # transLR.E3 is the best model.

summary(transLR.E3)

transLR.siteE.TSF3 = glm(transitionLR ~ density + I(density^2), data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
transLR.siteE.TSF2 = mean(data.droso$transitionLR[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "two")], na.rm = T)


# Plotting the observed data vs the predictions of the model

transLR.siteE.obs = aggregate(transitionLR ~ density, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three"), ], FUN = mean)

transLR.siteE.pred = data.frame(density = transLR.siteE.obs[, 1])
transLR.siteE.pred$pred = predict(transLR.siteE.TSF3, newdata = transLR.siteE.pred, type = "response")


# Obs vs. Pred

plot(transLR.siteE.obs$transitionLR, transLR.siteE.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted transition probability", main = "Observed vs. Predicted values of LR to LR transition - Site E")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transLR.siteE.new.data = expand.grid(density = seq(0, 60))


# 95% confidence intervals

transLR.siteE.predict = predict(transLR.siteE.TSF3, transLR.siteE.new.data, se.fit = T)
transLR.siteE.new.data$pred = inv.logit(transLR.siteE.predict$fit)
transLR.siteE.new.data$lwr = inv.logit(transLR.siteE.predict$fit - (1.96 * transLR.siteE.predict$se.fit))
transLR.siteE.new.data$upr = inv.logit(transLR.siteE.predict$fit + (1.96 * transLR.siteE.predict$se.fit))


# Plotting the predictions

plot.transLR.siteE = ggplot(transLR.siteE.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Transition probability") + 
 ggtitle("Large to small adult transition - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_LRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transLR.siteE

dev.off()


## 2.4.3. Site F ----
# --------------

table(data.droso$transitionLR[which(data.droso$site == "siteF" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "LR")]) # There are only 0s in TSF3. We can use a model for TSF2 and use the observations for TSF3.

transLR.F1 = glm(transitionLR ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "two", ], family = binomial)
transLR.F2 = glm(transitionLR ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "two", ], family = binomial)
transLR.F3 = glm(transitionLR ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "two", ], family = binomial)

AICctab(transLR.F1, transLR.F2, base = T) # transLR.F1 and transLR.F2 are in the 2 dAIC range. transLR.F1 is the simplest model.

summary(transLR.F1)
transLR.siteF.TSF2 = transLR.F1
transLR.siteF.TSF3 = 0
transLR.siteF = aggregate(transitionLR ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "LR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

transLR.siteF.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

transLR.siteF.predict = predict(transLR.siteF.TSF2, se.fit = T)
transLR.siteF.new.data$pred = unique(inv.logit(transLR.siteF.predict$fit))
transLR.siteF.new.data$lwr = unique(inv.logit(transLR.siteF.predict$fit - (1.96 * transLR.siteF.predict$se.fit)))
transLR.siteF.new.data$upr = unique(inv.logit(transLR.siteF.predict$fit + (1.96 * transLR.siteF.predict$se.fit)))


# Plotting the predictions

plot.transLR.siteF = ggplot(transLR.siteF.new.data, aes(x = 1, y = pred)) + 
   geom_point(size = 3) + 
   geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
   xlab("") + 
   ylab("Transition probability") + 
   ggtitle("Large to small adult transition - Site F, TSF2") + 
   theme_bw() + 
   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_LRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transLR.siteF

dev.off()


## 2.5. Flowering probability (fl) ----
# --------------------------------

## 2.5.1. Site C - Small reproductive ----
# -----------------------------------

table(data.droso$fl[which(data.droso$site == "siteC" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "SR")]) 

floweringprobSR.C1 = glm(fl ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C2 = glm(fl ~ TSF, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C3 = glm(fl ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C4 = glm(fl ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C5 = glm(fl ~ TSF + density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C6 = glm(fl ~ TSF + density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C7 = glm(fl ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C8 = glm(fl ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)
floweringprobSR.C8 = glm(fl ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one", ], family = binomial)


AICctab(floweringprobSR.C1, floweringprobSR.C2, floweringprobSR.C3, floweringprobSR.C4, floweringprobSR.C5, floweringprobSR.C6, floweringprobSR.C7, floweringprobSR.C8, base = T) # floweringprobSR.C2 is the best model.

summary(floweringprobSR.C2)
floweringprobSR.siteC = floweringprobSR.C2


# Plotting the observed data vs the predictions of the model

floweringprobSR.siteC.obs = aggregate(fl ~ TSF, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF != "one"), ], FUN = mean)

floweringprobSR.siteC.pred = data.frame(TSF = floweringprobSR.siteC.obs[, 1])
floweringprobSR.siteC.pred$pred = predict(floweringprobSR.C2, newdata = floweringprobSR.siteC.pred, type = "response")


# Obs vs. Pred

plot(floweringprobSR.siteC.obs$fl, floweringprobSR.siteC.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted transition probability", main = "Observed vs. Predicted values of SR flowering probability - Site C")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

floweringprobSR.siteC.new.data = expand.grid(TSF = c("two", "three"))


# 95% confidence intervals

floweringprobSR.siteC.predict = predict(floweringprobSR.siteC, floweringprobSR.siteC.new.data, se.fit = T)
floweringprobSR.siteC.new.data$pred = inv.logit(floweringprobSR.siteC.predict$fit)
floweringprobSR.siteC.new.data$lwr = inv.logit(floweringprobSR.siteC.predict$fit - (1.96 * floweringprobSR.siteC.predict$se.fit))
floweringprobSR.siteC.new.data$upr = inv.logit(floweringprobSR.siteC.predict$fit + (1.96 * floweringprobSR.siteC.predict$se.fit))


# Plotting the predictions

plot.floweringprobSR.siteC = ggplot(floweringprobSR.siteC.new.data, aes(x = TSF, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("Time Since Fire (TSF)") + 
 ylab("Flowering probability") + 
 ggtitle("Small adult flowering probability - Site C") + 
 scale_x_discrete(labels = c("TSF2", "TSF3")) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_SRFlowerProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobSR.siteC

dev.off()


## 2.5.2. Site C - Large reproductive ----
# -----------------------------------

table(data.droso$fl[which(data.droso$site == "siteC" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "LR")]) # There a few observations for TSF2. We can model TSF3 and use the observations for TSF2.

floweringprobLR.C1 = glm(fl ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
floweringprobLR.C2 = glm(fl ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
floweringprobLR.C3 = glm(fl ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)

AICctab(floweringprobLR.C1, floweringprobLR.C2, floweringprobLR.C3, base = T) # floweringprobLR.C1 is the best model.

summary(floweringprobLR.C1)
floweringprobLR.siteC.TSF3 = floweringprobLR.C1
floweringprobLR.siteC.TSF2 = mean(data.droso$fl[which(data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "two")])
floweringprobLR.siteC = aggregate(fl ~ TSF, data = data.droso[which(data.droso$site == "siteC" & data.droso$stage == "LR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

floweringprobLR.siteC.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

floweringprobLR.siteC.predict = predict(floweringprobLR.siteC.TSF3, floweringprobLR.siteC.new.data, se.fit = T)
floweringprobLR.siteC.new.data$pred = inv.logit(floweringprobLR.siteC.predict$fit)
floweringprobLR.siteC.new.data$lwr = inv.logit(floweringprobLR.siteC.predict$fit - (1.96 * floweringprobLR.siteC.predict$se.fit))
floweringprobLR.siteC.new.data$upr = inv.logit(floweringprobLR.siteC.predict$fit + (1.96 * floweringprobLR.siteC.predict$se.fit))


# Plotting the predictions

plot.floweringprobLR.siteC = ggplot(floweringprobLR.siteC.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Flowering probability") + 
 ggtitle("Large adult flowering probability - Site C, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_LRFlowerProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobLR.siteC

dev.off()


## 2.5.2. Site E - Small reproductive ----
# -----------------------------------

table(data.droso$fl[which(data.droso$site == "siteE" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "SR")]) # There a very few observations in TSF2. We can use a model for TSF3 and the observations for TSF2.

floweringprobSR.E1 = glm(fl ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)
floweringprobSR.E2 = glm(fl ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)
floweringprobSR.E3 = glm(fl ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = binomial)

AICctab(floweringprobSR.E1, floweringprobSR.E2, floweringprobSR.E3, base = T) # floweringprobSR.E1, floweringprobSR.E2, and floweringprobSR.E3 are in the 2 dAIC range. floweringprobSR.E1 is the simplest model.

summary(floweringprobSR.E1)
floweringprobSR.siteE.TSF3 = floweringprobSR.E1
floweringprobSR.siteE.TSF2 = mean(data.droso$fl[which(data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "two")], na.rm = T)
floweringprobSR.siteE = aggregate(fl ~ TSF, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "SR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

floweringprobSR.siteE.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

floweringprobSR.siteE.predict = predict(floweringprobSR.siteE.TSF3, floweringprobSR.siteE.new.data, se.fit = T)
floweringprobSR.siteE.new.data$pred = inv.logit(floweringprobSR.siteE.predict$fit)
floweringprobSR.siteE.new.data$lwr = inv.logit(floweringprobSR.siteE.predict$fit - (1.96 * floweringprobSR.siteE.predict$se.fit))
floweringprobSR.siteE.new.data$upr = inv.logit(floweringprobSR.siteE.predict$fit + (1.96 * floweringprobSR.siteE.predict$se.fit))


# Plotting the predictions

plot.floweringprobSR.siteE = ggplot(floweringprobSR.siteE.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Flowering probability") + 
 ggtitle("Small adult flowering probability - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_SRFlowerProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobSR.siteE

dev.off()


## 2.5.3 Site E - Large reproductive ----
# ----------------------------------

table(data.droso$fl[which(data.droso$site == "siteE" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "LR")]) # There are few observations in TSF2. We can use a model for TSF3 and the observations for TSF2.

floweringprobLR.E1 = glm(fl ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
floweringprobLR.E2 = glm(fl ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)
floweringprobLR.E3 = glm(fl ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = binomial)

AICctab(floweringprobLR.E1, floweringprobLR.E2, floweringprobLR.E3, base = T) # floweringprobLR.E1, floweringprobLR.E2, and floweringprobLR.E3 are in the 2 dAIC range. floweringprobLR.E1 is the simplest models.

summary(floweringprobLR.E1)
floweringprobLR.siteE.TSF3 = floweringprobLR.E1
floweringprobLR.siteE.TSF2 = mean(data.droso$fl[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "two")], FUN = mean)
floweringprobLR.siteE = aggregate(fl ~ TSF, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "LR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

floweringprobLR.siteE.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

floweringprobLR.siteE.predict = predict(floweringprobLR.siteE.TSF3, floweringprobLR.siteE.new.data, se.fit = T)
floweringprobLR.siteE.new.data$pred = inv.logit(floweringprobLR.siteE.predict$fit)
floweringprobLR.siteE.new.data$lwr = inv.logit(floweringprobLR.siteE.predict$fit - (1.96 * floweringprobLR.siteE.predict$se.fit))
floweringprobLR.siteE.new.data$upr = inv.logit(floweringprobLR.siteE.predict$fit + (1.96 * floweringprobLR.siteE.predict$se.fit))


# Plotting the predictions

plot.floweringprobLR.siteE = ggplot(floweringprobLR.siteE.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Flowering probability") + 
 ggtitle("Large adult flowering probability - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_LRFlowerProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobLR.siteE

dev.off()


## 2.5.4. Site F - Small reproductive ----
# -----------------------------------

table(data.droso$fl[which(data.droso$site == "siteF" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "SR")]) # There are few observations in all TSFs. We use the observations.

floweringprobSR.siteF = aggregate(fl ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "SR"), ], FUN = mean)


# Plotting the predictions

plot.floweringprobSR.siteF = ggplot(floweringprobSR.siteF, aes(x = TSF, y = fl)) + 
 geom_point(size = 3) + 
 xlab("Time Since Fire (TSF)") + 
 ylab("Flowering probability") + 
 ggtitle("Small adult flowering probability - Site F") + 
 scale_x_discrete(labels = c("TSF2", "TSF3")) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_SRFlowerProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobSR.siteF

dev.off()


## 2.5.5. Site F - Large reproductive ----
# -----------------------------------

table(data.droso$fl[which(data.droso$site == "siteF" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "LR")])

floweringprobLR.F1 = glm(fl ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F2 = glm(fl ~ TSF, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F3 = glm(fl ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F4 = glm(fl ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F5 = glm(fl ~ TSF + density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F6 = glm(fl ~ TSF + density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F7 = glm(fl ~ TSF + density + TSF:density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F8 = glm(fl ~ TSF + density + density2 + TSF:density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)
floweringprobLR.F9 = glm(fl ~ TSF + density + density2 + TSF:density + TSF:density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF != "one", ], family = binomial)

AICctab(floweringprobLR.F1, floweringprobLR.F2, floweringprobLR.F3, floweringprobLR.F4, floweringprobLR.F5, floweringprobLR.F6, floweringprobLR.F7, floweringprobLR.F8, floweringprobLR.F9, base = T) # floweringprobLR.F1, floweringprobLR.F2, floweringprobLR.F3, floweringprobLR.F4, floweringprobLR.F5, and floweringprobLR.F6 are in the 2 dAIC range. floweringprobLR.F1 is the simplest model.

summary(floweringprobLR.F1)
floweringprobLR.siteF = floweringprobLR.F1


# Computing the predictions of the model and the 95 % CI

floweringprobLR.siteF.new.data = data.frame(pred = NA, lwr = NA, upr = NA)


# 95% confidence intervals

floweringprobLR.siteF.predict = predict(floweringprobLR.siteF, floweringprobLR.siteF.new.data, se.fit = T)
floweringprobLR.siteF.new.data$pred = inv.logit(floweringprobLR.siteF.predict$fit)
floweringprobLR.siteF.new.data$lwr = inv.logit(floweringprobLR.siteF.predict$fit - (1.96 * floweringprobLR.siteF.predict$se.fit))
floweringprobLR.siteF.new.data$upr = inv.logit(floweringprobLR.siteF.predict$fit + (1.96 * floweringprobLR.siteF.predict$se.fit))


# Plotting the predictions

plot.floweringprobLR.siteF = ggplot(floweringprobLR.siteF.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Flowering probability") + 
 ggtitle("Large adult flowering probability - Site F") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_LRFlowerProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobLR.siteF

dev.off()


## 2.6. Number of flowering stalks (fs) ----
# -------------------------------------

## 2.6.1. TSF2 ----
# ------------

table(data.droso$fs, data.droso$TSF, data.droso$site) # For all sites (except 2 outlier individuals in site E), individuals only have 1 flowering stalk in TSF2

nbfsTSF2 = 1


## 2.6.2. TSF3 ----
# ------------

## 2.6.2.1. Site C - Small reproductive ----
# -------------------------------------

table(data.droso$fs[which(data.droso$site == "siteC" & data.droso$stage == "SR" & data.droso$TSF == "three")]) # Only 1 stalk for all individuals.

nbfsSR.siteC.TSF3 = 1


## 2.6.2.2. Site C - Large reproductive ----
# -------------------------------------

table(data.droso$fs[which(data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three")])

nbfsLR.C1 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfsLR.C2 = glm(fs ~ density, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfsLR.C3 = glm(fs ~ density + density2, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfsLR.C1, nbfsLR.C2, nbfsLR.C3, base = T) # nbfsLS.C1 is the best model.

summary(nbfsLR.C1) # Model underdispersed

nbfsLR.siteC.TSF3 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = quasipoisson)

summary(nbfsLR.siteC.TSF3)


# Computing the predictions of the model and the 95 % CI

nbfsLR.siteC.predict = predict(nbfsLR.siteC.TSF3, se.fit = T)
nbfsLR.siteC.new.data = data.frame(pred = unique(exp(nbfsLR.siteC.predict$fit)), lwr = unique(exp(nbfsLR.siteC.predict$fit - (1.96 * nbfsLR.siteC.predict$se.fit))), upr = unique(exp(nbfsLR.siteC.predict$fit + (1.96 * nbfsLR.siteC.predict$se.fit))))

nbfsLR.siteC.new.data = rbind(nbfsLR.siteC.new.data, c(1, NA, NA))
nbfsLR.siteC.new.data$TSF = c("TSF3", "TSF2")


# Plotting the predictions

plot.nbfsLR.siteC = ggplot(nbfsLR.siteC.new.data, aes(x = TSF, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("Time Since Fire (TSF)") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Large adult number of flowering stalks - Site C") + 
 ylim(c(0, 3)) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteC_LRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsLR.siteC

dev.off()


## 2.6.2.3. Site E - Small reproductive ----
# -------------------------------------

table(data.droso$fs[which(data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three")])

nbfsSR.E1 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = poisson)
nbfsSR.E2 = glm(fs ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = poisson)
nbfsSR.E3 = glm(fs ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfsSR.E1, nbfsSR.E2, nbfsSR.E3, base = T) # nbfsLS.E1 and nbfsLR.E2 are in the 2 dAIC range. nbfsLR.E1 is the simplest model, but we can also check the other to see whether it brings something to the predictions.

summary(nbfsSR.E1) # Underdispersed

nbfsSR.siteE.TSF3 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = quasipoisson)

summary(nbfsSR.siteE.TSF3)


# Computing the predictions of the model and the 95 % CI

nbfsSR.siteE.predict = predict(nbfsSR.siteE.TSF3, se.fit = T)
nbfsSR.siteE.new.data = data.frame(pred = unique(exp(nbfsSR.siteE.predict$fit)), 
                                   lwr = unique(exp(nbfsSR.siteE.predict$fit - (1.96 * nbfsSR.siteE.predict$se.fit))), 
                                   upr = unique(exp(nbfsSR.siteE.predict$fit + (1.96 * nbfsSR.siteE.predict$se.fit))))

nbfsSR.siteE.new.data = rbind(nbfsSR.siteE.new.data, c(1, NA, NA))
nbfsSR.siteE.new.data$TSF = c("TSF3", "TSF2")


# Plotting the predictions

plot.nbfsSR.siteE = ggplot(nbfsSR.siteE.new.data, aes(x = TSF, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("Time Since Fire (TSF)") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Small adult number of flowering stalks - Site E") + 
 ylim(c(0, 2)) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_SRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsSR.siteE

dev.off()


## 2.6.2.4. Site E - Large reproductive ----
# -------------------------------------

table(data.droso$fs[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three")])

nbfsLR.E1 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfsLR.E2 = glm(fs ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfsLR.E3 = glm(fs ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfsLR.E1, nbfsLR.E2, nbfsLR.E3, base = T) # nbfsLR.E1, nbfsLR.E2, and nbfsLR.E3 are in the 2 dAIC range. nbfsLR.E1 is the simplest model.

summary(nbfsLR.E1) # Underdispersed

nbfsLR.siteE.TSF3 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = quasipoisson)

summary(nbfsLR.siteE.TSF3)


# Computing the predictions of the model and the 95 % CI

nbfsLR.siteE.predict = predict(nbfsLR.siteE.TSF3, se.fit = T)
nbfsLR.siteE.new.data = data.frame(pred = unique(exp(nbfsLR.siteE.predict$fit)), 
                                   lwr = unique(exp(nbfsLR.siteE.predict$fit - (1.96 * nbfsLR.siteE.predict$se.fit))), 
                                   upr = unique(exp(nbfsLR.siteE.predict$fit + (1.96 * nbfsLR.siteE.predict$se.fit))))

nbfsLR.siteE.new.data = rbind(nbfsLR.siteE.new.data, c(1, NA, NA))
nbfsLR.siteE.new.data$TSF = c("TSF3", "TSF2")


# Plotting the predictions

plot.nbfsLR.siteE = ggplot(nbfsLR.siteE.new.data, aes(x = TSF, y = pred)) + 
   geom_point(size = 3) + 
   geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
   xlab("Time Since Fire (TSF)") + 
   ylab("Number of flowering stalks") + 
   ggtitle("Large adult number of flowering stalks - Site E") + 
   ylim(c(0, 3)) + 
   theme_bw() + 
   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_LRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsLR.siteE

dev.off()


## 2.6.2.5. Site F - Small reproductive ----
# -------------------------------------

table(data.droso$fs[which(data.droso$site == "siteF" & data.droso$stage == "SR" & data.droso$TSF == "three")]) # Only 1s

nbfsSR.siteF.TSF3 = 1


## 2.6.2.6. Site F - Large reproductive ----
# -------------------------------------

table(data.droso$fs[which(data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three")])

nbfsLR.F1 = glm(fs ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfsLR.F2 = glm(fs ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfsLR.F3 = glm(fs ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfsLR.F1, nbfsLR.F2, nbfsLR.F3, base = T) # nbfsLS.F3 is the best model.

summary(nbfsLR.F3) # Overdispersed

nbfsLR.siteF.TSF3 = glm(fs ~ density + I(density^2), data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = quasipoisson)

summary(nbfsLR.siteF.TSF3)


# Plotting the observed data vs the predictions of the model

nbfsLR.siteF.obs = aggregate(fs ~ density, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three"), ], FUN = mean)

nbfsLR.siteF.pred = data.frame(density = nbfsLR.siteF.obs[, 1])
nbfsLR.siteF.pred$pred = predict(nbfsLR.siteC.TSF3, newdata = nbfsLR.siteF.pred, type = "response")


# Obs vs. Pred

plot(nbfsLR.siteF.obs$fs, nbfsLR.siteF.pred$pred, pch = 16, ylim = c(0, 6), xlim = c(0, 6), xlab = "Observed", ylab = "Predicted ", main = "Observed vs. Predicted values of LR number of flowering stalks - Site F")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

nbfsLR.siteF.new.data = expand.grid(density = seq(0, 60))


# 95% confidence intervals

nbfsLR.siteF.predict = predict(nbfsLR.siteF.TSF3, nbfsLR.siteF.new.data, se.fit = T)
nbfsLR.siteF.new.data$pred = exp(nbfsLR.siteF.predict$fit)
nbfsLR.siteF.new.data$lwr = exp(nbfsLR.siteF.predict$fit - (1.96 * nbfsLR.siteF.predict$se.fit))
nbfsLR.siteF.new.data$upr = exp(nbfsLR.siteF.predict$fit + (1.96 * nbfsLR.siteF.predict$se.fit))


# Plotting the predictions

plot.nbfsLR.siteF = ggplot(nbfsLR.siteF.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Large adult number of flowering stalks - Site F, TSF3") + 
 ylim(c(0, 8)) + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_LRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsLR.siteF

dev.off()


## 2.7. Number of flowers per stalk (fps) ----
# ---------------------------------------

## 2.7.1. Site C - Small reproductive ----
# -----------------------------------

table(data.droso$fps[which(data.droso$site == "siteC" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "SR")]) # Few observations in all TSFs. We use the observations.

nbfpsSR.siteC = aggregate(fps ~ TSF, mean, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "SR", ], na.rm = T)


## 2.7.2. Site C - Large reproductive ----
# -----------------------------------

table(data.droso$fps[which(data.droso$site == "siteC" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteC" & data.droso$stage == "LR")]) # Few observations in all TSFs. We use the observations.

nbfpsLR.siteC = aggregate(fps ~ TSF, mean, data = data.droso[data.droso$site == "siteC" & data.droso$stage == "LR", ], na.rm = T)


## 2.7.3. Site E - Small reproductive ----
# -----------------------------------

table(data.droso$fps[which(data.droso$site == "siteE" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "SR")]) # Few observations in TSF2. We use a model for TSF3 and the observations for TSF2.

nbfpsSR.E1 = glm(fps ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = poisson)
nbfpsSR.E2 = glm(fps ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = poisson)
nbfpsSR.E3 = glm(fps ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfpsSR.E1, nbfpsSR.E2, nbfpsSR.E3, base = T) # nbfpsSR.E1 and nbfpsSR.E2 are in the 2 dAIC range. nbfpsSR.E1 is the simplest model.

summary(nbfpsSR.E1) # Underdispersed

nbfpsSR.siteE.TSF3 = glm(fps ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three", ], family = quasipoisson)

summary(nbfpsSR.siteE.TSF3)

nbfpsSR.siteE.TSF2 = mean(data.droso$fps[which(data.droso$site == "siteE" & data.droso$stage == "SR" & data.droso$TSF == "three")], na.rm = T)
nbfpsSR.siteE = aggregate(fps ~ TSF, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "SR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

nbfpsSR.siteE.new.data = data.frame(pred = NA, upr = NA, lwr = NA)


# 95% confidence intervals

nbfpsSR.siteE.predict = predict(nbfpsSR.siteE.TSF3, nbfpsSR.siteE.new.data, se.fit = T)
nbfpsSR.siteE.new.data$pred = exp(nbfpsSR.siteE.predict$fit)
nbfpsSR.siteE.new.data$lwr = exp(nbfpsSR.siteE.predict$fit - (1.96 * nbfpsSR.siteE.predict$se.fit))
nbfpsSR.siteE.new.data$upr = exp(nbfpsSR.siteE.predict$fit + (1.96 * nbfpsSR.siteE.predict$se.fit))


# Plotting the predictions

plot.nbfpsSR.siteE = ggplot(nbfpsSR.siteE.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Small adult number of flowers per stalk - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteE_SRNbFlowPerStalk.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsSR.siteE

dev.off()


## 2.7.4. Site E - Large reproductive ----
# -----------------------------------

table(data.droso$fps[which(data.droso$site == "siteE" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteE" & data.droso$stage == "LR")]) # Few observations in TSF2. We use a model for TSF3 and the observations for TSF2.

nbfpsLR.E1 = glm(fps ~ 1, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfpsLR.E2 = glm(fps ~ density, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfpsLR.E3 = glm(fps ~ density + density2, data = data.droso[data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfpsLR.E1, nbfpsLR.E2, nbfpsLR.E3, base = T) # nbfpsLR.E1 and nbfpsLR.E2 are in the 2 dAIC range. nbfpsLR.E1 is the simplest model.

summary(nbfpsLR.E1)

nbfpsLR.siteE.TSF3 = nbfpsLR.E1
nbfpsLR.siteE.TSF2 = mean(data.droso$fps[which(data.droso$site == "siteE" & data.droso$stage == "LR" & data.droso$TSF == "two")], na.rm = T)
nbfpsLR.siteE = aggregate(fps ~ TSF, data = data.droso[which(data.droso$site == "siteE" & data.droso$stage == "LR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

nbfpsLR.siteE.new.data = data.frame(pred = NA, upr = NA, lwr = NA)


# 95% confidence intervals

nbfpsLR.siteE.predict = predict(nbfpsLR.siteE.TSF3, nbfpsLR.siteE.new.data, se.fit = T)
nbfpsLR.siteE.new.data$pred = exp(nbfpsLR.siteE.predict$fit)
nbfpsLR.siteE.new.data$lwr = exp(nbfpsLR.siteE.predict$fit - (1.96 * nbfpsLR.siteE.predict$se.fit))
nbfpsLR.siteE.new.data$upr = exp(nbfpsLR.siteE.predict$fit + (1.96 * nbfpsLR.siteE.predict$se.fit))


# Plotting the predictions

plot.nbfpsLR.siteE = ggplot(nbfpsLR.siteE.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Large adult number of flowers per stalk - Site E, TSF3") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "siteE_LRNbFlowPerStalk.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsLR.siteE

dev.off()


## 2.7.5. Site F - Small reproductive ----
# -----------------------------------

table(data.droso$fps[which(data.droso$site == "siteF" & data.droso$stage == "SR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "SR")]) # Few observations in all TSFs. We use the observations.

nbfpsSR.siteF = aggregate(fps ~ TSF, mean, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "SR", ], na.rm = T)


## 2.7.6. Site F - Large reproductive ----
# -----------------------------------

table(data.droso$fps[which(data.droso$site == "siteF" & data.droso$stage == "LR")], data.droso$TSF[which(data.droso$site == "siteF" & data.droso$stage == "LR")]) # Few observations in TSF2. We use a model for TSF3 and the observations for TSF2.

nbfpsLR.F1 = glm(fps ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfpsLR.F2 = glm(fps ~ density, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)
nbfpsLR.F3 = glm(fps ~ density + density2, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = poisson)

AICctab(nbfpsLR.F1, nbfpsLR.F2, nbfpsLR.F3, base = T) # nbfpsLR.F1 is the best model.

summary(nbfpsLR.F1) # Underdispersed

nbfpsLR.siteF.TSF3 = glm(fps ~ 1, data = data.droso[data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "three", ], family = quasipoisson)

summary(nbfpsLR.siteF.TSF3)

nbfpsLR.siteF.TSF2 = mean(data.droso$fps[which(data.droso$site == "siteF" & data.droso$stage == "LR" & data.droso$TSF == "two")], na.rm = T)
nbfpsLR.siteF = aggregate(fps ~ TSF, data = data.droso[which(data.droso$site == "siteF" & data.droso$stage == "LR"), ], FUN = mean)


# Computing the predictions of the model and the 95 % CI

nbfpsLR.siteF.new.data = data.frame(pred = NA, upr = NA, lwr = NA)


# 95% confidence intervals

nbfpsLR.siteF.predict = predict(nbfpsLR.siteF.TSF3, nbfpsLR.siteF.new.data, se.fit = T)
nbfpsLR.siteF.new.data$pred = exp(nbfpsLR.siteF.predict$fit)
nbfpsLR.siteF.new.data$lwr = exp(nbfpsLR.siteF.predict$fit-(1.96*nbfpsLR.siteF.predict$se.fit))
nbfpsLR.siteF.new.data$upr = exp(nbfpsLR.siteF.predict$fit + (1.96*nbfpsLR.siteF.predict$se.fit))


# Plotting the predictions

plot.nbfpsLR.siteF = ggplot(nbfpsLR.siteF.new.data, aes(x = 1, y = pred)) + 
 geom_point(size = 3) + 
 geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
 xlab("") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Large adult number of flowers per stalk - Site F") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "SiteF_LRNbFlowPerStalk.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsLR.siteF

dev.off()




###########################################################################
#
# 3. Creating tables with mean densities and vital rates ----
#
###########################################################################

## 3.1. Computing per site yearly density ----
# ---------------------------------------

mean.density.per.square.siteC.TSF1 = density.init.siteC.TSF1
mean.density.per.square.siteE.TSF1 = density.init.siteE.TSF1
mean.density.per.square.siteF.TSF1 = density.init.siteF.TSF1

mean.density.per.square.siteC.TSF2 = density.init.siteC.TSF2
mean.density.per.square.siteE.TSF2 = density.init.siteE.TSF2
mean.density.per.square.siteF.TSF2 = density.init.siteF.TSF2

mean.density.per.square.siteC.TSF3 = density.init.siteC.TSF3
mean.density.per.square.siteE.TSF3 = density.init.siteE.TSF3
mean.density.per.square.siteF.TSF3 = density.init.siteF.TSF3

mean.density.per.site.per.TSF = expand.grid(site = unique(data.droso$site),
                                            TSF = unique(data.droso$TSF))
mean.density.per.site.per.TSF$mean.density = c(density.init.siteC.TSF1, 
                                               density.init.siteE.TSF1,
                                               density.init.siteF.TSF1,
                                               density.init.siteC.TSF2, 
                                               density.init.siteE.TSF2,
                                               density.init.siteF.TSF2,
                                               density.init.siteC.TSF3, 
                                               density.init.siteE.TSF3,
                                               density.init.siteF.TSF3)
   

## 3.2. Creating dataset of vital rates per site and TSF ----
# ------------------------------------------------------

vr.per.TSF.per.site = expand.grid(site = c("siteC", "siteE", "siteF"), 
                                  TSF = c("TSF1", "TSF2", "TSF3"))


## 3.2.1. Seedling survival ----
# -------------------------

vr.per.TSF.per.site$survSD = NA
vr.per.TSF.per.site$survSD[which(vr.per.TSF.per.site$site == "siteC")] = c(predict(survSD.siteC.TSF1, newdata = data.frame(density = mean.density.per.square.siteC.TSF1), type = "response"), survSD.siteC.TSF2, survSD.siteC.TSF3)

vr.per.TSF.per.site$survSD[which(vr.per.TSF.per.site$site == "siteE")] = c(predict(survSD.siteE.TSF1_3, newdata = data.frame(TSF = c("one"), density = mean.density.per.square.siteE.TSF1), type = "response"), survSD.siteE.TSF2, predict(survSD.siteE.TSF1_3, newdata = data.frame(TSF = c("three"), density = mean.density.per.square.siteE.TSF3), type = "response"))

vr.per.TSF.per.site$survSD[which(vr.per.TSF.per.site$site == "siteF")] = c(survSD.siteF.TSF1, survSD.siteF.TSF2, predict(survSD.siteF.TSF3, newdata = data.frame(density = mean.density.per.square.siteF.TSF3), type = "response"))


## 3.2.2. Juvenile survival ----
# -------------------------

vr.per.TSF.per.site$survJ = NA
vr.per.TSF.per.site$survJ[which(vr.per.TSF.per.site$site == "siteC")] = c(unique(predict(survJ.siteC.TSF1, type = "response")), survJ.siteC.TSF2, survJ.siteC.TSF3)

vr.per.TSF.per.site$survJ[which(vr.per.TSF.per.site$site == "siteE")] = predict(survJ.siteE, newdata = data.frame(TSF = c("one", "two", "three"), density = c(mean.density.per.square.siteE.TSF1, mean.density.per.square.siteE.TSF2, mean.density.per.square.siteE.TSF3)), type = "response")

vr.per.TSF.per.site$survJ[which(vr.per.TSF.per.site$site == "siteF")] = c(unique(predict(survJ.siteF.TSF1_3, type = "response")), survJ.siteF.TSF2, unique(predict(survJ.siteF.TSF1_3, type = "response")))


## 3.2.3. Small reproductive survival ----
# -----------------------------------

vr.per.TSF.per.site$survSR = NA
vr.per.TSF.per.site$survSR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, predict(survSR.siteC, newdata = data.frame(TSF = c("two", "three"), density = c(mean.density.per.square.siteC.TSF2, mean.density.per.square.siteC.TSF3)), type = "response"))

vr.per.TSF.per.site$survSR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, survSR.siteE.TSF2, predict(survSR.siteE.TSF3, newdata = data.frame(density = mean.density.per.square.siteE.TSF3), type = "response"))

vr.per.TSF.per.site$survSR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, survSR.siteF$surv)


## 3.2.4. Large reproductive survival ----
# -----------------------------------

vr.per.TSF.per.site$survLR = NA
vr.per.TSF.per.site$survLR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, survLR.siteC.TSF2, unique(predict(survLR.siteC.TSF3, type = "response")))

vr.per.TSF.per.site$survLR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, survLR.siteE.TSF2, predict(survLR.siteE.TSF3, newdata = data.frame(density = mean.density.per.square.siteE.TSF3), type = "response"))

vr.per.TSF.per.site$survLR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, predict(survLR.siteF, newdata = data.frame(TSF = c("two", "three")), type = "response"))


## 3.2.5. J-LR transition ----
# -----------------------

vr.per.TSF.per.site$transJ = NA
vr.per.TSF.per.site$transJ[which(vr.per.TSF.per.site$site == "siteC")] = c(unique(predict(transJ.siteC.TSF1, type = "response")), transJ.siteC.TSF2, transJ.siteC.TSF3)

vr.per.TSF.per.site$transJ[which(vr.per.TSF.per.site$site == "siteE")] = predict(transJ.siteE, newdata = data.frame(TSF = c("one", "two", "three"), density = c(mean.density.per.square.siteC.TSF1, mean.density.per.square.siteE.TSF2, mean.density.per.square.siteE.TSF3)), type = "response")

vr.per.TSF.per.site$transJ[which(vr.per.TSF.per.site$site == "siteF")] = c(unique(predict(transJ.siteF.TSF1, type = "response")), transJ.siteF.TSF2, transJ.siteF.TSF3)


## 3.2.6. SR-LR transition ----
# ------------------------

vr.per.TSF.per.site$transSR = NA
vr.per.TSF.per.site$transSR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, transSR.siteC$transitionSR)

vr.per.TSF.per.site$transSR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, transSR.siteE.TSF2, unique(predict(transSR.siteE.TSF3, type = "response")))

vr.per.TSF.per.site$transSR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, transSR.siteF$transitionSR)


## 3.2.7. LR-SR transition ----
# ------------------------

vr.per.TSF.per.site$transLR = NA
vr.per.TSF.per.site$transLR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, transLR.siteC.TSF2, predict(transLR.siteC.TSF3, newdata = data.frame(density = c(mean.density.per.square.siteC.TSF3)), type = "response"))

vr.per.TSF.per.site$transLR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, transLR.siteE.TSF2, predict(transLR.siteE.TSF3, newdata = data.frame(density = mean.density.per.square.siteE.TSF3), type = "response"))

vr.per.TSF.per.site$transLR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, unique(predict(transLR.siteF.TSF2, type = "response")), transLR.siteF.TSF3)


## 3.2.8. SR probability of flowering ----
# -----------------------------------

vr.per.TSF.per.site$floweringprobSR = NA
vr.per.TSF.per.site$floweringprobSR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, predict(floweringprobSR.siteC, newdata = data.frame(TSF = c("two", "three")), type = "response"))

vr.per.TSF.per.site$floweringprobSR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, floweringprobSR.siteE.TSF2, unique(predict(floweringprobSR.siteE.TSF3, type = "response")))

vr.per.TSF.per.site$floweringprobSR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, floweringprobSR.siteF$fl)


## 3.2.9. LR probability of flowering ----
# -----------------------------------

vr.per.TSF.per.site$floweringprobLR = NA
vr.per.TSF.per.site$floweringprobLR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, floweringprobLR.siteC.TSF2, unique(predict(floweringprobLR.siteC.TSF3, type = "response")))

vr.per.TSF.per.site$floweringprobLR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, floweringprobLR.siteE.TSF2, unique(predict(floweringprobLR.siteE.TSF3, type = "response")))

vr.per.TSF.per.site$floweringprobLR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, rep(unique(predict(floweringprobLR.siteF, type = "response")), 2))


## 3.2.10. SR number of flowering stalks ----
# --------------------------------------

vr.per.TSF.per.site$nbfsSR = NA
vr.per.TSF.per.site$nbfsSR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, nbfsTSF2, nbfsSR.siteC.TSF3)

vr.per.TSF.per.site$nbfsSR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, nbfsTSF2, unique(predict(nbfsSR.siteE.TSF3, type = "response")))

vr.per.TSF.per.site$nbfsSR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, nbfsTSF2, nbfsSR.siteF.TSF3)


## 3.2.11. LR number of flowering stalks ----
# --------------------------------------

vr.per.TSF.per.site$nbfsLR = NA
vr.per.TSF.per.site$nbfsLR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, nbfsTSF2, unique(predict(nbfsLR.siteC.TSF3, type = "response")))

vr.per.TSF.per.site$nbfsLR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, nbfsTSF2, unique(predict(nbfsLR.siteE.TSF3)))

vr.per.TSF.per.site$nbfsLR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, nbfsTSF2, unique(predict(nbfsLR.siteF.TSF3, newdata = data.frame(density = mean.density.per.square.siteE.TSF3), type = "response")))


## 3.2.12 SR number of flowers per stalk ----
# --------------------------------------

vr.per.TSF.per.site$nbfpsSR = NA
vr.per.TSF.per.site$nbfpsSR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, nbfpsSR.siteC$fps)

vr.per.TSF.per.site$nbfpsSR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, nbfpsSR.siteE.TSF2, unique(predict(nbfpsSR.siteE.TSF3, type = "response")))

vr.per.TSF.per.site$nbfpsSR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, nbfpsSR.siteF$fps)


## 3.2.13. LR number of flowers per stalk ----
# ---------------------------------------

vr.per.TSF.per.site$nbfpsLR = NA
vr.per.TSF.per.site$nbfpsLR[which(vr.per.TSF.per.site$site == "siteC")] = c(NA, nbfpsLR.siteC$fps)

vr.per.TSF.per.site$nbfpsLR[which(vr.per.TSF.per.site$site == "siteE")] = c(NA, nbfpsLR.siteE.TSF2, unique(predict(nbfpsLR.siteE.TSF3, type = "response")))

vr.per.TSF.per.site$nbfpsLR[which(vr.per.TSF.per.site$site == "siteF")] = c(NA, nbfpsLR.siteF.TSF2, unique(predict(nbfpsLR.siteF.TSF3, type = "response")))




###########################################################################
#
# 4. Saving models and tables ----
#
###########################################################################

## 4.1. Saving updated dataset, number of squares, vital rates, and density table ----
# -------------------------------------------------------------------------------

write.csv(data.droso, "../RData/dataDroso_2019_formatted.csv", row.names = F)
write.csv(nbsquares.per.site.per.TSF, "NbSquares_PerSites_TSF1_3.csv", row.names = F)
write.csv(vr.per.TSF.per.site, "DPVR_TSF1_to_3.csv", row.names = F)
write.csv(mean.density.per.site.per.TSF, "MeanDensities.csv", row.names = F)


## 4.2. Saving models ----
# -------------------

## 4.2.1. Seedling survival ----
# -------------------------

save(survSD.siteC.TSF1, file = "GLM_survSD_siteC_TSF1.RData")

save(survSD.siteE.TSF1_3, file = "GLM_survSD_siteE_TSF1_3.RData")

save(survSD.siteF.TSF3, file = "GLM_survSD_siteF_TSF3.RData")


## 4.2.2. Juvenile survival ----
# -------------------------

save(survJ.siteC.TSF1, file = "GLM_survJ_siteC_TSF1.RData")

save(survJ.siteE, file = "GLM_survJ_siteE.RData")

save(survJ.siteF.TSF1_3, file = "GLM_survJ_siteF_TSF1_3.RData")


## 4.2.3. Small reproductive survival ----
# -----------------------------------

save(survSR.siteC, file = "GLM_survSR_siteC.RData")

save(survSR.siteE.TSF3, file = "GLM_survSR_siteE_TSF3.RData")


## 4.2.4. Large reproductive survival ----
# -----------------------------------

save(survLR.siteC.TSF3, file = "GLM_survLR_siteC_TSF3.RData")

save(survLR.siteE.TSF3, file = "GLM_survLR_siteE_TSF3.RData")

save(survLR.siteF, file = "GLM_survLR_siteF.RData")


## 4.2.5. J-LR transition ----
# -----------------------

save(transJ.siteC.TSF1, file = "GLM_transJ_siteC_TSF1.RData")

save(transJ.siteE, file = "GLM_transJ_siteE.RData")

save(transJ.siteF.TSF1, file = "GLM_transJ_siteF_TSF1.RData")


## 4.2.6. SR-LR transition ----
# ------------------------

save(transSR.siteE.TSF3, file = "GLM_transSR_siteE_TSF3.RData")


## 4.2.7. LR-SR transition ----
# ------------------------

save(transLR.siteC.TSF3, file = "GLM_transLR_siteC_TSF3.RData")

save(transLR.siteE.TSF3, file = "GLM_transLR_siteE_TSF3.RData")

save(transLR.siteF.TSF2, file = "GLM_transLR_siteF_TSF2.RData")


## 4.2.8. SR flowering probability ----
# --------------------------------

save(floweringprobSR.siteC, file = "GLM_floweringprobSR_siteC.RData")

save(floweringprobSR.siteE.TSF3, file = "GLM_floweringprobSR_siteE_TSF3.RData")

## 4.2.9. LR flowering probability ----
# --------------------------------

save(floweringprobLR.siteC.TSF3, file = "GLM_floweringprobLR_siteC_TSF3.RData")

save(floweringprobLR.siteE.TSF3, file = "GLM_floweringprobLR_siteE_TSF3.RData")

save(floweringprobLR.siteF, file = "GLM_floweringprobLR_siteF.RData")


## 4.2.10. SR number of flowering stalks  ----
# -------------------------------------

save(nbfsSR.siteE.TSF3, file = "GLM_nbfsSR_siteE_TSF3.RData")


## 4.2.11. LR number of flowering stalks ----
# --------------------------------------

save(nbfsLR.siteC.TSF3, file = "GLM_nbfsLR_siteC_TSF3.RData")

save(nbfsLR.siteE.TSF3, file = "GLM_nbfsLR_siteE_TSF3.RData")

save(nbfsLR.siteF.TSF3, file = "GLM_nbfsLR_siteF_TSF3.RData")


## 4.2.12. SR number of flowers per stalk ----
# ---------------------------------------

save(nbfpsSR.siteE.TSF3, file = "GLM_nbfpsSR_siteE_TSF3.RData")


## 4.2.13. LR number of flowers per stalk ----
# ---------------------------------------

save(nbfpsLR.siteE.TSF3, file = "GLM_nbfpsLR_siteE_TSF3.RData")

save(nbfpsLR.siteF.TSF3, file = "GLM_nbfpsLR_siteF_TSF3.RData")

