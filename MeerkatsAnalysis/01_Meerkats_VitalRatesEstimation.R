############################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script uses the capture-recapture data of a meerkat population in the Kalahari desert, collected between 1998 and 2016. 
# The aim of this script is to model the survival, transitions, emigration, and recruitment of four life-history stages: juvenile, subadult, helper, and dominant. We model as well the population range. 
#
# Author: Eva Conquet
#
###########################################################################


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
  library(lubridate)
  library(lme4)
  library(ggplot2)
  library(MASS)
  library(nlme)
  library(bbmle)
  library(optimx)
  library(multcomp)
  library(MuMIn)
}

load.librairies()


## 1.3. Loading and prepating data ----
# --------------------------------

data.meerkats = read.csv("MeerkatsData.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Creating subsets for the different stages

juveniles = data.meerkats[data.meerkats$stage == "J", ]
subadults = data.meerkats[data.meerkats$stage == "S", ]
helpers   = data.meerkats[data.meerkats$stage == "H", ]
dominants = data.meerkats[data.meerkats$stage == "D", ]




###########################################################################
#
# 2. Fitting Generalized Linear Mixed Models (GLMMs) 
# to estimate the vital rates ----
#
###########################################################################

## 2.1. Survival (surv) ----
# ---------------------

## 2.1.1. Juveniles (J) ----
# ---------------------

# Finding best random effect

surv.J1 = glmer(surv ~ season + (1|year), data = juveniles, family = binomial)
surv.J2 = glmer(surv ~ season + (1 + season|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.J1, null = glmer(surv ~ 1 + (1|year), data = juveniles, family = binomial)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept only.
r.squaredGLMM(surv.J2, null = glmer(surv ~ 1 + (1 + season|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))) 


# Finding best fixed effect

surv.J4 = glmer(surv ~ 1 + (1|year), data = juveniles, family = binomial)

surv.J5 = glmer(surv ~ season + (1|year), data = juveniles, family = binomial)
surv.J6 = glmer(surv ~ density + (1|year), data = juveniles, family = binomial)

surv.J7 = glmer(surv ~ season + density + (1|year), data = juveniles, family = binomial)
surv.J8 = glmer(surv ~ season + density + season:density + (1|year), data = juveniles, family = binomial)

surv.J9 = glmer(surv ~ season + density + density2 + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.J10 = glmer(surv ~ season + density + density2 + season:density + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.J11 = glmer(surv ~ season + density + density2 + season:density + season:density2 + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.J12 = glmer(surv ~ density + density2 + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(surv.J4, surv.J5, surv.J6, surv.J7, surv.J8, surv.J9, surv.J10, surv.J11, surv.J12, base = T) # The best model is surv.J11. There is an effect of season, density and squared density, as well as an effect of the interaction between both season and density, and season and squared density.

summary(surv.J11) 

survJ = surv.J11
survJ.simple = glmer(surv ~ season + density + I(density^2) + season:density + season:I(density^2) + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Plotting the observed data vs the predictions of the model

survJ.obs.rain = aggregate(surv ~ year + density, data = data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "rain", ], FUN = mean)
survJ.obs.dry = aggregate(surv ~ year + density, data = data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "dry", ], FUN = mean)

  
  # Rain

survJ.pred.rain = data.frame(density = seq(min(survJ.obs.rain$density), max(survJ.obs.rain$density), length.out = length(survJ.obs.rain$density)), 
                             pred = predict(survJ.simple, newdata = data.frame(density = seq(min(survJ.obs.rain$density), max(survJ.obs.rain$density), length.out = length(survJ.obs.rain$density)), season = "rain"), type = "response", re.form = NA))


  # Obs vs. Pred

plot(survJ.obs.rain$surv, survJ.pred.rain$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of juveniles (J) survival - Rain season")
abline(a = 0, b = 1, col = "red")


  # Dry

survJ.pred.dry = data.frame(density = seq(min(survJ.obs.dry$density), max(survJ.obs.dry$density), length.out = length(survJ.obs.dry$density)), 
                            pred = predict(survJ.simple, newdata = data.frame(density = seq(min(survJ.obs.dry$density), max(survJ.obs.dry$density), length.out = length(survJ.obs.dry$density)), season = "dry"), type = "response", re.form = NA))


  # Obs vs. Pred

plot(survJ.obs.dry$surv, survJ.pred.dry$pred, 
     pch = 16, xlim = c(0.6, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of juveniles (J) survival - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survJ.new.data = expand.grid(year = unique(data.meerkats$year),
                             season = unique(data.meerkats$season),
                             density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = 100))


# 95% confidence intervals

bootstrap.survJ <- bootMer(survJ.simple, FUN = function(x) predict(x, survJ.new.data, re.form = NA),
                           nsim = 100)
survJ.new.data$lwr = inv.logit(apply(bootstrap.survJ$t, 2, quantile, 0.025, na.rm = T))
survJ.new.data$upr = inv.logit(apply(bootstrap.survJ $t, 2, quantile, 0.975, na.rm = T))
survJ.new.data$pred = predict(survJ.simple, newdata = survJ.new.data, type = "response", re.form = NA)


# Plotting the predictions 

plot.survJ = ggplot(survJ.new.data, aes(x = density, y = pred, group = season, colour = season)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), linetype = 0, alpha = 0.2) +
  xlab("Population density (indiv./km2)") +
  ylab("Survival probability") +
  ggtitle("Juvenile survival") +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Rain", "Dry"),
                    values = c(cbbPalette[2], cbbPalette[1])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Rain", "Dry"),
                      values = c(cbbPalette[2], cbbPalette[1])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "JuvenileSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survJ

dev.off()


## 2.1.2. Subadults (S) ----
# ---------------------

# Finding best random effect

surv.S1 = glmer(surv ~ season + (1|year), data = subadults, family = binomial)
surv.S2 = glmer(surv ~ season + (1 + season|year), data = subadults, family = binomial)


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.S1, null = glmer(surv ~ 1 + (1|year), data = subadults, family = binomial)) 
r.squaredGLMM(surv.S2, null = glmer(surv ~ 1 + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

surv.S4 = glmer(surv ~ 1 + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.S5 = glmer(surv ~ season + (1 + season|year), data = subadults, family = binomial)
surv.S6 = glmer(surv ~ density + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.S7 = glmer(surv ~ season + density + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.S8 = glmer(surv ~ season + density + season:density + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.S9 = glmer(surv ~ season + density + density2 + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.S10 = glmer(surv ~ season + density + density2 + season:density + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.S11 = glmer(surv ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.S12 = glmer(surv ~ density + density2 + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(surv.S4, surv.S5, surv.S6, surv.S7, surv.S8, surv.S9, surv.S10, surv.S11, surv.S12, base = T) # Best model is surv.S6. There is an effect of density.

summary(surv.S6) 

survS = surv.S6

# Plotting the observed data vs the predictions of the model

survS.obs.rain = aggregate(surv ~ year + density, data = data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "rain", ], FUN = mean)
survS.obs.dry = aggregate(surv ~ year + density, data = data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "dry", ], FUN = mean)


  # Rain

survS.pred.rain = data.frame(density = seq(min(survS.obs.rain$density), max(survS.obs.rain$density), length.out = length(survS.obs.rain$density)), 
                             pred = predict(survS, newdata = data.frame(density = seq(min(survS.obs.rain$density), max(survS.obs.rain$density), length.out = length(survS.obs.rain$density)), season = "rain"), type = "response", re.form = NA))


# Obs vs. Pred

plot(survS.obs.rain$surv, survS.pred.rain$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of subadults (S) survival - Rain season (w/o season)")
abline(a = 0, b = 1, col = "red")


  # Dry

survS.pred.dry = data.frame(density = seq(min(survS.obs.dry$density), max(survS.obs.dry$density), length.out = length(survS.obs.dry$density)),
                            pred = predict(survS, newdata = data.frame(density = seq(min(survS.obs.dry$density), max(survS.obs.dry$density), length.out = length(survS.obs.dry$density)), season = "dry"), type = "response", re.form = NA))


# Obs vs. Pred

plot(survS.obs.dry$surv, survS.pred.dry$pred, 
     pch = 16, xlim = c(0.6, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of subadults (S) survival - Dry season (w/o season)")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survS.new.data = expand.grid(year = levels(data.meerkats$year),
                             season = unique(data.meerkats$season),
                             density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = 100))


# 95% confidence intervals
bootstrap.survS <- bootMer(survS, FUN = function(x) predict(x, survS.new.data, re.form = NA),
                           nsim = 1000)
survS.new.data$lwr = inv.logit(apply(bootstrap.survS$t, 2, quantile, 0.025, na.rm = T))
survS.new.data$upr = inv.logit(apply(bootstrap.survS$t, 2, quantile, 0.975, na.rm = T))
survS.new.data$pred = predict(survS, newdata = survS.new.data, type = "response", re.form = NA)


# Plotting the predictions
 
plot.survS = ggplot(survS.new.data, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Population density (indiv./km2)") +
  ylab("Survival probability") +
  ggtitle("Juvenile survival") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "SubadultSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survS

dev.off()


## 2.1.3. Helpers (H) ----
# -------------------

# Finding best random effect

surv.H1 = glmer(surv ~ season + (1|year), data = helpers, family = binomial)
surv.H2 = glmer(surv ~ season + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.H1, null = glmer(surv ~ 1 + (1|year), data = helpers, family = binomial)) 
r.squaredGLMM(surv.H2, null = glmer(surv ~ 1 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

surv.H4 = glmer(surv ~ 1 + (1 + season|year), data = helpers, family = binomial)

surv.H5 = glmer(surv ~ season + (1 + season|year), data = helpers, family = binomial)
surv.H6 = glmer(surv ~ density + (1 + season|year), data = helpers, family = binomial)

surv.H7 = glmer(surv ~ season + density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.H8 = glmer(surv ~ season + density + season:density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.H9 = glmer(surv ~ season + density + density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.H10 = glmer(surv ~ season + density + density2 + season:density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.H11 = glmer(surv ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.H12 = glmer(surv ~ density + density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(surv.H4, surv.H5, surv.H6, surv.H7, surv.H8, surv.H9, surv.H10, surv.H11, surv.H12, base = T) # Best model is surv.H5. There is an effect of season.

summary(surv.H5) 

survH = surv.H5


# Plotting the observed data vs the predictions of the model

survH.obs.rain = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
survH.obs.dry = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


  # Rain

survH.pred.rain = data.frame(year = levels(data.meerkats$year), 
                             pred = predict(survH, newdata = data.frame(year = levels(data.meerkats$year), season = "rain"), type = "response"))


# Obs vs. Pred

plot(survH.obs.rain$surv, survH.pred.rain$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of helpers (H) survival - Rain season")
abline(a = 0, b = 1, col = "red")


  # Dry

survH.pred.dry = data.frame(year = levels(data.meerkats$year),
                            pred = predict(survH, newdata = data.frame(year = levels(data.meerkats$year), season = "dry"), type = "response"))


# Obs vs. Pred

plot(survH.obs.dry$surv, survH.pred.dry$pred,
     pch = 16, xlim = c(0.6, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of helpers (H) survival - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survH.new.data = expand.grid(year = levels(data.meerkats$year),
                             season = unique(data.meerkats$season))


# 95% confidence intervals

bootstrap.survH <- bootMer(survH, FUN = function(x) predict(x, survH.new.data, re.form = NA),
                           nsim = 100)
survH.new.data$lwr = inv.logit(apply(bootstrap.survH$t, 2, quantile, 0.025, na.rm = T))
survH.new.data$upr = inv.logit(apply(bootstrap.survH$t, 2, quantile, 0.975, na.rm = T))
survH.new.data$pred = predict(survH, newdata = survH.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survH = ggplot(survH.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 10, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, position = position_dodge(width = 0.4)) +
  xlab("Season") +
  ylab("Helper survival probability") +
  #ggtitle("Helper survival") +
  scale_x_discrete(labels = c("Dry", "Rain")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "HelperSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survH

dev.off()


## 2.1.4. Dominants (D) ----
# ---------------------

# Finding best random effect

surv.D1 = glmer(surv ~ season + (1|year), data = dominants, family = binomial)
surv.D2 = glmer(surv ~ season + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.D1, null = glmer(surv ~ 1 + (1|year), data = dominants, family = binomial)) 
r.squaredGLMM(surv.D2, null = glmer(surv ~ 1 + (1 + season|year), data = dominants, family = binomial,control = glmerControl(optimizer = "bobyqa", calc.derivs = F))) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

surv.D4 = glmer(surv ~ 1 + (1 + season|year), data = dominants, family = binomial)

surv.D5 = glmer(surv ~ season + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.D6 = glmer(surv ~ density + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.D7 = glmer(surv ~ season + density + (1 + season|year), data = dominants, family = binomial)
surv.D8 = glmer(surv ~ season + density + season:density + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.D9 = glmer(surv ~ season + density + density2 + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.D10 = glmer(surv ~ season + density + density2 + season:density + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
surv.D11 = glmer(surv ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = dominants, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

surv.D12 = glmer(surv ~ density + density2 + (1 + season|year), data = dominants, family = binomial,control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(surv.D4, surv.D5, surv.D6, surv.D7, surv.D8, surv.D9, surv.D10, surv.D11, surv.D12, base = T)

summary(surv.D4) 

survD = surv.D4 


# Plotting the observed data vs the predictions of the model

survD.obs.rain = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "rain", ], FUN = mean)
survD.obs.dry = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "dry", ], FUN = mean)


# Rain

survD.pred.rain = data.frame(year = levels(data.meerkats$year), 
                             pred = predict(survD, newdata = data.frame(year = levels(data.meerkats$year), season = "rain"), type = "response"))


# Obs vs. Pred

plot(survD.obs.rain$surv, survD.pred.rain$pred,
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of dominants (D) survival - Rain season")
abline(a = 0, b = 1, col = "red")


# Dry

survD.pred.dry = data.frame(year = levels(data.meerkats$year),
                            pred = predict(survD, newdata = data.frame(year = levels(data.meerkats$year), season = "dry"), type = "response"))

# Obs vs. Pred

plot(survD.obs.dry$surv, survD.pred.dry$pred,
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of dominants (D) survival - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survD.new.data = expand.grid(year = levels(data.meerkats$year),
                             season = levels(data.meerkats$season))


# 95% confidence intervals

bootstrap.survD <- bootMer(survD, FUN = function(x) predict(x, survD.new.data),
                           nsim = 1000)
survD.new.data$lwr = inv.logit(apply(bootstrap.survD$t, 2, quantile, 0.025, na.rm = T))
survD.new.data$upr = inv.logit(apply(bootstrap.survD$t, 2, quantile, 0.975, na.rm = T))
survD.new.data$pred = predict(survD, newdata = survD.new.data, type = "response")


# Plotting the predictions

plot.survD = ggplot(survD.new.data, aes(x = as.numeric(year), y = pred, group = season, colour = season)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Survival probability") +
  ggtitle("Dominant survival") +
  scale_x_discrete(labels = as.numeric(survD.new.data$year)) +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Rain", "Dry"),
                    values = c(cbbPalette[2], cbbPalette[1])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Rain", "Dry"),
                      values = c(cbbPalette[2], cbbPalette[1])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "DominantSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survD

dev.off()



## 2.2. Emigration (emig) ----
# -----------------------

# Finding best random effect

emig1 = glmer(emig ~ season + (1|year), data = helpers, family = binomial)
emig2 = glmer(emig~ season + (1 + season|year), data = helpers, family = binomial)


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(emig1, null = glmer(emig ~ 1 + (1|year), data = helpers, family = binomial)) 
r.squaredGLMM(emig2, null = glmer(emig ~ 1 + (1 + season|year), data = helpers, family = binomial)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

emig4 = glmer(emig ~ 1 + (1 + season|year), data = helpers, family = binomial)

emig5 = glmer(emig ~ season + (1 + season|year), data = helpers, family = binomial)
emig6 = glmer(emig ~ density + (1 + season|year), data = helpers, family = binomial)

emig7 = glmer(emig ~ season + density + (1 + season|year), data = helpers, family = binomial)
emig8 = glmer(emig ~ season + density + season:density + (1 + season|year), data = helpers, family = binomial)

emig9 = glmer(emig ~ season + density + density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
emig10 = glmer(emig ~ season + density + density2 + season:density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
emig11 = glmer(emig ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

emig12 = glmer(emig ~ density + density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(emig4, emig5, emig6, emig7, emig8, emig9, emig10, emig11, emig12, base = T) # The best model is emig11. There is an effect of season, density and squared density, as well as an effect of the interaction between both season and density, and season and squared density.

summary(emig6)

emig = emig6


# Plotting the observed data vs the predictions of the model

emig.obs.rain = aggregate(emig ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
emig.obs.dry = aggregate(emig ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


# Rain

emig.pred.rain = expand.grid(season = "rain", 
                             year = levels(data.meerkats$year), 
                             density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = length(data.meerkats$density)))
emig.pred.rain$pred = predict(emig, newdata = emig.pred.rain, type = "response")
emig.pred.rain = aggregate(pred ~ year, data = emig.pred.rain, FUN = mean)


# Obs vs. Pred

plot(emig.obs.rain$emig, emig.pred.rain$pred,
     pch = 16, ylim = c(0, 1), xlab = "Observed emigration probability", ylab = "Predicted emigration probability",
     main = "Observed vs. Predicted values of helpers (H) emigration - Rain season")
abline(a = 0, b = 1, col = "red")


# Dry

emig.pred.dry = expand.grid(season = "dry",
                            year = levels(data.meerkats$year),
                            density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = length(data.meerkats$density)))
emig.pred.dry$pred = predict(emig, newdata = emig.pred.dry, type = "response")
emig.pred.dry = aggregate(pred ~ year, data = emig.pred.dry, FUN = mean)


# Obs vs. Pred

plot(emig.obs.dry$emig, emig.pred.dry$pred, 
     pch = 16, ylim = c(0, 1), xlab = "Observed emigration probability", ylab = "Predicted emigration probability",
     main="Observed vs. Predicted values of helpers (H) emigration - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

emig.new.data = expand.grid(year = levels(data.meerkats$year),
                            season = levels(data.meerkats$season),
                            density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = 100))


# 95% confidence intervals

bootstrap.emig <- bootMer(emig, FUN = function(x) predict(x, emig.new.data, re.form = NA),
                          nsim = 1000)
emig.new.data$lwr = inv.logit(apply(bootstrap.emig$t, 2, quantile, 0.025, na.rm = T))
emig.new.data$upr = inv.logit(apply(bootstrap.emig$t, 2, quantile, 0.975, na.rm = T))
emig.new.data$pred = predict(emig, newdata = emig.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.emig = ggplot(emig.new.data, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Population density (indiv./km2)") +
  ylab("Emigration probability") +
  ggtitle("Helper emigration") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "HelperEmigration.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.emig

dev.off()


## 2.3. Transition from helper to dominant (transition) ----
# -----------------------------------------------------

# Finding best random effect

transition1 = glmer(transition ~ season + (1|year), data = helpers, family = binomial)
transition2 = glmer(transition ~ season + (1 + season|year), data = helpers, family = binomial)


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(transition1, null = glmer(transition ~ 1 + (1|year), data = helpers, family = binomial)) 
r.squaredGLMM(transition2, null = glmer(transition ~ 1 + (1 + season|year), data = helpers, family = binomial)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

transition4 = glmer(transition ~ 1 + (1 + season|year), data = helpers, family = binomial)

transition5 = glmer(transition ~ season + (1 + season|year), data = helpers, family = binomial)
transition6 = glmer(transition ~ density  + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

transition7 = glmer(transition ~ season + density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
transition8 = glmer(transition ~ season + density + season:density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

transition9 = glmer(transition ~ season + density + density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
transition10 = glmer(transition ~ season + density + density2 + season:density + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
transition11 = glmer(transition ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

transition12 = glmer(transition ~ density + density2 + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(transition4, transition5, transition6, transition7, transition8, transition9, transition10, transition11, transition12, base = T) # The best model is transition11. There is an effect of season, density and squared density, as well as an effect of the interaction between both season and density, and season and squared density.

summary(transition6)

transition = transition6


# Plotting the observed data vs the predictions of the model

transition.obs.rain = aggregate(transition ~ year + density, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
transition.obs.dry = aggregate(transition ~ year + density, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


# Rain

transition.pred.rain = data.frame(density = seq(min(transition.obs.rain$density), max(transition.obs.rain$density), length.out = length(transition.obs.rain$density)),
                                  pred = predict(transition, newdata = data.frame(density = seq(min(transition.obs.rain$density), max(transition.obs.rain$density), length.out = length(transition.obs.rain$density)), season = "rain"), type = "response", re.form = NA))


# Obs vs. Pred

plot(transition.obs.rain$transition, transition.pred.rain$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed transition probability", ylab = "Predicted transition probability",
     main="Observed vs. Predicted values of H-D transition - Rain season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transition.new.data = expand.grid(year = levels(data.meerkats$year),
                                  season = levels(data.meerkats$season),
                                  density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = 100))


# 95% confidence intervals

bootstrap.transition <- bootMer(transition, FUN = function(x) predict(x, transition.new.data, re.form = NA),
                                nsim = 1000)
transition.new.data$lwr = inv.logit(apply(bootstrap.transition$t, 2, quantile, 0.025, na.rm = T))
transition.new.data$upr = inv.logit(apply(bootstrap.transition$t, 2, quantile, 0.975, na.rm = T))
transition.new.data$pred = predict(transition, newdata = transition.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.transition = ggplot(transition.new.data, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Population density (indiv./km2)") +
  ylab("Transition probability") +
  ggtitle("Helper-Dominant transition") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "HelperTransition.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.transition

dev.off()


## 2.4. Recruitment (pups) ----
# ------------------------

## 2.4.1. Helpers ----
# ---------------

# Finding best random effect

recruitment.H1 = glmer(pups ~ season + (1|year), data = helpers, family = poisson)
recruitment.H2 = glmer(pups~ season + (1 + season|year), data = helpers, family = poisson)


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(recruitment.H1, null = glmer(pups ~ 1 + (1|year), data = helpers, family = poisson)) 
r.squaredGLMM(recruitment.H2, null = glmer(pups ~ 1 + (1 + season|year), data = helpers, family = poisson)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

recruitment.H4 = glmer(pups ~ 1 + (1 + season|year), data = helpers, family = poisson)

recruitment.H5 = glmer(pups ~ season + (1 + season|year), data = helpers, family = poisson)
recruitment.H6 = glmer(pups ~ density + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.H7 = glmer(pups ~ season + density + (1 + season|year), data = helpers, family = poisson,control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
recruitment.H8 = glmer(pups ~ season + density + season:density + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.H9 = glmer(pups ~ season + density + density2 + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
recruitment.H10 = glmer(pups ~ season + density + density2 + season:density + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
recruitment.H11 = glmer(pups ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.H12 = glmer(pups ~ density + density2 + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(recruitment.H4, recruitment.H5, recruitment.H6, recruitment.H7, recruitment.H8, recruitment.H9, recruitment.H10, recruitment.H11, recruitment.H12, base = T) # The best model is recruitment.H11. There is an effect of season, density and squared density, as well as an effect of the interaction between both season and density, and season and squared density.

summary(recruitment.H6)

recruitment.H = recruitment.H6

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(recruitment.H)

recruitmentH.quasiP = MASS::glmmPQL(pups ~ density, random =  ~ 1+season|year, family = poisson, data = helpers)

summary(recruitmentH.quasiP)

recruitment.H = recruitmentH.quasiP


# Plotting the observed data vs the predictions of the model

recruitment.H.obs.rain = aggregate(pups ~ year + density, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
recruitment.H.obs.dry = aggregate(pups ~ year + density, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


# Rain  

recruitment.H.pred.rain = data.frame(density = seq(min(recruitment.H.obs.rain$density), max(recruitment.H.obs.rain$density), length.out = length(recruitment.H.obs.rain$density)),
                                     pred = predict(recruitment.H, newdata = expand.grid(density = seq(min(recruitment.H.obs.rain$density), max(recruitment.H.obs.rain$density), length.out = length(recruitment.H.obs.rain$density)), season = unique(data.meerkats$season)), type = "response", re.form = NA, level = 0))


# Obs vs. Pred

plot(recruitment.H.obs.rain$pups, recruitment.H.pred.rain$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of helper (H) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")


# Dry

recruitment.H.pred.dry = data.frame(density = seq(min(recruitment.H.obs.dry$density), max(recruitment.H.obs.dry$density), length.out = length(recruitment.H.obs.dry$density)),
                                    pred = predict(recruitment.H, newdata = data.frame(density = seq(min(recruitment.H.obs.dry$density), max(recruitment.H.obs.dry$density), length.out = length(recruitment.H.obs.dry$density)), season = "dry"), type = "response", re.form = NA))


# Obs vs. Pred

plot(recruitment.H.obs.dry$pups, recruitment.H.pred.dry$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of helper (H) recruitment - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

recruitment.H.new.data = expand.grid(year = levels(data.meerkats$year),
                                     season = levels(data.meerkats$season),
                                     density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = 100))


# 95% confidence intervals

bootstrap.recruitment.H <- bootMer(recruitment.H, FUN = function(x) predict(x, recruitment.H.new.data, re.form = NA), nsim = 100)
recruitment.H.new.data$lwr = exp(apply(bootstrap.recruitment.H$t, 2, quantile, 0.025, na.rm = T))
recruitment.H.new.data$upr = exp(apply(bootstrap.recruitment.H$t, 2, quantile, 0.975, na.rm = T))
recruitment.H.new.data$pred = predict(recruitment.H, newdata = recruitment.H.new.data, type = "response", re.form = NA)

plot.recruitment.H = ggplot(recruitment.H.new.data, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Population density (indiv./km2)") +
  ylab("Number of pups") +
  ggtitle("Helper recruitment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l =0 )), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "HelperRecruitment.H.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.recruitment.H

dev.off()


## 2.4.2. Dominants ----
# -----------------

# Finding best random effect

recruitment.D1 = glmer(pups ~ season + (1|year), data = dominants, family = poisson)
recruitment.D2 = glmer(pups ~ season + (1 + season|year), data = dominants, family = poisson)


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(recruitment.D1, null = glmer(pups ~ 1 + (1|year), data = dominants, family = poisson)) 
r.squaredGLMM(recruitment.D2, null = glmer(pups ~ 1 + (1 + season|year), data = dominants, family = poisson)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

recruitment.D4 = glmer(pups ~ 1 + (1 + season|year), data = dominants, family = poisson)

recruitment.D5 = glmer(pups ~ season + (1 + season|year), data = dominants, family = poisson)
recruitment.D6 = glmer(pups ~ density + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.D7 = glmer(pups ~ season + density + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
recruitment.D8 = glmer(pups ~ season + density + season:density + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.D9 = glmer(pups ~ season + density + density2 + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
recruitment.D10 = glmer(pups ~  season + density + density2 + season:density + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))
recruitment.D11 = glmer(pups ~ season + density + density2 + season:density + season:density2 + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.D12 = glmer(pups ~ density + density2 + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(recruitment.D4, recruitment.D5, recruitment.D6, recruitment.D7, recruitment.D8, recruitment.D9,  recruitment.D10,recruitment.D11, recruitment.D12, base = T) # The best model is recruitment.D11. There is an effect of season, density and squared density, as well as an effect of the interaction between both season and density, and season and squared density.

summary(recruitment.D8)
summary(recruitment.D10)

recruitment.D = glmer(pups ~ season + density + I(density^2) + season:density + season:I(density^2) + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
recruitment.D.2 = glmer(pups ~ season + density + I(density^2) + season:density + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

overdisp_fun(recruitment.D)

recruitmentD.quasiP = MASS::glmmPQL(pups ~ season + density + I(density^2) + season:density + season:I(density^2), random =  ~ 1+season|year, family = poisson, data = dominants, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

summary(recruitmentD.quasiP)

recruitment.D = recruitmentD.quasiP


# Plotting the observed data vs the predictions of the model

recruitment.D.obs.rain = aggregate(pups ~ year + density, data = data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "rain", ], FUN = mean)
recruitment.D.obs.dry = aggregate(pups ~ year + density, data = data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "dry", ], FUN = mean)


  # Rain

recruitment.D.pred.rain = data.frame(density = seq(min(recruitment.D.obs.rain$density), 20, length.out = length(recruitment.D.obs.rain$density)),
                                     pred = predict(recruitment.D, newdata = data.frame(density = seq(min(recruitment.D.obs.rain$density), 20, length.out = length(recruitment.D.obs.rain$density)), season = "rain"), type = "response", re.form = NA))

recruitment.D.pred.rain2 = data.frame(density = seq(min(recruitment.D.obs.rain$density), 20, length.out = length(recruitment.D.obs.rain$density)),
                                      pred = predict(recruitment.D.2, newdata = data.frame(density = seq(min(recruitment.D.obs.rain$density), 20, length.out = length(recruitment.D.obs.rain$density)), season = "rain"), type = "response", re.form = NA))

recruitment.D.pred.rain3 = data.frame(density = seq(min(recruitment.D.obs.rain$density), 20, length.out = length(recruitment.D.obs.rain$density)),
                                      pred = predict(recruitment.D8, newdata = data.frame(density = seq(min(recruitment.D.obs.rain$density), 20, length.out = length(recruitment.D.obs.rain$density)), season = "rain"), type = "response", re.form = NA))


# Obs vs. Pred

par(mfrow = c(1, 2))
 
plot(recruitment.D.obs.rain$pups, recruitment.D.pred.rain$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")

plot(recruitment.D.obs.rain$pups, recruitment.D.pred.rain2$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")

plot(recruitment.D.obs.rain$pups, recruitment.D.pred.rain3$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")
layout(1)


  # Dry

recruitment.D.pred.dry = data.frame(density = seq(min(recruitment.D.obs.dry$density), 20, length.out = length(recruitment.D.obs.dry$density)),
                                    pred = predict(recruitment.D, newdata = data.frame(density = seq(min(recruitment.D.obs.dry$density), 20, length.out = length(recruitment.D.obs.dry$density)), season = "dry"), type = "response", re.form = NA))

recruitment.D.pred.dry2 = data.frame(density = seq(min(recruitment.D.obs.dry$density), 20, length.out = length(recruitment.D.obs.dry$density)),
                                     pred = predict(recruitment.D8, newdata = data.frame(density = seq(min(recruitment.D.obs.dry$density), 20, length.out = length(recruitment.D.obs.dry$density)), season = "dry"), type = "response", re.form = NA))

recruitment.D.pred.dry3 = data.frame(density = seq(min(recruitment.D.obs.dry$density), 20, length.out = length(recruitment.D.obs.dry$density)),
                                     pred = predict(recruitment.D.2, newdata = data.frame(density = seq(min(recruitment.D.obs.dry$density), 20, length.out = length(recruitment.D.obs.dry$density)), season = "dry"), type = "response", re.form = NA))


# Obs vs. Pred

par(mfrow = c(1, 3))

plot(recruitment.D.obs.dry$pups, recruitment.D.pred.dry$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Dry season")
abline(a = 0, b = 1, col = "red")

plot(recruitment.D.obs.dry$pups, recruitment.D.pred.dry2$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Dry season")
abline(a = 0, b = 1, col = "red")

plot(recruitment.D.obs.dry$pups, recruitment.D.pred.dry3$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Dry season")
abline(a = 0, b = 1, col = "red")
layout(1)


# The two best models according to the AICc are not biologically relevant (increasing recruitment with population density), with no inflection point. We compare the predictions of the two best models and the most biologically relevant.

# Best models (recruitment.D8)

# Computing the predictions of the model and the 95 % CI
recruitment.D8.new.data = expand.grid(year = levels(data.meerkats$year),
                                      season = levels(data.meerkats$season),
                                      density = seq(min(data.meerkats$density),15,length.out = 100))


# 95% confidence intervals

bootstrap.recruitment.D8 <- bootMer(recruitment.D8, FUN = function(x) predict(x, recruitment.D8.new.data, re.form = NA),
                                    nsim = 1000)
recruitment.D8.new.data$lwr = exp(apply(bootstrap.recruitment.D8$t, 2, quantile, 0.025, na.rm = T))
recruitment.D8.new.data$upr = exp(apply(bootstrap.recruitment.D8$t, 2, quantile ,0.975, na.rm = T))
recruitment.D8.new.data$pred = predict(recruitment.D8, newdata = recruitment.D8.new.data, type = "response", re.form = NA)

plot.recruitment.D8 = ggplot(recruitment.D8.new.data, aes(x = density, y = pred, group = season, colour = season)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), linetype = 0, alpha = 0.2) +
  xlab("") +
  ylab("Number of pups") +
  ggtitle("") +
  ylim(0, 10) +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[1], cbbPalette[2])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Wet", "Dry"),
                      values = c(cbbPalette[1], cbbPalette[2])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l =0 )), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(5, "lines"))


# Second best model (recruitment.D10)

# Computing the predictions of the model and the 95 % CI

recruitment.D10.new.data = expand.grid(year = levels(data.meerkats$year),
                                       season = levels(data.meerkats$season),
                                       density = seq(min(data.meerkats$density), 15, length.out = 100))


# 95% confidence intervals

bootstrap.recruitment.D10 <- bootMer(recruitment.D.2, FUN = function(x) predict(x, recruitment.D10.new.data, re.form = NA), 
                                     nsim = 1000)
recruitment.D10.new.data$lwr = exp(apply(bootstrap.recruitment.D10$t, 2, quantile, 0.025, na.rm = T))
recruitment.D10.new.data$upr = exp(apply(bootstrap.recruitment.D10$t, 2, quantile, 0.975, na.rm = T))
recruitment.D10.new.data$pred = predict(recruitment.D.2, newdata = recruitment.D10.new.data, type = "response", re.form = NA)

plot.recruitment.D10 = ggplot(recruitment.D10.new.data, aes(x = density, y = pred, group = season, colour = season)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), linetype = 0, alpha = 0.2) +
  xlab(expression(bold("Pop. density (indiv./"*km^2*")"))) +
  ylab("") +
  ggtitle("") +
  ylim(0, 10) +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[1], cbbPalette[2])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Wet", "Dry"),
                      values = c(cbbPalette[1], cbbPalette[2])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l =0 )), 
        axis.text.y = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(5, "lines"))


# Most biologically relevant model

# Computing the predictions of the model and the 95 % CI

recruitment.D.new.data = expand.grid(year = levels(data.meerkats$year),
                                     season = levels(data.meerkats$season),
                                     density = seq(min(data.meerkats$density), 15, length.out = 1000))


# 95% confidence intervals

bootstrap.recruitment.D <- bootMer(recruitment.D, FUN = function(x) predict(x, recruitment.D.new.data, re.form = NA),
                                   nsim = 1000)
recruitment.D.new.data$lwr = exp(apply(bootstrap.recruitment.D$t, 2, quantile, 0.025, na.rm = T))
recruitment.D.new.data$upr = exp(apply(bootstrap.recruitment.D$t, 2, quantile, 0.975, na.rm = T))
recruitment.D.new.data$pred = predict(recruitment.D, newdata = recruitment.D.new.data, type = "response", re.form = NA)


plot.recruitment.D.FigS4 = ggplot(recruitment.D.new.data, aes(x = density, y = pred, group = season, colour = season)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = season), linetype = 0, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  ylim(0, 10) +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[1], cbbPalette[2])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Wet", "Dry"),
                      values = c(cbbPalette[1], cbbPalette[2])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "none", 
        legend.key.height = unit(4, "lines"), 
        legend.key.width = unit(5, "lines"))


# Figure S4.1 

meerkats.fig.S4.1 = ggarrange(plot.recruitment.D8, 
                              plot.recruitment.D10, 
                              plot.recruitment.D.FigS4, 
                              labels = c("(A)", "(B)", "(C)"), ncol = 3, nrow = 1, 
                              font.label = list(size = 26), hjust = -0.1, 
                              common.legend = T, legend = "bottom")

png(filename = "Figure_S4_1.png",
    width = 5500,
    height = 2000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

meerkats.fig.S4.1

dev.off()


## 2.5. Population range (to get population density) ----
# --------------------------------------------------

pr = aggregate(home.range ~ season + year, mean, data = data.meerkats, na.rm = T)
table(pr$year) # Check that we have two seasons per year


# Computing the total number of females
pr$totF = aggregate(data = data.meerkats, ID ~ season + year, function(x) length(unique(x)))$ID


# Computing the number of dominants
dom = droplevels(data.meerkats[data.meerkats$stage == "D", ])
pr$domNum = aggregate(data = dom, ID ~ season + year, function(x) length(unique(x)))$ID


# Computing the average group size - Dividing the total number of individuals (2*number of females) by the number of dominants (1 dominant pair per group)
pr$avG = (2 * pr$totF) / pr$domNum


# We fit the log of population range as response variable - as we want to avoid ever getting negative values.
hist(pr$home.range)
hist(log(pr$home.range))

pr$home.range = log(pr$home.range)


# Finding the best fixed effect. We use a RE on the intercept only because we have only one observation per season-year combination.

pop.range.null = lmer(home.range ~ 1 + (1|year), data = pr)
pop.range.full = lmer(home.range ~ season + totF + I(totF^2) + domNum + I(domNum^2) + avG + I(avG^2) + season:totF + season:domNum + season:avG + (1|year), data = pr, na.action = "na.fail")

pop.range1 = lmer(home.range ~ domNum + (1|year), data = pr)

MuMIn::dredge(pop.range.full)


# We take the exponential of the prediction in order to have the real population range value and not only the log value, as our response variable is on the log scale.

plot(1:length(unique(pr$year)), exp(aggregate(home.range ~ year, data = pr, mean)[, 2]), type = "l")
lines(1:length(unique(pr$year)), exp(predict(pop.range.null, newdata = data.frame(year = unique(pr$year)))), col = "red")
lines(1:length(unique(pr$year)), exp(predict(pop.range1, newdata = data.frame(domNum = aggregate(domNum ~ year, data = pr, mean)[, 2], year = unique(pr$year)))), col = "blue")


# Plotting the observed data vs the predictions of the model

pr.obs = aggregate(home.range ~ year + domNum, data = pr, FUN = mean)

pr.pred = data.frame(year = pr.obs$year,
                     domNum = pr.obs$domNum,
                     pred = predict(pop.range1, newdata = pr.obs[, c(1, 2)], type = "response"))


# Obs vs. Pred

plot(pr.obs$home.range, pr.pred$pred, 
     pch = 16, xlab = "Observed population range", ylab = "Predicted population range",
     main = "Observed vs. Predicted values of population range - Rain season")
abline(a = 0, b = 1, col = "red")


# Checking for temporal autocorrelation.

pr.autocorr <- lme(home.range ~ domNum, data = pr, 
                   random = ~ 1|year,
                   correlation = corAR1(),
                   na.action = na.omit)
summary(pr.autocorr) # The phi value is 0, there is little evidence for autocorrelation.

pop.range = pop.range1
summary(pop.range)




###########################################################################
#
# 3. Saving the vital rates ----
#
###########################################################################

vr.per.year = data.frame(year = unique(data.meerkats$year))
mean.density = mean(data.meerkats$density)
mean.density.squared = mean.density^2


# Survival of juveniles

surv.J.per.year.dry = predict(survJ, newdata = data.frame(year = vr.per.year$year, density = mean.density, density2 = mean.density.squared, season = "dry"), type = "response")
surv.J.per.year.wet = predict(survJ, newdata = data.frame(year = vr.per.year$year, density = mean.density, density2 = mean.density.squared, season = "rain"), type = "response")
vr.per.year$survJ.dry = surv.J.per.year.dry
vr.per.year$survJ.wet = surv.J.per.year.wet


# Survival of subadults

surv.S.per.year.dry = predict(survS, newdata = data.frame(year = vr.per.year$year, density = mean.density, season = "dry"), type = "response")
surv.S.per.year.wet = predict(survS, newdata = data.frame(year = vr.per.year$year, density = mean.density, season = "rain"), type = "response")
vr.per.year$survS.dry = surv.S.per.year.dry
vr.per.year$survS.wet = surv.S.per.year.wet


# Survival of helpers

surv.H.per.year.dry = predict(survH, newdata = data.frame(year = vr.per.year$year, season = "dry"), type = "response")
surv.H.per.year.wet = predict(survH, newdata = data.frame(year = vr.per.year$year, season = "rain"), type = "response")
vr.per.year$survH.dry = surv.H.per.year.dry
vr.per.year$survH.wet = surv.H.per.year.wet


# Survival of dominants

surv.D.per.year.dry = predict(survD, newdata = data.frame(year = vr.per.year$year, season = "dry"), type = "response")
surv.D.per.year.wet = predict(survD, newdata = data.frame(year = vr.per.year$year, season = "rain"), type = "response")
vr.per.year$survD.dry = surv.D.per.year.dry
vr.per.year$survD.wet = surv.D.per.year.wet


# Transition H -> D

trans.H.D.per.year.dry = predict(transition, newdata = data.frame(year = vr.per.year$year, density = mean.density, season = "dry"), type = "response")
trans.H.D.per.year.wet = predict(transition, newdata = data.frame(year = vr.per.year$year, density = mean.density, season = "rain"), type = "response")
vr.per.year$transHD.dry = trans.H.D.per.year.dry
vr.per.year$transHD.wet = trans.H.D.per.year.wet


# Recruitment of helpers

coefs.recruitment.H = coef(recruitment.H)

recruit.H.per.year.dry = exp(coefs.recruitment.H$`(Intercept)` + (coefs.recruitment.H$density * mean.density))
recruit.H.per.year.wet = exp(coefs.recruitment.H$`(Intercept)` + (coefs.recruitment.H$density * mean.density) + coefs.recruitment.H$seasonrain)
vr.per.year$recruitH.dry = recruit.H.per.year.dry
vr.per.year$recruitH.wet = recruit.H.per.year.wet


# Recruitment of dominants

coefs.recruitment.D = coef(recruitment.D)

recruit.D.per.year.dry = exp(coefs.recruitment.D$`(Intercept)` + (coefs.recruitment.D$density * mean.density) + (coefs.recruitment.D$`I(density^2)` * mean.density.squared))
recruit.D.per.year.wet = exp(coefs.recruitment.D$`(Intercept)` + (coefs.recruitment.D$density * mean.density) + coefs.recruitment.D$seasonrain + (coefs.recruitment.D$`I(density^2)` * mean.density.squared) + (coefs.recruitment.D$`seasonrain:density` * mean.density) + (coefs.recruitment.D$`seasonrain:I(density^2)` * mean.density.squared))
vr.per.year$recruitD.dry = recruit.D.per.year.dry
vr.per.year$recruitD.wet = recruit.D.per.year.wet


# Emigration of helpers

emig.per.year.dry = predict(emig, newdata = data.frame(year = vr.per.year$year, density = mean.density, season = "dry"), type = "response")
emig.per.year.wet = predict(emig, newdata = data.frame(year = vr.per.year$year, density = mean.density, season = "rain"), type = "response")
vr.per.year$emig.dry = emig.per.year.dry
vr.per.year$emig.wet = emig.per.year.wet




###########################################################################
#
# 4. Saving files and data
#
###########################################################################

write.csv(vr.per.year, "Meerkats_VitalRatesEstimates.csv", row.names = F)

save(survJ, file = "GLMM_survJ.RData")
save(survS, file = "GLMM_survS.RData")
save(survH, file = "GLMM_survH.RData")
save(survD, file = "GLMM_survD.RData")
save(recruitment.H, file = "GLMM_recruitH.RData")
save(recruitment.D, file = "GLMM_recruitD.RData")
save(transition, file = "GLMM_transHD.RData")
save(emig, file = "GLMM_emig.RData")
save(pop.range, file = "GLMM_poprange.RData")
