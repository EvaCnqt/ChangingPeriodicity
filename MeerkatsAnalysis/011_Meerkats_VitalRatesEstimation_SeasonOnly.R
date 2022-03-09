############################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script uses the capture-recapture data of a meerkat population in the Kalahari desert, collected between 1998 and 2016. 
# The aim of this script is to model the survival, transitions, emigration, and recruitment of four life-history stages: juvenile, subadult, helper, and dominant. We model as well the population range. 
#
# We focus here only on the effect of the season on vital rates. That is, 
# based on the original models, we remove the effect of density and keep
# only season for models that already had a season effect or just an
# intercept for models without a season effect.
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
  library(lme4)
  library(ggplot2)
  library(MASS)
  library(nlme)
  library(bbmle)
  library(optimx)
  library(MuMIn)
}

load.librairies()


## 1.3. Loading functions ----
# -----------------------

# Function to test for overdispersion in GLMMs.
# Taken from material from Prof. Ben Bolker (accessible at https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


# Function to compute 95 % CIs for glmmPQL models
# Taken from Prof. Marc Girondot (accessible at https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

easyPredCI <- function(model, newdata=NULL, alpha=0.05) {
  # Marc Girondot - 2016-01-09
  if (is.null(newdata)) {
    if (any(class(model)=="glmerMod")) newdata <- model@frame
    if (any(class(model)=="glmmPQL") | any(class(model)=="glm")) newdata <- model$data
    if (any(class(model)=="glmmadmb")) newdata <- model$frame
  }
  
  ## baseline prediction, on the linear predictor scale:
  pred0 <- predict(model, re.form=NA, newdata=newdata)
  ## fixed-effects model matrix for new data
  if (any(class(model)=="glmmadmb")) {
    X <- model.matrix(delete.response(model$terms), newdata)
  } else {
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
                      newdata)
  }
  
  if (any(class(model)=="glm")) {
    # Marc Girondot - 2016-01-09
    # Note that beta is not used
    beta <- model$coefficients
  } else {
    beta <- fixef(model) ## fixed-effects coefficients
  }
  
  V <- vcov(model)     ## variance-covariance matrix of beta
  
  # Marc Girondot - 2016-01-09
  if (any(!(colnames(V) %in% colnames(X)))) {
    dfi <- matrix(data = rep(0, dim(X)[1]*sum(!(colnames(V) %in% colnames(X)))), nrow = dim(X)[1])
    colnames(dfi) <- colnames(V)[!(colnames(V) %in% colnames(X))]
    X <- cbind(X, dfi)
  }
  
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  
  ## inverse-link function
  # Marc Girondot - 2016-01-09
  if (any(class(model)=="glmmPQL") | any(class(model)=="glm")) linkinv <- model$family$linkinv
  if (any(class(model)=="glmerMod")) linkinv <- model@resp$family$linkinv
  if (any(class(model)=="glmmadmb")) linkinv <- model$ilinkfun
  
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(lwr=pred0-crit*pred.se,
                upr=pred0+crit*pred.se))
}


## 1.4. Loading and prepating data ----
# --------------------------------

data.meerkats = read.csv("MeerkatsData.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Turning year into a factor

data.meerkats$year = as.factor(data.meerkats$year)

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

# Original model:
# survJ = glmer(surv ~ season + density + density2 + season:density + season:density2 + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

survJ = glmer(surv ~ season + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

# Plotting the observed data vs the predictions of the model

survJ.obs.rain = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "rain", ], FUN = mean)
survJ.obs.dry = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "dry", ], FUN = mean)

  
  # Rain

survJ.pred.rain = data.frame(year = unique(data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "rain", ]$year), pred = predict(survJ, newdata = data.frame(season = "rain", year = unique(data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "rain", ]$year)), type = "response"))


  # Obs vs. Pred

plot(survJ.obs.rain$surv, survJ.pred.rain$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of juveniles (J) survival - Rain season")
abline(a = 0, b = 1, col = "red")


  # Dry

survJ.pred.dry = data.frame(year = unique(data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "dry", ]$year), pred = predict(survJ, newdata = data.frame(season = "dry", year = unique(data.meerkats[data.meerkats$stage == "J" & data.meerkats$season == "dry", ]$year)), type = "response"))


  # Obs vs. Pred

plot(survJ.obs.dry$surv, survJ.pred.dry$pred, 
     pch = 16, xlim = c(0.6, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of juveniles (J) survival - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survJ.new.data = expand.grid(year = unique(data.meerkats$year),
                             season = unique(data.meerkats$season))


# 95% confidence intervals

bootstrap.survJ <- bootMer(survJ, FUN = function(x) predict(x, survJ.new.data, re.form = NA),
                           nsim = 1000)
survJ.new.data$lwr = inv.logit(apply(bootstrap.survJ$t, 2, quantile, 0.025, na.rm = T))
survJ.new.data$upr = inv.logit(apply(bootstrap.survJ $t, 2, quantile, 0.975, na.rm = T))
survJ.new.data$pred = predict(survJ, newdata = survJ.new.data, type = "response", re.form = NA)


# Plotting the predictions 

plot.survJ = ggplot(survJ.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_x_discrete(labels = c("Rain", "Dry")) +
  xlab("Season") +
  ylab("Survival probability") +
  ggtitle("Juvenile survival") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# survS = glmer(surv ~ density + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

survS = glmer(surv ~ 1 + (1 + season|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

# Plotting the observed data vs the predictions of the model

survS.obs.rain = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "rain", ], FUN = mean)
survS.obs.dry = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "dry", ], FUN = mean)


  # Rain

survS.pred.rain = data.frame(year = unique(data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "rain", ]$year), 
                             pred = predict(survS, newdata = data.frame(year = unique(data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "rain", ]$year), season = "rain"), type = "response"))


# Obs vs. Pred

plot(survS.obs.rain$surv, survS.pred.rain$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of subadults (S) survival - Rain season (w/o season)")
abline(a = 0, b = 1, col = "red")


  # Dry

survS.pred.dry = data.frame(year = unique(data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "dry", ]$year), 
                            pred = predict(survS, newdata = data.frame(year = unique(data.meerkats[data.meerkats$stage == "S" & data.meerkats$season == "dry", ]$year), season = "dry"), type = "response"))


# Obs vs. Pred

plot(survS.obs.dry$surv, survS.pred.dry$pred, 
     pch = 16, xlim = c(0.6, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of subadults (S) survival - Dry season (w/o season)")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survS.new.data = expand.grid(year = unique(data.meerkats$year),
                             season = unique(data.meerkats$season))


# 95% confidence intervals
bootstrap.survS <- bootMer(survS, FUN = function(x) predict(x, survS.new.data, re.form = NA),
                           nsim = 1000)
survS.new.data$lwr = inv.logit(apply(bootstrap.survS$t, 2, quantile, 0.025, na.rm = T))
survS.new.data$upr = inv.logit(apply(bootstrap.survS$t, 2, quantile, 0.975, na.rm = T))
survS.new.data$pred = predict(survS, newdata = survS.new.data, type = "response", re.form = NA)


# Plotting the predictions
 
plot.survS = ggplot(survS.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  xlab("Population density (indiv./km2)") +
  ylab("Survival probability") +
  ggtitle("Juvenile survival") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# survH = glmer(surv ~ season + (1 + season|year), data = helpers, family = binomial)

survH = glmer(surv ~ season + (1 + season|year), data = helpers, family = binomial)


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
                           nsim = 1000)
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
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# survD = glmer(surv ~ 1 + (1 + season|year), data = dominants, family = binomial) 

survD = glmer(surv ~ 1 + (1 + season|year), data = dominants, family = binomial)

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
                             season = unique(data.meerkats$season))


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
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# emig = glmer(emig ~ density + (1 + season|year), data = helpers, family = binomial)

emig = glmer(emig ~ 1 + (1 + season|year), data = helpers, family = binomial)


# Plotting the observed data vs the predictions of the model

emig.obs.rain = aggregate(emig ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
emig.obs.dry = aggregate(emig ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


# Rain

emig.pred.rain = expand.grid(season = "rain", 
                             year = levels(data.meerkats$year))
emig.pred.rain$pred = predict(emig, newdata = emig.pred.rain, type = "response")
emig.pred.rain = aggregate(pred ~ year, data = emig.pred.rain, FUN = mean)


# Obs vs. Pred

plot(emig.obs.rain$emig, emig.pred.rain$pred,
     pch = 16, ylim = c(0, 1), xlab = "Observed emigration probability", ylab = "Predicted emigration probability",
     main = "Observed vs. Predicted values of helpers (H) emigration - Rain season")
abline(a = 0, b = 1, col = "red")


# Dry

emig.pred.dry = expand.grid(season = "dry",
                            year = levels(data.meerkats$year))
emig.pred.dry$pred = predict(emig, newdata = emig.pred.dry, type = "response")
emig.pred.dry = aggregate(pred ~ year, data = emig.pred.dry, FUN = mean)


# Obs vs. Pred

plot(emig.obs.dry$emig, emig.pred.dry$pred, 
     pch = 16, ylim = c(0, 1), xlab = "Observed emigration probability", ylab = "Predicted emigration probability",
     main="Observed vs. Predicted values of helpers (H) emigration - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

emig.new.data = expand.grid(year = levels(data.meerkats$year),
                            season = unique(data.meerkats$season))


# 95% confidence intervals

bootstrap.emig <- bootMer(emig, FUN = function(x) predict(x, emig.new.data, re.form = NA),
                          nsim = 1000)
emig.new.data$lwr = inv.logit(apply(bootstrap.emig$t, 2, quantile, 0.025, na.rm = T))
emig.new.data$upr = inv.logit(apply(bootstrap.emig$t, 2, quantile, 0.975, na.rm = T))
emig.new.data$pred = predict(emig, newdata = emig.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.emig = ggplot(emig.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  scale_x_discrete(labels = c("Rain", "Dry")) +
  xlab("Season") +
  ylab("Emigration probability") +
  ggtitle("Helper emigration") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# transition = glmer(transition ~ density  + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

transition = glmer(transition ~ 1  + (1 + season|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Plotting the observed data vs the predictions of the model

transition.obs.rain = aggregate(transition ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
transition.obs.dry = aggregate(transition ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


# Rain

transition.pred.rain = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ]$year),
                                  pred = predict(transition, newdata = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ]$year), season = "rain"), type = "response"))


# Obs vs. Pred

plot(transition.obs.rain$transition, transition.pred.rain$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed transition probability", ylab = "Predicted transition probability",
     main="Observed vs. Predicted values of H-D transition - Rain season")
abline(a = 0, b = 1, col = "red")


# Dry

transition.pred.dry = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ]$year),
                                  pred = predict(transition, newdata = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ]$year), season = "dry"), type = "response"))


# Obs vs. Pred

plot(transition.obs.dry$transition, transition.pred.dry$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed transition probability", ylab = "Predicted transition probability",
     main="Observed vs. Predicted values of H-D transition - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transition.new.data = expand.grid(year = levels(data.meerkats$year),
                                  season = unique(data.meerkats$season))


# 95% confidence intervals

bootstrap.transition <- bootMer(transition, FUN = function(x) predict(x, transition.new.data, re.form = NA),
                                nsim = 1000)
transition.new.data$lwr = inv.logit(apply(bootstrap.transition$t, 2, quantile, 0.025, na.rm = T))
transition.new.data$upr = inv.logit(apply(bootstrap.transition$t, 2, quantile, 0.975, na.rm = T))
transition.new.data$pred = predict(transition, newdata = transition.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.transition = ggplot(transition.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_x_discrete(labels = c("Rain", "Dry")) +
  xlab("Season") +
  ylab("Transition probability") +
  ggtitle("Helper-Dominant transition") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# recruitment.H = glmer(pups ~ density + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

recruitment.H = glmer(pups ~ 1 + (1 + season|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(recruitment.H)

recruitmentH.quasiP = MASS::glmmPQL(pups ~ 1, random =  ~ 1+season|year, family = poisson, data = helpers)

summary(recruitmentH.quasiP)

recruitment.H = recruitmentH.quasiP


# Plotting the observed data vs the predictions of the model

recruitment.H.obs.rain = aggregate(pups ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ], FUN = mean)
recruitment.H.obs.dry = aggregate(pups ~ year, data = data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ], FUN = mean)


# Rain  

recruitment.H.pred.rain = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ]$year),
                                     pred = predict(recruitment.H, newdata = expand.grid(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "rain", ]$year), season = unique(data.meerkats$season)), type = "response", level = 0)[1:19])


# Obs vs. Pred

plot(recruitment.H.obs.rain$pups, recruitment.H.pred.rain$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of helper (H) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")


# Dry

recruitment.H.pred.dry = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ]$year),
                                    pred = predict(recruitment.H, newdata = data.frame(year = unique(data.meerkats[data.meerkats$stage == "H" & data.meerkats$season == "dry", ]$year), season = "dry"), type = "response", level = 0))


# Obs vs. Pred

plot(recruitment.H.obs.dry$pups, recruitment.H.pred.dry$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of helper (H) recruitment - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
recruitment.H.new.data = expand.grid(year = unique(data.meerkats$year),
                                      season = unique(data.meerkats$season))
recruitment.H.new.data$lwr = NA
recruitment.H.new.data$upr = NA
recruitment.H.new.data$pred = NA

# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function from Prof. Marc Girondot (https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

easyPredCI <- function(model, newdata=NULL, alpha=0.05) {
  # Marc Girondot - 2016-01-09
  if (is.null(newdata)) {
    if (any(class(model)=="glmerMod")) newdata <- model@frame
    if (any(class(model)=="glmmPQL") | any(class(model)=="glm")) newdata <- model$data
    if (any(class(model)=="glmmadmb")) newdata <- model$frame
  }
  
  ## baseline prediction, on the linear predictor scale:
  pred0 <- predict(model, re.form=NA, newdata=newdata)
  ## fixed-effects model matrix for new data
  if (any(class(model)=="glmmadmb")) {
    X <- model.matrix(delete.response(model$terms), newdata)
  } else {
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
                      newdata)
  }
  
  if (any(class(model)=="glm")) {
    # Marc Girondot - 2016-01-09
    # Note that beta is not used
    beta <- model$coefficients
  } else {
    beta <- fixef(model) ## fixed-effects coefficients
  }
  
  V <- vcov(model)     ## variance-covariance matrix of beta
  
  # Marc Girondot - 2016-01-09
  if (any(!(colnames(V) %in% colnames(X)))) {
    dfi <- matrix(data = rep(0, dim(X)[1]*sum(!(colnames(V) %in% colnames(X)))), nrow = dim(X)[1])
    colnames(dfi) <- colnames(V)[!(colnames(V) %in% colnames(X))]
    X <- cbind(X, dfi)
  }
  
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  
  ## inverse-link function
  # Marc Girondot - 2016-01-09
  if (any(class(model)=="glmmPQL") | any(class(model)=="glm")) linkinv <- model$family$linkinv
  if (any(class(model)=="glmerMod")) linkinv <- model@resp$family$linkinv
  if (any(class(model)=="glmmadmb")) linkinv <- model$ilinkfun
  
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(lwr=pred0-crit*pred.se,
                upr=pred0+crit*pred.se))
}

recruitment.H.new.data$lwr[which(recruitment.H.new.data$season == "rain")] = apply(recruitment.H.new.data[which(recruitment.H.new.data$season == "rain"), ], 1, FUN = function(x) easyPredCI(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[1, 1])
recruitment.H.new.data$lwr[which(recruitment.H.new.data$season == "dry")] = apply(recruitment.H.new.data[which(recruitment.H.new.data$season == "dry"), ], 1, FUN = function(x) easyPredCI(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[2, 1])

recruitment.H.new.data$upr[which(recruitment.H.new.data$season == "rain")] = apply(recruitment.H.new.data[which(recruitment.H.new.data$season == "rain"), ], 1, FUN = function(x) easyPredCI(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[1, 2])
recruitment.H.new.data$upr[which(recruitment.H.new.data$season == "dry")] = apply(recruitment.H.new.data[which(recruitment.H.new.data$season == "dry"), ], 1, FUN = function(x) easyPredCI(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[2, 2])

recruitment.H.new.data$pred[which(recruitment.H.new.data$season == "rain")] = apply(recruitment.H.new.data[which(recruitment.H.new.data$season == "rain"), ], 1, FUN = function(x) predict(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")), type = "response")[1])
recruitment.H.new.data$pred[which(recruitment.H.new.data$season == "dry")] = apply(recruitment.H.new.data[which(recruitment.H.new.data$season == "dry"), ], 1, FUN = function(x) predict(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")), type = "response")[2])

recruitment.H.RE.avg = data.frame(season = aggregate(pred ~ season, data = recruitment.H.new.data, mean)$season,
                                   pred = aggregate(pred ~ season, data = recruitment.H.new.data, mean)$pred,
                                   lwr = aggregate(lwr ~ season, data = recruitment.H.new.data, mean)$lwr,
                                   upr = aggregate(upr ~ season, data = recruitment.H.new.data, mean)$upr)

plot.recruitment.H = ggplot(recruitment.H.RE.avg, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_x_discrete(labels = c("Rain", "Dry")) +
  xlab("Season") +
  ylab("Number of pups") +
  ggtitle("Helper recruitment") +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[1], cbbPalette[2])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Wet", "Dry"),
                      values = c(cbbPalette[1], cbbPalette[2])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l =0 )), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), 
        legend.text = element_text(size = 40), 
        legend.title = element_text(size = 40), 
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

# Original model:
# recruitment.D = glmer(pups ~ season + density + I(density^2) + season:density + season:I(density^2) + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

recruitment.D = glmer(pups ~ season + (1 + season|year), data = dominants, family = poisson, control = glmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(recruitment.D)

recruitmentD.quasiP = MASS::glmmPQL(pups ~ season, random =  ~ 1+season|year, family = poisson, data = dominants)

summary(recruitmentD.quasiP)

recruitment.D = recruitmentD.quasiP

# Plotting the observed data vs the predictions of the model

recruitment.D.obs.rain = aggregate(pups ~ year, data = data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "rain", ], FUN = mean)
recruitment.D.obs.dry = aggregate(pups ~ year, data = data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "dry", ], FUN = mean)


  # Rain

recruitment.D.pred.rain = data.frame(year = unique(data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "rain", ]$year),
                                     pred = predict(recruitment.D, newdata = expand.grid(year = unique(data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "rain", ]$year), season = unique(data.meerkats$season)), type = "response", level = 0)[1:19])

# Obs vs. Pred

plot(recruitment.D.obs.rain$pups, recruitment.D.pred.rain$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")


  # Dry

recruitment.D.pred.dry = data.frame(year = unique(data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "dry", ]$year),
                                    pred = predict(recruitment.D, newdata = expand.grid(year = unique(data.meerkats[data.meerkats$stage == "D" & data.meerkats$season == "dry", ]$year), season = unique(data.meerkats$season)), type = "response", level = 0)[20:38])

# Obs vs. Pred

plot(recruitment.D.obs.dry$pups, recruitment.D.pred.dry$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Dry season")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
recruitment.D.new.data = expand.grid(year = unique(data.meerkats$year),
                                     season = unique(data.meerkats$season))
recruitment.D.new.data$lwr = NA
recruitment.D.new.data$upr = NA
recruitment.D.new.data$pred = NA

# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function from Prof. Marc Girondot (https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

recruitment.D.new.data$lwr[which(recruitment.D.new.data$season == "rain")] = apply(recruitment.D.new.data[which(recruitment.D.new.data$season == "rain"), ], 1, FUN = function(x) easyPredCI(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[1, 1])
recruitment.D.new.data$lwr[which(recruitment.D.new.data$season == "dry")] = apply(recruitment.D.new.data[which(recruitment.D.new.data$season == "dry"), ], 1, FUN = function(x) easyPredCI(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[2, 1])

recruitment.D.new.data$upr[which(recruitment.D.new.data$season == "rain")] = apply(recruitment.D.new.data[which(recruitment.D.new.data$season == "rain"), ], 1, FUN = function(x) easyPredCI(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[1, 2])
recruitment.D.new.data$upr[which(recruitment.D.new.data$season == "dry")] = apply(recruitment.D.new.data[which(recruitment.D.new.data$season == "dry"), ], 1, FUN = function(x) easyPredCI(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")))[2, 2])

recruitment.D.new.data$pred[which(recruitment.D.new.data$season == "rain")] = apply(recruitment.D.new.data[which(recruitment.D.new.data$season == "rain"), ], 1, FUN = function(x) predict(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")), type = "response")[1])
recruitment.D.new.data$pred[which(recruitment.D.new.data$season == "dry")] = apply(recruitment.D.new.data[which(recruitment.D.new.data$season == "dry"), ], 1, FUN = function(x) predict(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), season = c("rain", "dry")), type = "response")[2])

recruitment.D.RE.avg = data.frame(season = aggregate(pred ~ season, data = recruitment.D.new.data, mean)$season,
                                    pred = aggregate(pred ~ season, data = recruitment.D.new.data, mean)$pred,
                                    lwr = aggregate(lwr ~ season, data = recruitment.D.new.data, mean)$lwr,
                                    upr = aggregate(upr ~ season, data = recruitment.D.new.data, mean)$upr)

plot.recruitment.D.FigS4 = ggplot(recruitment.D.RE.avg, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_x_discrete(labels = c("Rain", "Dry")) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  ylim(0, 5) +
  scale_fill_manual(name = "Season",
                    breaks = c("rain", "dry"),
                    labels = c("Wet", "Dry"),
                    values = c(cbbPalette[1], cbbPalette[2])) +
  scale_colour_manual(name = "Season",
                      breaks = c("rain", "dry"),
                      labels = c("Wet", "Dry"),
                      values = c(cbbPalette[1], cbbPalette[2])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9), 
        legend.position = "none", 
        legend.key.height = unit(2, "lines"), 
        legend.key.width = unit(2, "lines"))




###########################################################################
#
# 3. Saving the vital rates ----
#
###########################################################################

vr.per.year = data.frame(year = unique(data.meerkats$year))


# Survival of juveniles

surv.J.per.year.dry = predict(survJ, newdata = data.frame(year = vr.per.year$year, season = "dry"), type = "response")
surv.J.per.year.wet = predict(survJ, newdata = data.frame(year = vr.per.year$year, season = "rain"), type = "response")
vr.per.year$survJ.dry = surv.J.per.year.dry
vr.per.year$survJ.wet = surv.J.per.year.wet


# Survival of subadults

surv.S.per.year.dry = predict(survS, newdata = data.frame(year = vr.per.year$year, season = "dry"), type = "response")
surv.S.per.year.wet = predict(survS, newdata = data.frame(year = vr.per.year$year, season = "rain"), type = "response")
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

trans.H.D.per.year.dry = predict(transition, newdata = data.frame(year = vr.per.year$year, season = "dry"), type = "response")
trans.H.D.per.year.wet = predict(transition, newdata = data.frame(year = vr.per.year$year, season = "rain"), type = "response")
vr.per.year$transHD.dry = trans.H.D.per.year.dry
vr.per.year$transHD.wet = trans.H.D.per.year.wet


# Recruitment of helpers

coefs.recruitment.H = coef(recruitment.H)

recruit.H.per.year.dry = exp(coefs.recruitment.H$`(Intercept)`)
recruit.H.per.year.wet = exp(coefs.recruitment.H$`(Intercept)` + coefs.recruitment.H$seasonrain)
vr.per.year$recruitH.dry = recruit.H.per.year.dry
vr.per.year$recruitH.wet = recruit.H.per.year.wet


# Recruitment of dominants

coefs.recruitment.D = coef(recruitment.D)

recruit.D.per.year.dry = exp(coefs.recruitment.D$`(Intercept)`)
recruit.D.per.year.wet = exp(coefs.recruitment.D$`(Intercept)` + coefs.recruitment.D$seasonrain)
vr.per.year$recruitD.dry = recruit.D.per.year.dry
vr.per.year$recruitD.wet = recruit.D.per.year.wet


# Emigration of helpers

emig.per.year.dry = predict(emig, newdata = data.frame(year = vr.per.year$year, season = "dry"), type = "response")
emig.per.year.wet = predict(emig, newdata = data.frame(year = vr.per.year$year, season = "rain"), type = "response")
vr.per.year$emig.dry = emig.per.year.dry
vr.per.year$emig.wet = emig.per.year.wet




###########################################################################
#
# 4. Saving files and data
#
###########################################################################

write.csv(vr.per.year, "Meerkats_VitalRatesEstimates_SeasonOnly.csv", row.names = F)

save(survJ, file = "GLMM_survJ_SeasonOnly.RData")
save(survS, file = "GLMM_survS_SeasonOnly.RData")
save(survH, file = "GLMM_survH_SeasonOnly.RData")
save(survD, file = "GLMM_survD_SeasonOnly.RData")
save(recruitment.H, file = "GLMM_recruitH_SeasonOnly.RData")
save(recruitment.D, file = "GLMM_recruitD_SeasonOnly.RData")
save(transition, file = "GLMM_transHD_SeasonOnly.RData")
save(emig, file = "GLMM_emig_SeasonOnly.RData")
