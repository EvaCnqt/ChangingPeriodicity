##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
#
# This script uses the capture-recapture data of a meerkat population in the Kalahari desert, 
# collected between 1998 and 2016. 
# The aim of this script is to model the survival, transitions, emigration, and recruitment of four 
# life-history stages: juvenile, subadult, helper, and dominant. We model as well the population range. 
#
# We focus here only on the effect of density on vital rates. That is, 
# based on the original models, we remove the effect of season and keep
# only density for models that already had a density effect or just an
# intercept for models without a density effect.
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


## 1.4. Loading and preparing data ----
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
# 2. Fitting Generalized Linear Mixed Models (GLMMs) ----
# to estimate the vital rates
#
###########################################################################

## 2.1. Survival (surv) ----
# ---------------------

## 2.1.1. Juveniles (J) ----
# ---------------------

# Original model:
# survJ = glmer(surv ~ season + density + I(density^2) + season:density + season:I(density^2) + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

survJ = glmer(surv ~ density + I(density^2) + (1|year), data = juveniles, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Plotting the observed data vs the predictions of the model

survJ.obs = aggregate(surv ~ year + density, data = juveniles, FUN = mean)

survJ.pred = data.frame(density = seq(min(survJ.obs$density), max(survJ.obs$density), length.out = length(survJ.obs$density)), 
                        pred = predict(survJ, newdata = data.frame(density = seq(min(survJ.obs$density), max(survJ.obs$density), length.out = length(survJ.obs$density))), type = "response", re.form = NA))


  # Obs vs. Pred

plot(survJ.obs$surv, survJ.pred$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of juveniles (J) survival")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survJ.new.data = expand.grid(year = unique(data.meerkats$year),
                             density = seq(min(data.meerkats$density), max(data.meerkats$density), length.out = 100))


# 95% confidence intervals

bootstrap.survJ <- bootMer(survJ, FUN = function(x) predict(x, survJ.new.data, re.form = NA),
                           nsim = 1000)
survJ.new.data$lwr = inv.logit(apply(bootstrap.survJ$t, 2, quantile, 0.025, na.rm = T))
survJ.new.data$upr = inv.logit(apply(bootstrap.survJ $t, 2, quantile, 0.975, na.rm = T))
survJ.new.data$pred = predict(survJ, newdata = survJ.new.data, type = "response", re.form = NA)


# Plotting the predictions 

plot.survJ = ggplot(survJ.new.data, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
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

survS = glmer(surv ~ density + (1|year), data = subadults, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

# Plotting the observed data vs the predictions of the model

survS.obs = aggregate(surv ~ year + density, data = data.meerkats[data.meerkats$stage == "S", ], FUN = mean)

survS.pred = data.frame(density = seq(min(survS.obs$density), max(survS.obs$density), length.out = length(survS.obs$density)), 
                             pred = predict(survS, newdata = data.frame(density = seq(min(survS.obs$density), max(survS.obs$density), length.out = length(survS.obs$density))), type = "response", re.form = NA))


# Obs vs. Pred

plot(survS.obs$surv, survS.pred$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.6, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of subadults (S) survival")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survS.new.data = expand.grid(year = unique(data.meerkats$year),
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
  ggtitle("Subadult survival") +
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

survH = glmer(surv ~ 1 + (1|year), data = helpers, family = binomial)


# Plotting the observed data vs the predictions of the model

survH.obs = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "H", ], FUN = mean)

survH.pred = data.frame(year = levels(data.meerkats$year), 
                             pred = predict(survH, newdata = data.frame(year = levels(data.meerkats$year)), type = "response"))


# Obs vs. Pred

plot(survH.obs$surv, survH.pred$pred, 
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of helpers (H) survival")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survH.new.data = expand.grid(year = levels(data.meerkats$year))


# 95% confidence intervals

bootstrap.survH <- bootMer(survH, FUN = function(x) predict(x, survH.new.data),
                           nsim = 1000)
survH.new.data$lwr = inv.logit(apply(bootstrap.survH$t, 2, quantile, 0.025, na.rm = T))
survH.new.data$upr = inv.logit(apply(bootstrap.survH$t, 2, quantile, 0.975, na.rm = T))
survH.new.data$pred = predict(survH, newdata = survH.new.data, type = "response")


# Plotting the predictions

plot.survH = ggplot(survH.new.data, aes(x = as.numeric(year), y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  scale_x_discrete(labels = as.numeric(survH.new.data$year)) +
  xlab("Year") +
  ylab("Helper survival probability") +
  #ggtitle("Helper survival") +
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

survD = glmer(surv ~ 1 + (1|year), data = dominants, family = binomial)


# Plotting the observed data vs the predictions of the model

survD.obs = aggregate(surv ~ year, data = data.meerkats[data.meerkats$stage == "D", ], FUN = mean)

survD.pred = data.frame(year = levels(data.meerkats$year), 
                        pred = predict(survD, newdata = data.frame(year = levels(data.meerkats$year)), type = "response"))


# Obs vs. Pred

plot(survD.obs$surv, survD.pred$pred,
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability",
     main = "Observed vs. Predicted values of dominants (D) survival")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survD.new.data = expand.grid(year = levels(data.meerkats$year))


# 95% confidence intervals

bootstrap.survD <- bootMer(survD, FUN = function(x) predict(x, survD.new.data),
                           nsim = 1000)
survD.new.data$lwr = inv.logit(apply(bootstrap.survD$t, 2, quantile, 0.025, na.rm = T))
survD.new.data$upr = inv.logit(apply(bootstrap.survD$t, 2, quantile, 0.975, na.rm = T))
survD.new.data$pred = predict(survD, newdata = survD.new.data, type = "response")


# Plotting the predictions

plot.survD = ggplot(survD.new.data, aes(x = as.numeric(year), y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Survival probability") +
  ggtitle("Dominant survival") +
  scale_x_discrete(labels = as.numeric(survD.new.data$year)) +
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

emig = glmer(emig ~ density + (1|year), data = helpers, family = binomial)


# Plotting the observed data vs the predictions of the model

emig.obs = aggregate(emig ~ year + density, data = data.meerkats[data.meerkats$stage == "H", ], FUN = mean)

emig.pred = data.frame(year = levels(emig.obs$year), 
                       density = unique(emig.obs$density))
emig.pred$pred = predict(emig, newdata = emig.pred, type = "response")

# Obs vs. Pred

plot(emig.obs$emig, emig.pred$pred,
     pch = 16, ylim = c(0, 1), xlab = "Observed emigration probability", ylab = "Predicted emigration probability",
     main = "Observed vs. Predicted values of helpers (H) emigration")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

emig.new.data = expand.grid(year = levels(data.meerkats$year),
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

transition = glmer(transition ~ density  + (1|year), data = helpers, family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Plotting the observed data vs the predictions of the model

transition.obs = aggregate(transition ~ year + density, data = data.meerkats[data.meerkats$stage == "H", ], FUN = mean)

transition.pred = data.frame(density = seq(min(transition.obs$density), max(transition.obs$density), length.out = length(transition.obs$density)),
                             pred = predict(transition, newdata = data.frame(density = seq(min(transition.obs$density), max(transition.obs$density), length.out = length(transition.obs$density))), type = "response", re.form = NA))


# Obs vs. Pred

plot(transition.obs$transition, transition.pred$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed transition probability", ylab = "Predicted transition probability",
     main="Observed vs. Predicted values of H-D transition")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transition.new.data = expand.grid(year = levels(data.meerkats$year),
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

recruitment.H = glmer(pups ~ density + (1|year), data = helpers, family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(recruitment.H)

recruitmentH.quasiP = MASS::glmmPQL(pups ~ density, random =  ~ 1|year, family = poisson, data = helpers)

summary(recruitmentH.quasiP)

recruitment.H = recruitmentH.quasiP


# Plotting the observed data vs the predictions of the model

recruitment.H.obs = aggregate(pups ~ year + density, data = data.meerkats[data.meerkats$stage == "H", ], FUN = mean)

recruitment.H.pred = data.frame(density = seq(min(recruitment.H.obs$density), max(recruitment.H.obs$density), length.out = length(recruitment.H.obs$density)),
                                pred = predict(recruitment.H, newdata = expand.grid(density = seq(min(recruitment.H.obs$density), max(recruitment.H.obs$density), length.out = length(recruitment.H.obs$density))), type = "response", re.form = NA, level = 0))


# Obs vs. Pred

plot(recruitment.H.obs$pups, recruitment.H.pred$pred,
     pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2), xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of helper (H) recruitment")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
recruitment.H.new.data = expand.grid(year = unique(data.meerkats$year),
                                     density = seq(min(data.meerkats$density),15,length.out = 100))
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

recruitment.H.new.data$lwr = apply(recruitment.H.new.data, 1, FUN = function(x) easyPredCI(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), density = as.numeric(x[2])))[, 1])

recruitment.H.new.data$upr = apply(recruitment.H.new.data, 1, FUN = function(x) easyPredCI(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), density = as.numeric(x[2])))[, 2])

recruitment.H.new.data$pred = apply(recruitment.H.new.data, 1, FUN = function(x) predict(recruitment.H, newdata = data.frame(year = as.numeric(x[1]), density = as.numeric(x[2])), type = "response")[1])

recruitment.H.RE.avg = data.frame(density = aggregate(pred ~ density, data = recruitment.H.new.data, mean)$density,
                                   pred = aggregate(pred ~ density, data = recruitment.H.new.data, mean)$pred,
                                   lwr = aggregate(lwr ~ density, data = recruitment.H.new.data, mean)$lwr,
                                   upr = aggregate(upr ~ density, data = recruitment.H.new.data, mean)$upr)

plot.recruitment.H = ggplot(recruitment.H.RE.avg, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Population density (indiv./km2)") +
  ylab("Number of pups") +
  ggtitle("Helper recruitment") +
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

recruitment.D = glmer(pups ~ density + I(density^2) + (1|year), data = dominants, family = poisson, control = glmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(recruitment.D)

recruitmentD.quasiP = MASS::glmmPQL(pups ~ density + I(density^2), random =  ~ 1|year, family = poisson, data = dominants)

summary(recruitmentD.quasiP)

recruitment.D = recruitmentD.quasiP

# Plotting the observed data vs the predictions of the model

recruitment.D.obs = aggregate(pups ~ year + density, data = data.meerkats[data.meerkats$stage == "D", ], FUN = mean)

recruitment.D.pred = data.frame(density = seq(min(recruitment.D.obs$density), 20, length.out = length(recruitment.D.obs$density)),
                                pred = predict(recruitment.D, newdata = expand.grid(density = seq(min(recruitment.D.obs$density), 20, length.out = length(recruitment.D.obs$density))), type = "response", re.form = NA, level = 0))

# Obs vs. Pred

plot(recruitment.D.obs$pups, recruitment.D.pred$pred,
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment",
     main = "Observed vs. Predicted values of dominant (D) recruitment - Rain season")
abline(a = 0, b = 1, col = "red")


# Most biologically relevant model

# Computing the predictions of the model and the 95 % CI
recruitment.D.new.data = expand.grid(year = unique(data.meerkats$year),
                                     density = seq(min(data.meerkats$density),15,length.out = 100))
recruitment.D.new.data$lwr = NA
recruitment.D.new.data$upr = NA
recruitment.D.new.data$pred = NA

# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function from Prof. Marc Girondot (https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

recruitment.D.new.data$lwr = apply(recruitment.D.new.data, 1, FUN = function(x) easyPredCI(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), density = as.numeric(x[2])))[, 1])

recruitment.D.new.data$upr = apply(recruitment.D.new.data, 1, FUN = function(x) easyPredCI(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), density = as.numeric(x[2])))[, 2])

recruitment.D.new.data$pred = apply(recruitment.D.new.data, 1, FUN = function(x) predict(recruitment.D, newdata = data.frame(year = as.numeric(x[1]), density = as.numeric(x[2])), type = "response")[1])

recruitment.D.RE.avg = data.frame(density = aggregate(pred ~ density, data = recruitment.D.new.data, mean)$density,
                                    pred = aggregate(pred ~ density, data = recruitment.D.new.data, mean)$pred,
                                    lwr = aggregate(lwr ~ density, data = recruitment.D.new.data, mean)$lwr,
                                    upr = aggregate(upr ~ density, data = recruitment.D.new.data, mean)$upr)

plot.recruitment.D.FigS4 = ggplot(recruitment.D.RE.avg, aes(x = density, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  ylim(0, 20) +
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

surv.J.per.year = predict(survJ, newdata = data.frame(year = vr.per.year$year, density = mean.density, density2 = mean.density.squared), type = "response")
vr.per.year$survJ = surv.J.per.year


# Survival of subadults

surv.S.per.year = predict(survS, newdata = data.frame(year = vr.per.year$year, density = mean.density), type = "response")
vr.per.year$survS = surv.S.per.year


# Survival of helpers

surv.H.per.year = predict(survH, newdata = data.frame(year = vr.per.year$year), type = "response")
vr.per.year$survH = surv.H.per.year


# Survival of dominants

surv.D.per.year = predict(survD, newdata = data.frame(year = vr.per.year$year), type = "response")
vr.per.year$survD = surv.D.per.year


# Transition H -> D

trans.H.D.per.year = predict(transition, newdata = data.frame(year = vr.per.year$year, density = mean.density), type = "response")
vr.per.year$transHD = trans.H.D.per.year


# Recruitment of helpers

coefs.recruitment.H = coef(recruitment.H)

recruit.H.per.year = exp(coefs.recruitment.H$`(Intercept)` + (coefs.recruitment.H$density * mean.density))
vr.per.year$recruitH = recruit.H.per.year


# Recruitment of dominants

coefs.recruitment.D = coef(recruitment.D)

recruit.D.per.year = exp(coefs.recruitment.D$`(Intercept)` + (coefs.recruitment.D$density * mean.density) + (coefs.recruitment.D$`I(density^2)` * mean.density.squared))
vr.per.year$recruitD = recruit.D.per.year


# Emigration of helpers

emig.per.year = predict(emig, newdata = data.frame(year = vr.per.year$year, density = mean.density), type = "response")
vr.per.year$emig = emig.per.year




###########################################################################
#
# 4. Saving files and data ----
#
###########################################################################

write.csv(vr.per.year, "Meerkats_VitalRatesEstimates_DensityOnly.csv", row.names = F)

save(survJ, file = "GLMM_survJ_DensityOnly.RData")
save(survS, file = "GLMM_survS_DensityOnly.RData")
save(survH, file = "GLMM_survH_DensityOnly.RData")
save(survD, file = "GLMM_survD_DensityOnly.RData")
save(recruitment.H, file = "GLMM_recruitH_DensityOnly.RData")
save(recruitment.D, file = "GLMM_recruitD_DensityOnly.RData")
save(transition, file = "GLMM_transHD_DensityOnly.RData")
save(emig, file = "GLMM_emig_DensityOnly.RData")
save(pop.range, file = "GLMM_poprange_DensityOnly.RData")
