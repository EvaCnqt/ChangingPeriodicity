##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
#
# This script uses the data of a dewy pine population in South-Eastern Spain collected 
# between 2012 and 2019. 
# The aim of this script is to model the survival, transitions, and reproductive rates of four 
# life-history stages: seedlings, juveniles, and small and large reproductive adults. 
# We model the rates in the last stochastic post-fire habitat, TSF>3, in natural populations 
# (i.e., low grazing) and populations exposed to human disturbances (i.e. high grazing).
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
 library(lme4)
 library(bbmle)
 library(MASS)
 library(boot)
 library(ggplot2)
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

data.grazing = read.csv("DewyPinesTSF3plus_Data.csv")
head(data.grazing)


# Turning the year into a factor to use it as a random effect
data.grazing$time = as.factor(data.grazing$time)


# Compute average abundance per square (i.e. population density)

data.low.grazing = data.grazing[which(data.grazing$LS == "LG"), ]
data.high.grazing = data.grazing[which(data.grazing$LS == "HG"), ]
data.low.grazing$IDsquare = NA
data.high.grazing$IDsquare = NA
data.low.grazing$IDsquare = paste(data.low.grazing$transect, data.low.grazing$subQuadrat, data.low.grazing$site, sep = "_")
data.high.grazing$IDsquare = paste(data.high.grazing$transect, data.high.grazing$subQuadrat, data.high.grazing$site, sep = "_")
nbsquares.LG = aggregate(IDsquare ~ time, data = data.low.grazing, FUN = function(x) length(unique(x)))
nbsquares.HG = aggregate(IDsquare ~ time, data = data.high.grazing, FUN = function(x) length(unique(x)))

density.per.square = aggregate(density ~ IDsquare + time, data = data.low.grazing, function(x) unique(x))
yearly.density.per.square.LG = aggregate(density ~ time, data = density.per.square, function(x) sum(x))
yearly.density.per.square.LG$density = yearly.density.per.square.LG$density/nbsquares.LG$IDsquare
yearly.density.per.square.LG = yearly.density.per.square.LG[- nrow(yearly.density.per.square.LG), ]

density.per.square = aggregate(density ~ IDsquare + time, data = data.high.grazing, function(x) unique(x))
yearly.density.per.square.HG = aggregate(density ~ time, data = density.per.square, function(x) sum(x))
yearly.density.per.square.HG$density = yearly.density.per.square.HG$density/nbsquares.HG$IDsquare
yearly.density.per.square.HG = yearly.density.per.square.HG[- nrow(yearly.density.per.square.HG), ]

nbsquares.TSF3plus = data.frame()




###########################################################################
#
# 2. Fitting Generalized Linear Mixed Models (GLMMs) ----
# to estimate the vital rates for stochastic TSF>3 
#
###########################################################################

## 2.1. Survival (surv) ----
# ---------------------

## 2.1.1. Little grazed (LG) - Seedlings (SD) ----
# -------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "LG" & data.grazing$stage == "SD")])

# Finding best random effect

survSD.LG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "LG", ], family = binomial)


# Testing the part of variance explained by the RE

r.squaredGLMM(survSD.LG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "LG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding best fixed effect

survSD.LG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "LG", ], family = binomial)
survSD.LG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "LG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survSD.LG1, survSD.LG2, survSD.LG3, base = T) # survSD.LG2 and survSD.LG3 are in the 2 dAIC. survSD.LG2 is the simplest model.

summary(survSD.LG2) 
survSD.LG = survSD.LG2


# Plotting the observed data vs the predictions of the model

survSD.LG.obs = aggregate(surv ~ density + time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "SD"), ], FUN = mean)

survSD.LG.pred = survSD.LG.obs[, c(1, 2)]
survSD.LG.pred$pred = predict(survSD.LG, newdata = survSD.LG.pred, type = "response")


# Obs vs. Pred

plot(survSD.LG.obs$surv, survSD.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SD) survival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survSD.LG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.survSD.LG <- bootMer(survSD.LG, FUN = function(x) predict(x, survSD.LG.new.data, re.form = NA), 
  nsim = 100)
survSD.LG.new.data$lwr = inv.logit(apply(bootstrap.survSD.LG$t, 2, quantile, 0.025, na.rm = T))
survSD.LG.new.data$upr = inv.logit(apply(bootstrap.survSD.LG$t, 2, quantile, 0.975, na.rm = T))
survSD.LG.new.data$pred = predict(survSD.LG, newdata = survSD.LG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survSD.LG = ggplot(survSD.LG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Seedling survival - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_SDSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSD.LG

dev.off()


## 2.1.2. Little grazed (LG) - Juveniles (J) ----
# ------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "LG" & data.grazing$stage == "J")])


# Finding the best random effect

survJ.LG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "LG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(survJ.LG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "LG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survJ.LG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "LG", ], family = binomial)
survJ.LG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "LG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survJ.LG1, survJ.LG2, survJ.LG3, base = T, nobs = nrow(data.grazing[data.grazing$stage == "J" & data.grazing$LS == "LG", ])) # survJ.LG2 and survJ.LG3 are in 2 dAIC. survJ.LG2 is the simplest model.

summary(survJ.LG2) 
survJ.LG = survJ.LG2


# Plotting the observed data vs the predictions of the model
survJ.LG.obs = aggregate(surv ~ density + time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "J"), ], FUN = mean)

survJ.LG.pred = survJ.LG.obs[, c(1, 2)]
survJ.LG.pred$pred = predict(survJ.LG, newdata = survJ.LG.pred, type = "response")


# Obs vs. Pred

plot(survJ.LG.obs$surv, survJ.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of juveniles (J) survival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survJ.LG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.survJ.LG <- bootMer(survJ.LG, FUN = function(x) predict(x, survJ.LG.new.data, re.form = NA), 
  nsim = 100)
survJ.LG.new.data$lwr = inv.logit(apply(bootstrap.survJ.LG$t, 2, quantile, 0.025, na.rm = T))
survJ.LG.new.data$upr = inv.logit(apply(bootstrap.survJ.LG$t, 2, quantile, 0.975, na.rm = T))
survJ.LG.new.data$pred = predict(survJ.LG, newdata = survJ.LG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survJ.LG = ggplot(survJ.LG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Juvenile survival - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_JSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survJ.LG

dev.off()


## 2.1.3. Little grazed (LG) - Small reproductives (SR) ----
# -----------------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")])


# Find the best random effect

survSR.LG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(survSR.LG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survSR.LG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)
survSR.LG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survSR.LG1, survSR.LG2, survSR.LG3, base = T) # survSR.LG1 and survSR.LG2 are in 2 dAIC. survSR.LG1 is the simplest model.

summary(survSR.LG1)
survSR.LG = survSR.LG1


# Plotting the observed data vs the predictions of the model

survSR.LG.obs = aggregate(surv ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "SR"), ], FUN = mean)

survSR.LG.pred = data.frame(time = survSR.LG.obs[, 1])
survSR.LG.pred$pred = predict(survSR.LG, newdata = survSR.LG.pred, type = "response")


# Obs vs. Pred

plot(survSR.LG.obs$surv, survSR.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SR) survival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survSR.LG.new.data = expand.grid(time = rownames(coef(survSR.LG)$time))


# 95% confidence intervals
bootstrap.survSR.LG <- bootMer(survSR.LG, FUN = function(x) predict(x, survSR.LG.new.data), 
  nsim = 100)
survSR.LG.new.data$lwr = inv.logit(apply(bootstrap.survSR.LG$t, 2, quantile, 0.025, na.rm = T))
survSR.LG.new.data$upr = inv.logit(apply(bootstrap.survSR.LG$t, 2, quantile, 0.975, na.rm = T))
survSR.LG.new.data$pred = predict(survSR.LG, newdata = survSR.LG.new.data, type = "response")


# Plotting the predictions

plot.survSR.LG = ggplot(survSR.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Survival probability") + 
 ggtitle("Small reproductive survival - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_SRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSR.LG

dev.off()


## 2.1.4. Little grazed (LG) - Large reproductives (LR) ----
# -----------------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")])


# Finding the best random effect

survLR.LG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(survLR.LG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survLR.LG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)
survLR.LG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survLR.LG1, survLR.LG2, survLR.LG3, base = T) # survLR.LG1 and survLR.LG2 are in 2 dAIC. survLR.LG1 is the simplest model.

summary(survLR.LG1) 
survLR.LG = survLR.LG1


# Plotting the observed data vs the predictions of the model

survLR.LG.obs = aggregate(surv ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "LR"), ], FUN = mean)

survLR.LG.pred = data.frame(time = survLR.LG.obs[, 1])
survLR.LG.pred$pred = predict(survLR.LG1, newdata = survLR.LG.pred, type = "response")


# Obs vs. Pred

plot(survLR.LG.obs$surv, survLR.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of large reproductives (LR) survival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survLR.LG.new.data = expand.grid(time = rownames(coef(survLR.LG)$time))


# 95% confidence intervals
bootstrap.survLR.LG <- bootMer(survLR.LG, FUN = function(x) predict(x, survLR.LG.new.data), 
  nsim = 100)
survLR.LG.new.data$lwr = inv.logit(apply(bootstrap.survLR.LG$t, 2, quantile, 0.025, na.rm = T))
survLR.LG.new.data$upr = inv.logit(apply(bootstrap.survLR.LG$t, 2, quantile, 0.975, na.rm = T))
survLR.LG.new.data$pred = predict(survLR.LG, newdata = survLR.LG.new.data, type = "response")


# Plotting the predictions

plot.survLR.LG = ggplot(survLR.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Survival probability") + 
 ggtitle("Large reproductive survival - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_LRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survLR.LG

dev.off()


## 2.1.5. Highly grazed (HG) - Seedlings (SD) ----
# -------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "HG" & data.grazing$stage == "SD")])


# Finding the best random effect

survSD.HG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "HG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(survSD.HG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "HG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survSD.HG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "HG", ], family = binomial)
survSD.HG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survSD.HG1, survSD.HG2, survSD.HG3, base = T) # survSD.HG3 is the best model

summary(survSD.HG3) 
survSD.HG = glmer(surv ~ density + I(density^2) + (1|time), data = data.grazing[data.grazing$stage == "SD" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Plotting the observed data vs the predictions of the model

survSD.HG.obs = aggregate(surv ~ density + time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "SD"), ], FUN = mean)

survSD.HG.pred = survSD.HG.obs[, c(1, 2)]
survSD.HG.pred$pred = predict(survSD.HG, newdata = survSD.HG.pred, type = "response")


# Obs vs. Pred

plot(survSD.HG.obs$surv, survSD.HG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SD) survival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survSD.HG.new.data = expand.grid(density = seq(10, 15))


# 95% confidence intervals
bootstrap.survSD.HG <- bootMer(survSD.HG, FUN = function(x) predict(x, survSD.HG.new.data, re.form = NA), 
  nsim = 100) # Bootstraps fail
survSD.HG.new.data$lwr = inv.logit(apply(bootstrap.survSD.HG$t, 2, quantile, 0.025, na.rm = T))
survSD.HG.new.data$upr = inv.logit(apply(bootstrap.survSD.HG$t, 2, quantile, 0.975, na.rm = T))
survSD.HG.new.data$pred = predict(survSD.HG, newdata = survSD.HG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survSD.HG = ggplot(survSD.HG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 #geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Seedling survival - highly grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_SDSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSD.HG

dev.off()


## 2.1.6. Highly grazed (HG) - Juveniles (J) ----
# ------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "HG" & data.grazing$stage == "J")])


# Finding the best random effect

survJ.HG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "HG", ], family = binomial)

# Testing the part of variance explained by the RE
r.squaredGLMM(survJ.HG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "HG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survJ.HG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "HG", ], family = binomial)
survJ.HG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "J" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survJ.HG1, survJ.HG2, survJ.HG3, base = T) # survJ.HG1 and survJ.HG2 are in 2 dAIC. survJ.HG1 is the simplest model.

summary(survJ.HG1)
survJ.HG = survJ.HG1


# Plotting the observed data vs the predictions of the model

survJ.HG.obs = aggregate(surv ~ time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "J"), ], FUN = mean)

survJ.HG.pred = data.frame(time = survJ.HG.obs[, 1])
survJ.HG.pred$pred = predict(survJ.HG1, newdata = survJ.HG.pred, type = "response")


# Obs vs. Pred

plot(survJ.HG.obs$surv, survJ.HG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (J) survival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survJ.HG.new.data = expand.grid(time = rownames(coef(survJ.HG)$time))


# 95% confidence intervals
bootstrap.survJ.HG <- bootMer(survJ.HG, FUN = function(x) predict(x, survJ.HG.new.data), 
  nsim = 100)
survJ.HG.new.data$lwr = inv.logit(apply(bootstrap.survJ.HG$t, 2, quantile, 0.025, na.rm = T))
survJ.HG.new.data$upr = inv.logit(apply(bootstrap.survJ.HG$t, 2, quantile, 0.975, na.rm = T))
survJ.HG.new.data$pred = predict(survJ.HG, newdata = survJ.HG.new.data, type = "response")


# Plotting the predictions

plot.survJ.HG = ggplot(survJ.HG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Survival probability") + 
 ggtitle("Juvenile survival - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_JSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survJ.HG

dev.off()


## 2.1.7. Highly grazed (HG) - Small reproductives (SR) ----
# -----------------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "HG" & data.grazing$stage == "SR")])


# Finding the best random effect

survSR.HG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(survSR.HG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survSR.HG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial)
survSR.HG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survSR.HG1, survSR.HG2, survSR.HG3, base = T) # survSR.HG2 and survSR.HG3 are in 2 dAIC. survSR.HG2 is the simplest model.

summary(survSR.HG2)
survSR.HG = survSR.HG2


# Plotting the observed data vs the predictions of the model

survSR.HG.obs = aggregate(surv ~ density + time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "SR"), ], FUN = mean)

survSR.HG.pred = survSR.HG.obs[, c(1, 2)]
survSR.HG.pred$pred = predict(survSR.HG, newdata = survSR.HG.pred, type = "response")


# Obs vs. Pred

plot(survSR.HG.obs$surv, survSR.HG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (SR) survival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survSR.HG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.survSR.HG <- bootMer(survSR.HG, FUN = function(x) predict(x, survSR.HG.new.data, re.form = NA), 
  nsim = 100)
survSR.HG.new.data$lwr = inv.logit(apply(bootstrap.survSR.HG$t, 2, quantile, 0.025, na.rm = T))
survSR.HG.new.data$upr = inv.logit(apply(bootstrap.survSR.HG$t, 2, quantile, 0.975, na.rm = T))
survSR.HG.new.data$pred = predict(survSR.HG, newdata = survSR.HG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survSR.HG = ggplot(survSR.HG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Survival probability") + 
 ggtitle("Small reproductive survival - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_SRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survSR.HG

dev.off()


## 2.1.8. Highly grazed (HG) - Large reproductives (LR) ----
# -----------------------------------------------------

table(data.grazing$surv[which(data.grazing$LS == "HG" & data.grazing$stage == "SR")])


# Finding the best random effect

survLR.HG1 = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(survLR.HG1, null = glmer(surv ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

survLR.HG2 = glmer(surv ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial)
survLR.HG3 = glmer(surv ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(survLR.HG1, survLR.HG2, survLR.HG3, base = T) # All three models are in 2 dAIC. The simplest model is survLR.HG1.

summary(survLR.HG1) 
survLR.HG = survLR.HG1


# Plotting the observed data vs the predictions of the model

survLR.HG.obs = aggregate(surv ~ time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "LR"), ], FUN = mean)

survLR.HG.pred = data.frame(time = survLR.HG.obs[, 1])
survLR.HG.pred$pred = predict(survLR.HG1, newdata = survLR.HG.pred, type = "response")


# Obs vs. Pred

plot(survLR.HG.obs$surv, survLR.HG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", main = "Observed vs. Predicted values of seedlings (LR) survival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
survLR.HG.new.data = expand.grid(time = rownames(coef(survLR.HG)$time))


# 95% confidence intervals
bootstrap.survLR.HG <- bootMer(survLR.HG, FUN = function(x) predict(x, survLR.HG.new.data), 
  nsim = 100)
survLR.HG.new.data$lwr = inv.logit(apply(bootstrap.survLR.HG$t, 2, quantile, 0.025, na.rm = T))
survLR.HG.new.data$upr = inv.logit(apply(bootstrap.survLR.HG$t, 2, quantile, 0.975, na.rm = T))
survLR.HG.new.data$pred = predict(survLR.HG, newdata = survLR.HG.new.data, type = "response")


# Plotting the predictions

plot.survLR.HG = ggplot(survLR.HG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Survival probability") + 
 ggtitle("Large reproductive survival - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_LRSurvival.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.survLR.HG

dev.off()


## 2.2. Transition from juvenile (J) to either small reproductive (SR) or large reproductive (LR) (transitionJ) ----
# -------------------------------------------------------------------------------------------------------------

## 2.2.1. Little grazed (LG) ----
# --------------------------

table(data.grazing$transitionJ[which(data.grazing$LS == "LG" & data.grazing$stage == "J")])


# Finding the best random effect

transJ.LG1 = glmer(transitionJ ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "J", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(transJ.LG1, null = glmer(transitionJ ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "J", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

transJ.LG2 = glmer(transitionJ ~ density + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "J", ], family = binomial)
transJ.LG3 = glmer(transitionJ ~ density + density2 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "J", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(transJ.LG1, transJ.LG2, transJ.LG3, base = T) # transJ.LG1 and transJ.LG3 are in 2 dAIC. transJ.LG1 is the simplest model.

summary(transJ.LG1)
transJ.LG = transJ.LG1


# Plotting the observed data vs the predictions of the model

transJ.LG.obs = aggregate(transitionJ ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "J"), ], FUN = mean)

transJ.LG.pred = data.frame(time = transJ.LG.obs[, 1])
transJ.LG.pred$pred = predict(transJ.LG, newdata = transJ.LG.pred, type = "response")


# Obs vs. Pred

plot(transJ.LG.obs$transitionJ, transJ.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed transitionJival probability", ylab = "Predicted transitionJival probability", main = "Observed vs. Predicted values of seedlings (J) transitionJival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
transJ.LG.new.data = expand.grid(time = rownames(coef(transJ.LG)$time))


# 95% confidence intervals
bootstrap.transJ.LG <- bootMer(transJ.LG, FUN = function(x) predict(x, transJ.LG.new.data), 
  nsim = 100)
transJ.LG.new.data$lwr = inv.logit(apply(bootstrap.transJ.LG$t, 2, quantile, 0.025, na.rm = T))
transJ.LG.new.data$upr = inv.logit(apply(bootstrap.transJ.LG$t, 2, quantile, 0.975, na.rm = T))
transJ.LG.new.data$pred = predict(transJ.LG, newdata = transJ.LG.new.data, type = "response")


# Plotting the predictions

plot.transJ.LG = ggplot(transJ.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Transition probability") + 
 ggtitle("J-LR transition - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_JTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transJ.LG

dev.off()


## 2.2.2. Highly grazed (HG) ----
# --------------------------

table(data.grazing$transitionJ[which(data.grazing$LS == "HG" & data.grazing$stage == "J")])


# Finding the best random effect

transJ.HG1 = glmer(transitionJ ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "J", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(transJ.HG1, null = glmer(transitionJ ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "J", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

transJ.HG2 = glmer(transitionJ ~ density + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "J", ], family = binomial)
transJ.HG3 = glmer(transitionJ ~ density + density2 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "J", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(transJ.HG1, transJ.HG2, transJ.HG3, base = T) # transJ.HG2 and transJ.HG3 are in 2 dAIC. transJ.HG2 is the simplest model.

summary(transJ.HG2)
transJ.HG = transJ.HG2


# Plotting the observed data vs the predictions of the model

transJ.HG.obs = aggregate(transitionJ ~ density + time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "J"), ], FUN = mean)

transJ.HG.pred = transJ.HG.obs[, c(1, 2)]
transJ.HG.pred$pred = predict(transJ.HG2, newdata = transJ.HG.pred, type = "response")


# Obs vs. Pred

plot(transJ.HG.obs$transitionJ[which(transJ.HG.obs$time == "2017")], transJ.HG.pred$pred[which(transJ.HG.pred$time == "2017")], pch = 16, ylim = c(0, 1), xlab = "Observed transitionJival probability", ylab = "Predicted transitionJival probability", main = "Observed vs. Predicted values of seedlings (J) transitionJival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
transJ.HG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.transJ.HG <- bootMer(transJ.HG, FUN = function(x) predict(x, transJ.HG.new.data, re.form = NA), 
  nsim = 100)
transJ.HG.new.data$lwr = inv.logit(apply(bootstrap.transJ.HG$t, 2, quantile, 0.025, na.rm = T))
transJ.HG.new.data$upr = inv.logit(apply(bootstrap.transJ.HG$t, 2, quantile, 0.975, na.rm = T))
transJ.HG.new.data$pred = predict(transJ.HG, newdata = transJ.HG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.transJ.HG = ggplot(transJ.HG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Transition probability") + 
 ggtitle("J-LR transition - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_JTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transJ.HG

dev.off()


## 2.3. Transition from small reproductive (SR) to either small reproductive (SR = stasis) or large reproductive (LR = growth) (transitionSR) ----
# -------------------------------------------------------------------------------------------------------------------------------------------

## 2.3.1. Little grazed (LG) ----
# --------------------------

table(data.grazing$transitionSR[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")])


# Finding the best random effect

transSR.LG1 = glmer(transitionSR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "SR", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(transSR.LG1, null = glmer(transitionSR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "SR", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

transSR.LG2 = glmer(transitionSR ~ density + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "SR", ], family = binomial)
transSR.LG3 = glmer(transitionSR ~ density + density2 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "SR", ], family = binomial)

AICctab(transSR.LG1, transSR.LG2, transSR.LG3, base = T) # transSR.LG1 and transSR.LG2 are in 2 dAIC. transSR.LG1 is the simplest model.

summary(transSR.LG1)
transSR.LG = transSR.LG1


# Plotting the observed data vs the predictions of the model

transSR.LG.obs = aggregate(transitionSR ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "SR"), ], FUN = mean)

transSR.LG.pred = data.frame(time = transSR.LG.obs[, 1])
transSR.LG.pred$pred = predict(transSR.LG1, newdata = transSR.LG.pred, type = "response")


# Obs vs. Pred

plot(transSR.LG.obs$transitionSR, transSR.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed transitionSRival probability", ylab = "Predicted transitionSRival probability", main = "Observed vs. Predicted values of seedlings (J) transitionSRival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
transSR.LG.new.data = expand.grid(time = rownames(coef(transSR.LG)$time))


# 95% confidence intervals
bootstrap.transSR.LG <- bootMer(transSR.LG, FUN = function(x) predict(x, transSR.LG.new.data), 
  nsim = 100)
transSR.LG.new.data$lwr = inv.logit(apply(bootstrap.transSR.LG$t, 2, quantile, 0.025, na.rm = T))
transSR.LG.new.data$upr = inv.logit(apply(bootstrap.transSR.LG$t, 2, quantile, 0.975, na.rm = T))
transSR.LG.new.data$pred = predict(transSR.LG, newdata = transSR.LG.new.data, type = "response")

transSR.LG.pred = transSR.LG.new.data[, c(1, 4)]
transSR.LG.pred$time = as.character(transSR.LG.pred$time)
transSR.LG.pred = rbind(transSR.LG.pred, c(2017, 0))


# Plotting the predictions

plot.transSR.LG = ggplot(transSR.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Transition probability") + 
 ggtitle("SR-LR transition - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_SRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transSR.LG

dev.off()


## 2.3.2. Highly grazed (HG) ----
# --------------------------

table(data.grazing$transitionSR[which(data.grazing$LS == "HG" & data.grazing$stage == "SR")])


# Finding the best random effect

transSR.HG1 = glmer(transitionSR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "SR", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(transSR.HG1, null = glmer(transitionSR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "SR", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

transSR.HG2 = glmer(transitionSR ~ density + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "SR", ], family = binomial)
transSR.HG3 = glmer(transitionSR ~ density + density2 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "SR", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(transSR.HG1, transSR.HG2, transSR.HG3, base = T) # transSR.HG1 and transSR.HG3 are in 2 dAIC. transSR.HG1 is the simplest model.

summary(transSR.HG1)
transSR.HG = transSR.HG1


# Plotting the observed data vs the predictions of the model

transSR.HG.obs = aggregate(transitionSR ~ time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "SR"), ], FUN = mean)

transSR.HG.pred = data.frame(time = transSR.HG.obs[, 1])
transSR.HG.pred$pred = predict(transSR.HG1, newdata = transSR.HG.pred, type = "response")


# Obs vs. Pred

plot(transSR.HG.obs$transitionSR, transSR.HG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed transitionSRival probability", ylab = "Predicted transitionSRival probability", main = "Observed vs. Predicted values of seedlings (J) transitionSRival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
transSR.HG.new.data = expand.grid(time = rownames(coef(transSR.HG)$time))


# 95% confidence intervals
bootstrap.transSR.HG <- bootMer(transSR.HG, FUN = function(x) predict(x, transSR.HG.new.data), 
  nsim = 100)
transSR.HG.new.data$lwr = inv.logit(apply(bootstrap.transSR.HG$t, 2, quantile, 0.025, na.rm = T))
transSR.HG.new.data$upr = inv.logit(apply(bootstrap.transSR.HG$t, 2, quantile, 0.975, na.rm = T))
transSR.HG.new.data$pred = predict(transSR.HG, newdata = transSR.HG.new.data, type = "response")


# Plotting the predictions

plot.transSR.HG = ggplot(transSR.HG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Transition probability") + 
 ggtitle("SR-LR transition - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_SRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transSR.HG

dev.off()


## 2.4. Transition from large reproductive (LR) to either small reproductive (SR = shrinkage) or large reproductive (LR = stasis) (transitionLR) ----
# ----------------------------------------------------------------------------------------------------------------------------------------------

## 2.4.1. Little grazed (LG) ----
# --------------------------

table(data.grazing$transitionLR[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")])


# Finding the best random effect

transLR.LG1 = glmer(transitionLR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "LR", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(transLR.LG1, null = glmer(transitionLR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "LR", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

transLR.LG2 = glmer(transitionLR ~ density + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "LR", ], family = binomial)
transLR.LG3 = glmer(transitionLR ~ density + density2 + (1|time), data = data.grazing[data.grazing$LS == "LG" & data.grazing$stage == "LR", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(transLR.LG1, transLR.LG2, transLR.LG3, base = T) # transLR.LG1 and transLR.LG2 are in 2 dAIC. transLR.LG1 is the simplest model.

summary(transLR.LG1)
transLR.LG = transLR.LG1


# Plotting the observed data vs the predictions of the model

transLR.LG.obs = aggregate(transitionLR ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "LR"), ], FUN = mean)

transLR.LG.pred = data.frame(time = transLR.LG.obs[, 1])
transLR.LG.pred$pred = predict(transLR.LG, newdata = transLR.LG.pred, type = "response")


# Obs vs. Pred

plot(transLR.LG.obs$transitionLR, transLR.LG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed transitionLRival probability", ylab = "Predicted transitionLRival probability", main = "Observed vs. Predicted values of seedlings (J) transitionLRival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
transLR.LG.new.data = expand.grid(time = rownames(coef(transLR.LG)$time))


# 95% confidence intervals
bootstrap.transLR.LG <- bootMer(transLR.LG, FUN = function(x) predict(x, transLR.LG.new.data), 
  nsim = 100)
transLR.LG.new.data$lwr = inv.logit(apply(bootstrap.transLR.LG$t, 2, quantile, 0.025, na.rm = T))
transLR.LG.new.data$upr = inv.logit(apply(bootstrap.transLR.LG$t, 2, quantile, 0.975, na.rm = T))
transLR.LG.new.data$pred = predict(transLR.LG, newdata = transLR.LG.new.data, type = "response")


# Plotting the predictions

plot.transLR.LG = ggplot(transLR.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Transition probability") + 
 ggtitle("LR-SR transition - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_LRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transLR.LG

dev.off()


## 2.4.2. Highly grazed (HG) ----
# --------------------------

table(data.grazing$transitionLR[which(data.grazing$LS == "HG" & data.grazing$stage == "LR")])


# Finding the best random effect

transLR.HG1 = glmer(transitionLR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "LR", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Testing the part of variance explained by the RE
r.squaredGLMM(transLR.HG1, null = glmer(transitionLR ~ 1 + (1|time), data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "LR", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))) # The RE does not explain any part of the variance. We can thus model the vital rate with a GLM.


# Finding the best fixed effect

transLR.HG2 = glm(transitionLR ~ density, data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "LR", ], family = binomial)
transLR.HG3 = glm(transitionLR ~ density + density2, data = data.grazing[data.grazing$LS == "HG" & data.grazing$stage == "LR", ], family = binomial)

AICctab(transLR.HG1, transLR.HG2, transLR.HG3, base = T) # transLR.HG2 and transLR.HG3 are in 2 dAIC. transLR.HG2 is the simplest model.

summary(transLR.HG2)
transLR.HG = transLR.HG2


# Plotting the observed data vs the predictions of the model

transLR.HG.obs = aggregate(transitionLR ~ density + time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "LR"), ], FUN = mean)

transLR.HG.pred = transLR.HG.obs[, c(1, 2)]
transLR.HG.pred$pred = predict(transLR.HG, newdata = transLR.HG.pred, type = "response")


# Obs vs. Pred

plot(transLR.HG.obs$transitionLR[which(transLR.HG.obs$time == 2011)], transLR.HG.pred$pred[which(transLR.HG.pred$time == 2011)], pch = 16, ylim = c(0, 1), xlab = "Observed transitionLRival probability", ylab = "Predicted transitionLRival probability", main = "Observed vs. Predicted values of seedlings (J) transitionLRival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
transLR.HG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
transLR.HG.predict = predict(transLR.HG, newdata = transLR.HG.new.data, se.fit = T)
transLR.HG.new.data$lwr = inv.logit(transLR.HG.predict$fit - (1.96 * transLR.HG.predict$se.fit))
transLR.HG.new.data$upr = inv.logit(transLR.HG.predict$fit + (1.96 * transLR.HG.predict$se.fit))
transLR.HG.new.data$pred = inv.logit(transLR.HG.predict$fit)


# Plotting the predictions

plot.transLR.HG = ggplot(transLR.HG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Transition probability") + 
 ggtitle("LR-SR transition - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_LRTransition.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.transLR.HG

dev.off()


## 2.5. Probability of flowering (fl) ----
# -----------------------------------

## 2.5.1. Little grazed (LG) - Small reproductive (SR) ----
# ----------------------------------------------------

table(data.grazing$fl[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")])


# Finding the best random effect

floweringprobSR.LG1 = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(floweringprobSR.LG1, null = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

floweringprobSR.LG2 = glmer(fl ~ density + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)
floweringprobSR.LG3 = glmer(fl ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)

AICctab(floweringprobSR.LG1, floweringprobSR.LG2, floweringprobSR.LG3, base = T) # The best model is floweringprobSR.LG3.

summary(floweringprobSR.LG3) 
floweringprobSR.LG = glmer(fl ~ density + I(density^2) + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "LG", ], family = binomial)


# Plotting the observed data vs the predictions of the model
floweringprobSR.LG.obs = aggregate(fl ~ density + time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "SR"), ], FUN = mean)

floweringprobSR.LG.pred = floweringprobSR.LG.obs[, c(1, 2)]
floweringprobSR.LG.pred$pred = predict(floweringprobSR.LG, newdata = floweringprobSR.LG.pred, type = "response")


# Obs vs. Pred

plot(floweringprobSR.LG.obs$fl[which(floweringprobSR.LG.obs$time == 2011)], floweringprobSR.LG.pred$pred[which(floweringprobSR.LG.pred$time == 2011)], pch = 16, ylim = c(0, 1), xlab = "Observed flival probability", ylab = "Predicted flival probability", main = "Observed vs. Predicted values of seedlings (J) flival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
floweringprobSR.LG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.floweringprobSR.LG <- bootMer(floweringprobSR.LG, FUN = function(x) predict(x, floweringprobSR.LG.new.data, re.form = NA), 
   nsim = 100)
floweringprobSR.LG.new.data$lwr = inv.logit(apply(bootstrap.floweringprobSR.LG$t, 2, quantile, 0.025, na.rm = T))
floweringprobSR.LG.new.data$upr = inv.logit(apply(bootstrap.floweringprobSR.LG$t, 2, quantile, 0.975, na.rm = T))
floweringprobSR.LG.new.data$pred = predict(floweringprobSR.LG, newdata = floweringprobSR.LG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.floweringprobSR.LG = ggplot(floweringprobSR.LG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Flowering probability") + 
 ggtitle("Small reproductive flowering probability - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_LRFlowProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobSR.LG

dev.off()


## 2.5.2. Little grazed (LG) -  Large reproductive (LR) ----
# -----------------------------------------------------

table(data.grazing$fl[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")])


# Finding the best random effect

floweringprobLR.LG1 = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(floweringprobLR.LG1, null = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

floweringprobLR.LG2 = glmer(fl ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)
floweringprobLR.LG3 = glmer(fl ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = binomial)

AICctab(floweringprobLR.LG1, floweringprobLR.LG2, floweringprobLR.LG3, base = T) #The best model is floweringprobLR.LG2.

summary(floweringprobLR.LG2) 
floweringprobLR.LG = floweringprobLR.LG2


# Plotting the observed data vs the predictions of the model
floweringprobLR.LG.obs = aggregate(fl ~ density + time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "LR"), ], FUN = mean)

floweringprobLR.LG.pred = floweringprobLR.LG.obs[, c(1, 2)]
floweringprobLR.LG.pred$pred = predict(floweringprobLR.LG2, newdata = floweringprobLR.LG.pred, type = "response")


# Obs vs. Pred

plot(floweringprobLR.LG.obs$fl[which(floweringprobLR.LG.obs$time == 2011)], floweringprobLR.LG.pred$pred[which(floweringprobLR.LG.pred$time == 2011)], pch = 16, ylim = c(0, 1), xlab = "Observed flival probability", ylab = "Predicted flival probability", main = "Observed vs. Predicted values of seedlings (J) flival - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
floweringprobLR.LG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.floweringprobLR.LG <- bootMer(floweringprobLR.LG, FUN = function(x) predict(x, floweringprobLR.LG.new.data, re.form = NA), 
   nsim = 100)
floweringprobLR.LG.new.data$lwr = inv.logit(apply(bootstrap.floweringprobLR.LG$t, 2, quantile, 0.025, na.rm = T))
floweringprobLR.LG.new.data$upr = inv.logit(apply(bootstrap.floweringprobLR.LG$t, 2, quantile, 0.975, na.rm = T))
floweringprobLR.LG.new.data$pred = predict(floweringprobLR.LG, newdata = floweringprobLR.LG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.floweringprobLR.LG = ggplot(floweringprobLR.LG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Flowering probability") + 
 ggtitle("Large reproductive flowering probability - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_LRFlowProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobLR.LG

dev.off()


## 2.5.3. Highly grazed (HG) - Small reproductive (SR) ----
# ----------------------------------------------------

table(data.grazing$fl[which(data.grazing$LS == "HG" & data.grazing$stage == "SR")])


# Finding the best random effect

floweringprobSR.HG1 = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(floweringprobSR.HG1, null = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

floweringprobSR.HG2 = glmer(fl ~ density + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial)
floweringprobSR.HG3 = glmer(fl ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(floweringprobSR.HG1, floweringprobSR.HG2, floweringprobSR.HG3, base = T) # All three models are in 2 dAIC. floweringprobSR.HG1 is the simplest model.

summary(floweringprobSR.HG1)
floweringprobSR.HG = floweringprobSR.HG1


# Plotting the observed data vs the predictions of the model

floweringprobSR.HG.obs = aggregate(fl ~ time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "SR"), ], FUN = mean)

floweringprobSR.HG.pred = data.frame(time = floweringprobSR.HG.obs[, 1])
floweringprobSR.HG.pred$pred = predict(floweringprobSR.HG, newdata = floweringprobSR.HG.pred, type = "response")


# Obs vs. Pred

plot(floweringprobSR.HG.obs$fl, floweringprobSR.HG.pred$pred, pch = 16, ylim = c(0, 1), xlab = "Observed flival probability", ylab = "Predicted flival probability", main = "Observed vs. Predicted values of seedlings (J) flival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
floweringprobSR.HG.new.data = expand.grid(time = rownames(coef(floweringprobSR.HG)$time))


# 95% confidence intervals
bootstrap.floweringprobSR.HG <- bootMer(floweringprobSR.HG, FUN = function(x) predict(x, floweringprobSR.HG.new.data), 
   nsim = 100)
floweringprobSR.HG.new.data$lwr = inv.logit(apply(bootstrap.floweringprobSR.HG$t, 2, quantile, 0.025, na.rm = T))
floweringprobSR.HG.new.data$upr = inv.logit(apply(bootstrap.floweringprobSR.HG$t, 2, quantile, 0.975, na.rm = T))
floweringprobSR.HG.new.data$pred = predict(floweringprobSR.HG, newdata = floweringprobSR.HG.new.data, type = "response")


# Plotting the predictions 

plot.floweringprobSR.HG = ggplot(floweringprobSR.HG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Flowering probability") + 
 ggtitle("Small reproductive flowering probability - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_SRFlowProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobSR.HG

dev.off()


## 2.5.4. Highly grazed (HG) - Large reproductive (LR) ----
# ----------------------------------------------------

table(data.grazing$fl[which(data.grazing$LS == "HG" & data.grazing$stage == "LR")])


# Finding the best random effect

floweringprobLR.HG1 = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial)


# Testing the part of variance explained by the RE
r.squaredGLMM(floweringprobLR.HG1, null = glmer(fl ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

floweringprobLR.HG2 = glmer(fl ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial)
floweringprobLR.HG3 = glmer(fl ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(floweringprobLR.HG1, floweringprobLR.HG2, floweringprobLR.HG3, base = T) # floweringprobLR.HG2 and floweringprobLR.HG3 are in 2 dAIC. floweringprobLR.HG2 is the simplest model.

summary(floweringprobLR.HG2) 
floweringprobLR.HG = floweringprobLR.HG2


# Plotting the observed data vs the predictions of the model

floweringprobLR.HG.obs = aggregate(fl ~ density + time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "LR"), ], FUN = mean)

floweringprobLR.HG.pred = floweringprobLR.HG.obs[, c(1, 2)]
floweringprobLR.HG.pred$pred = predict(floweringprobLR.HG, newdata = floweringprobLR.HG.pred, type = "response")


# Obs vs. Pred

plot(floweringprobLR.HG.obs$fl[which(floweringprobLR.HG.obs$time == 2012)], floweringprobLR.HG.pred$pred[which(floweringprobLR.HG.pred$time == 2012)], pch = 16, ylim = c(0, 1), xlab = "Observed flival probability", ylab = "Predicted flival probability", main = "Observed vs. Predicted values of seedlings (J) flival - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
floweringprobLR.HG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# 95% confidence intervals
bootstrap.floweringprobLR.HG <- bootMer(floweringprobLR.HG, FUN = function(x) predict(x, floweringprobLR.HG.new.data, re.form = NA), 
   nsim = 100)
floweringprobLR.HG.new.data$lwr = inv.logit(apply(bootstrap.floweringprobLR.HG$t, 2, quantile, 0.025, na.rm = T))
floweringprobLR.HG.new.data$upr = inv.logit(apply(bootstrap.floweringprobLR.HG$t, 2, quantile, 0.975, na.rm = T))
floweringprobLR.HG.new.data$pred = predict(floweringprobLR.HG, newdata = floweringprobLR.HG.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.floweringprobLR.HG = ggplot(floweringprobLR.HG.new.data, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Flowering probability") + 
 ggtitle("Large reproductive flowering probability - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_LRFlowProb.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.floweringprobLR.HG

dev.off()


## 2.6. Number of flowering stalks (fs) ----
# -------------------------------------

## 2.6.1. Little grazed (LG) - Small reproductive (SR) ----
# ----------------------------------------------------

table(data.grazing$fs[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")], data.grazing$time[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")]) # There is not much variation, we can take the mean per year.

nbfsSR.LG = aggregate(fs ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "SR"), ], FUN = mean, na.rm = T)
nbfsSR.LG$time = as.numeric(as.character(nbfsSR.LG$time))

nbfsSR.LG = rbind(nbfsSR.LG, c(2016, mean(nbfsSR.LG$fs)))


# Plotting the predictions

plot.nbfsSR.LG = ggplot(nbfsSR.LG, aes(x = as.numeric(as.character(time)), y = fs)) + 
 geom_line(size = 3) + 
 xlab("Years") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Small reproductive number of flowering stalks - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_SRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsSR.LG

dev.off()


## 2.6.2. Little grazed (LG) - Large reproductive (LR) ----
# ----------------------------------------------------

table(data.grazing$fs[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")], data.grazing$time[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")]) 


# Finding the best random effect 

nbfsLR.LG1 = glmer(fs ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson)


# Testing the part of variance explained by the RE
r.squaredGLMM(nbfsLR.LG1, null = glmer(fs ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

nbfsLR.LG2 = glmer(fs ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson)
nbfsLR.LG3 = glmer(fs ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(nbfsLR.LG1, nbfsLR.LG2, nbfsLR.LG3, base = T) # All three models are in 2 dAIC. nbfsLR.LG1 is the simplest model.


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(nbfsLR.LG1)

nbfsLR.LG.quasiP = MASS::glmmPQL(fs ~ 1, random =  ~ 1|time, family = poisson, data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ])

summary(nbfsLR.LG.quasiP)
nbfsLR.LG = nbfsLR.LG.quasiP


# Plotting the observed data vs the predictions of the model

nbfsLR.LG.obs = aggregate(fs ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "LR"), ], FUN = mean)

nbfsLR.LG.pred = data.frame(time = nbfsLR.LG.obs[, 1])
nbfsLR.LG.pred$pred = predict(nbfsLR.LG, newdata = nbfsLR.LG.pred, type = "response")


# Obs vs. Pred

plot(nbfsLR.LG.obs$fs, nbfsLR.LG.pred$pred, pch = 16, xlab = "Observed LR number of flowering stalks", ylab = "Predicted LR number of flowering stalks", main = "Observed vs. Predicted values of LR number of flowering stalks - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
nbfsLR.LG.new.data = expand.grid(time = rownames(coef(nbfsLR.LG)))
nbfsLR.LG.new.data$lwr = NA
nbfsLR.LG.new.data$upr = NA
nbfsLR.LG.new.data$pred = NA

# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function 

nbfsLR.LG.new.data$lwr = apply(nbfsLR.LG.new.data, 1, FUN = function(x) easyPredCI(nbfsLR.LG, newdata = data.frame(time = as.numeric(x[1])))[1, 1])

nbfsLR.LG.new.data$upr = apply(nbfsLR.LG.new.data, 1, FUN = function(x) easyPredCI(nbfsLR.LG, newdata = data.frame(time = as.numeric(x[1])))[1, 2])

nbfsLR.LG.new.data$pred = apply(nbfsLR.LG.new.data, 1, FUN = function(x) predict(nbfsLR.LG, newdata = data.frame(time = as.numeric(x[1])), type = "response"))


nbfsLR.LG.pred = nbfsLR.LG.new.data
nbfsLR.LG.pred$time = as.numeric(as.character(nbfsLR.LG.pred$time))
nbfsLR.LG.pred = rbind(nbfsLR.LG.pred, c(2016, apply(nbfsLR.LG.pred[, c(2:4)], 2, mean)))


# Plotting the predictions

plot.nbfsLR.LG = ggplot(nbfsLR.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Large reproductive number of flowering stalks - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_LRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsLR.LG

dev.off()


## 2.6.3. Highly grazed (HG) - Small reproductive (SR) ----
# ---------------------------------------------------

table(data.grazing$fs[which(data.grazing$LS == "HG" & data.grazing$stage == "SR")], data.grazing$time[which(data.grazing$LS == "HG" & data.grazing$stage == "SR")]) # There is not much variation, we can take the mean per year.

nbfsSR.HG = aggregate(fs ~ time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "SR"), ], FUN = mean, na.rm = T)
nbfsSR.HG$time = as.numeric(as.character(nbfsSR.HG$time))


# Plotting the predictions

plot.nbfsSR.HG = ggplot(nbfsSR.HG, aes(x = as.numeric(as.character(time)), y = fs)) + 
 geom_line(size = 3) + 
 xlab("Years") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Small reproductive number of flowering stalks - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_SRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsSR.HG

dev.off()


## 2.6.4. Highly grazed (HG) - Large reproductive (LR) ----
# ----------------------------------------------------

table(data.grazing$fs[which(data.grazing$LS == "HG" & data.grazing$stage == "LR")], data.grazing$time[which(data.grazing$LS == "HG" & data.grazing$stage == "LR")]) # There is not much variation, we can take the mean per year.

nbfsLR.HG = aggregate(fs ~ time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "LR"), ], FUN = mean, na.rm = T)
nbfsLR.HG$time = as.numeric(as.character(nbfsLR.HG$time))


# Plotting the predictions

plot.nbfsLR.HG = ggplot(nbfsLR.HG, aes(x = as.numeric(as.character(time)), y = fs)) + 
 geom_line(size = 3) + 
 xlab("Years") + 
 ylab("Number of flowering stalks") + 
 ggtitle("Large reproductive number of flowering stalks - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_LRNbFlowStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfsLR.HG

dev.off()


## 2.7. Number of flowers per stalk (fps) ----
# ---------------------------------------

## 2.7.1. Little grazing (LG) - Small reproductive (SR) ----
# -----------------------------------------------------

table(data.grazing$fps[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")], data.grazing$time[which(data.grazing$LS == "LG" & data.grazing$stage == "SR")]) # There is not enough variation in most years. We take the mean per year.

nbfpsSR.LG = aggregate(fps ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "SR"), ], FUN = mean, na.rm = T)
nbfpsSR.LG$time = as.numeric(as.character(nbfpsSR.LG$time))

nbfpsSR.LG = rbind(nbfpsSR.LG, c(2016, mean(nbfpsSR.LG$fps)))


# Plotting the predictions 

plot.nbfpsSR.LG = ggplot(nbfpsSR.LG, aes(x = as.numeric(as.character(time)), y = fps)) + 
 geom_line(size = 3) + 
 xlab("Years") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Small reproductive number of flowers per stalk - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_SRNbFlowPerStalk.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsSR.LG

dev.off()


## 2.7.2. Little grazed (LG) - Large reproductive (LR) ----
# ----------------------------------------------------

table(data.grazing$fps[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")], data.grazing$time[which(data.grazing$LS == "LG" & data.grazing$stage == "LR")])


# Finding the best random effect

nbfpsLR.LG1 = glmer(fps ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson)


# Testing the part of variance explained by the RE
r.squaredGLMM(nbfpsLR.LG1, null = glmer(fps ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

nbfpsLR.LG2 = glmer(fps ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson)
nbfpsLR.LG3 = glmer(fps ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ], family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(nbfpsLR.LG1, nbfpsLR.LG2, nbfpsLR.LG3, base = T) # nbfpsLR.LG1 and nbfpsLR.LG2 are in 2 dAIC. nbfpsLR.LG1 is the simplest model.


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(nbfpsLR.LG1)

nbfpsLR.LG.quasiP = MASS::glmmPQL(fps ~ 1, random =  ~ 1|time, family = poisson, data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "LG", ])

summary(nbfpsLR.LG.quasiP)
nbfpsLR.LG = nbfpsLR.LG.quasiP


# Plotting the observed data vs the predictions of the model

nbfpsLR.LG.obs = aggregate(fps ~ time, data = data.grazing[which(data.grazing$LS == "LG" & data.grazing$stage == "LR"), ], FUN = mean)

nbfpsLR.LG.pred = data.frame(time = nbfpsLR.LG.obs[, 1])
nbfpsLR.LG.pred$pred = predict(nbfpsLR.LG, newdata = nbfpsLR.LG.pred, type = "response")


# Obs vs. Pred

plot(nbfpsLR.LG.obs$fps, nbfpsLR.LG.pred$pred, pch = 16, xlab = "Observed LR number of flowers per stalks", ylab = "Predicted LR number of flowers per stalks", main = "Observed vs. Predicted values of LR number of flowers per stalk - LG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
nbfpsLR.LG.new.data = expand.grid(time = rownames(coef(nbfpsLR.LG)))
nbfpsLR.LG.new.data$lwr = NA
nbfpsLR.LG.new.data$upr = NA
nbfpsLR.LG.new.data$pred = NA

# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function from Prof. Marc Girondot (https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

nbfpsLR.LG.new.data$lwr = apply(nbfpsLR.LG.new.data, 1, FUN = function(x) easyPredCI(nbfpsLR.LG, newdata = data.frame(time = as.numeric(x[1])))[1, 1])

nbfpsLR.LG.new.data$upr = apply(nbfpsLR.LG.new.data, 1, FUN = function(x) easyPredCI(nbfpsLR.LG, newdata = data.frame(time = as.numeric(x[1])))[1, 2])

nbfpsLR.LG.new.data$pred = apply(nbfpsLR.LG.new.data, 1, FUN = function(x) predict(nbfpsLR.LG, newdata = data.frame(time = as.numeric(x[1])), type = "response"))


nbfpsLR.LG.pred = nbfpsLR.LG.new.data
nbfpsLR.LG.pred$time = as.numeric(as.character(nbfpsLR.LG.pred$time))
nbfpsLR.LG.pred = rbind(nbfpsLR.LG.pred, c(2016, apply(nbfpsLR.LG.pred[, c(2:4)], 2, mean)))

# Plotting the predictions

plot.nbfpsLR.LG = ggplot(nbfpsLR.LG.new.data, aes(x = as.numeric(as.character(time)), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Years") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Large reproductive number of flowers per stalks - Little grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "LG_LRNbFlowPerStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsLR.LG

dev.off()


## 2.7.3. High grazing (HG) - Small reproductive (SR) ----
# ---------------------------------------------------

table(data.grazing$fps[data.grazing$LS == "HG" & data.grazing$stage == "SR"], data.grazing$time[data.grazing$LS == "HG" & data.grazing$stage == "SR"])


# Finding the best random effect

nbfpsSR.HG1 = glmer(fps ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = poisson)


# Testing the part of variance explained by the RE
r.squaredGLMM(nbfpsSR.HG1, null = glmer(fps ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = poisson)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

nbfpsSR.HG2 = glmer(fps ~ density + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = poisson)
nbfpsSR.HG3 = glmer(fps ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ], family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(nbfpsSR.HG1, nbfpsSR.HG2, nbfpsSR.HG3, base = T) # The best model is nbfpsSR.HG2.


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(nbfpsSR.HG2)

nbfpsSR.HG.quasiP = MASS::glmmPQL(fps ~ density, random =  ~ 1|time, family = poisson, data = data.grazing[data.grazing$stage == "SR" & data.grazing$LS == "HG", ])

summary(nbfpsSR.HG.quasiP)
nbfpsSR.HG = nbfpsSR.HG.quasiP


# Plotting the observed data vs the predictions of the model

nbfpsSR.HG.obs = aggregate(fps ~ density + time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "SR"), ], FUN = mean)

nbfpsSR.HG.pred = nbfpsSR.HG.obs[, c(1, 2)]
nbfpsSR.HG.pred$pred = predict(nbfpsSR.HG, newdata = nbfpsSR.HG.pred, type = "response")


# Obs vs. Pred

plot(nbfpsSR.HG.obs$fps[which(nbfpsSR.HG.obs$time == 2012)], nbfpsSR.HG.pred$pred[which(nbfpsSR.HG.pred$time == 2012)], pch = 16, xlab = "Observed SR number of flowers per stalk", ylab = "Predicted SR number of flowers per stalk", main = "Observed vs. Predicted values of SR number of flowers per stalk - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
nbfpsSR.HG.new.data = expand.grid(density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))


# Computing the predictions of the model and the 95 % CI
nbfpsSR.HG.new.data = expand.grid(time = rownames(coef(nbfpsSR.HG)),
                                  density = seq(min(data.grazing$density, na.rm = T), max(data.grazing$density, na.rm = T)))
nbfpsSR.HG.new.data$lwr = NA
nbfpsSR.HG.new.data$upr = NA
nbfpsSR.HG.new.data$pred = NA


# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function from Prof. Marc Girondot (https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

nbfpsSR.HG.new.data$lwr = apply(nbfpsSR.HG.new.data, 1, FUN = function(x) easyPredCI(nbfpsSR.HG, newdata = data.frame(time = as.numeric(x[1]), density = as.numeric(x[2])))[1, 1])

nbfpsSR.HG.new.data$upr = apply(nbfpsSR.HG.new.data, 1, FUN = function(x) easyPredCI(nbfpsSR.HG, newdata = data.frame(time = as.numeric(x[1]), density = as.numeric(x[2])))[1, 2])

nbfpsSR.HG.new.data$pred = apply(nbfpsSR.HG.new.data, 1, FUN = function(x) predict(nbfpsSR.HG, newdata = data.frame(time = as.numeric(x[1]), density = as.numeric(x[2])), type = "response"))

nbfpsSR.HG.RE.avg = data.frame(density = aggregate(pred ~ density, data = nbfpsSR.HG.new.data, mean)$density,
                                  pred = aggregate(pred ~ density, data = nbfpsSR.HG.new.data, mean)$pred,
                                  lwr = aggregate(lwr ~ density, data = nbfpsSR.HG.new.data, mean)$lwr,
                                  upr = aggregate(upr ~ density, data = nbfpsSR.HG.new.data, mean)$upr)


# Plotting the predictions

plot.nbfpsSR.HG = ggplot(nbfpsSR.HG.RE.avg, aes(x = density, y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Total abundance of dewy pines") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Small reproductive number of flowers per stalks - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_SRNbFlowPerStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsSR.HG

dev.off()


## 2.7.4. Highly grazed (HG) - Large reproductive (LR) ----
# ----------------------------------------------------

table(data.grazing$fps[data.grazing$LS == "HG" & data.grazing$stage == "LR"], data.grazing$time[data.grazing$LS == "HG" & data.grazing$stage == "LR"])


# Finding the best random effect

nbfpsLR.HG1 = glmer(fps ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = poisson)


# Testing the part of variance explained by the RE
r.squaredGLMM(nbfpsLR.HG1, null = glmer(fps ~ 1 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = poisson)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.


# Finding the best fixed effect

nbfpsLR.HG2 = glmer(fps ~ density + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = poisson)
nbfpsLR.HG3 = glmer(fps ~ density + density2 + (1|time), data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ], family = poisson, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))

AICctab(nbfpsLR.HG1, nbfpsLR.HG2, nbfpsLR.HG3, base = T) # nbfpsLR.HG1 and nbfpsLR.HG3 are in 2 dAIC. nbfpsLR.HG1 is the simplest model. 


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(nbfpsLR.HG1)

nbfpsLR.HG.quasiP = MASS::glmmPQL(fps ~ 1, random =  ~ 1|time, family = poisson, data = data.grazing[data.grazing$stage == "LR" & data.grazing$LS == "HG", ])

summary(nbfpsLR.HG.quasiP)
nbfpsLR.HG = nbfpsLR.HG.quasiP


# Plotting the observed data vs the predictions of the model

nbfpsLR.HG.obs = aggregate(fps ~  time, data = data.grazing[which(data.grazing$LS == "HG" & data.grazing$stage == "LR"), ], FUN = mean)

nbfpsLR.HG.pred = data.frame(time = nbfpsLR.HG.obs[, 1])
nbfpsLR.HG.pred$pred = predict(nbfpsLR.HG, newdata = nbfpsLR.HG.pred, type = "response")


# Obs vs. Pred

plot(nbfpsLR.HG.obs$fps, nbfpsLR.HG.pred$pred, pch = 16, xlab = "Observed LR number of flowers per stalks", ylab = "Predicted LR number of flowers per stalks", main = "Observed vs. Predicted values of LR number of flowers per stalk - HG site")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI
nbfpsLR.HG.new.data = expand.grid(time = rownames(coef(nbfpsLR.HG)))
nbfpsLR.HG.new.data$lwr = NA
nbfpsLR.HG.new.data$upr = NA
nbfpsLR.HG.new.data$pred = NA

# 95% confidence intervals: We cannot use the bootstrap here because the model is a glmmPQL. We thus use the easyPredCI function from Prof. Marc Girondot (https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

nbfpsLR.HG.new.data$lwr = apply(nbfpsLR.HG.new.data, 1, FUN = function(x) easyPredCI(nbfpsLR.HG, newdata = data.frame(time = as.numeric(x[1])))[1, 1])

nbfpsLR.HG.new.data$upr = apply(nbfpsLR.HG.new.data, 1, FUN = function(x) easyPredCI(nbfpsLR.HG, newdata = data.frame(time = as.numeric(x[1])))[1, 2])

nbfpsLR.HG.new.data$pred = apply(nbfpsLR.HG.new.data, 1, FUN = function(x) predict(nbfpsLR.HG, newdata = data.frame(time = as.numeric(x[1])), type = "response"))


# Plotting the predictions

plot.nbfpsLR.HG = ggplot(nbfpsLR.HG.new.data, aes(x = as.numeric(time), y = pred)) + 
 geom_line(size = 3) + 
 geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) + 
 xlab("Year") + 
 ylab("Number of flowers per stalk") + 
 ggtitle("Large reproductive number of flowers per stalks - High grazing") + 
 theme_bw() + 
 theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "HG_LRNbFlowPerStalks.png", 
 width = 4000, 
 height = 3000, 
 units = "px", 
 bg = "white", 
 res = 300, 
 type = "cairo")

plot.nbfpsLR.HG

dev.off()




###########################################################################
#
# 3. Creating tables with vital rates ----
#
###########################################################################

vr.per.year.TSF3plus = expand.grid(year = unique(data.grazing$time)[- length(unique(data.grazing$time))],
                                   grazing = c("low", "high"))


## 3.1. Seedling survival ----
# -----------------------

vr.per.year.TSF3plus$survSD = NA

vr.per.year.TSF3plus$survSD[which(vr.per.year.TSF3plus$grazing == "low")] = predict(survSD.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.LG$density), type = "response")

vr.per.year.TSF3plus$survSD[which(vr.per.year.TSF3plus$grazing == "high")] = predict(survSD.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.HG$density), type = "response")


## 3.2. Juvenile survival ----
# -----------------------

vr.per.year.TSF3plus$survJ = NA

vr.per.year.TSF3plus$survJ[which(vr.per.year.TSF3plus$grazing == "low")] = predict(survJ.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.LG$density), type = "response")

vr.per.year.TSF3plus$survJ[which(vr.per.year.TSF3plus$grazing == "high")] = predict(survJ.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")


## 3.3. Small reproductive survival ----
# ---------------------------------

vr.per.year.TSF3plus$survSR = NA

vr.per.year.TSF3plus$survSR[which(vr.per.year.TSF3plus$grazing == "low")] = predict(survSR.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")

vr.per.year.TSF3plus$survSR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(survSR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.HG$density), type = "response")


## 3.4. Large reproductive survival ----
# ---------------------------------

vr.per.year.TSF3plus$survLR = NA

vr.per.year.TSF3plus$survLR[which(vr.per.year.TSF3plus$grazing == "low")] = predict(survLR.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")

vr.per.year.TSF3plus$survLR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(survLR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")


## 3.5. J-LR transition ----
# ---------------------

vr.per.year.TSF3plus$transJ = NA

vr.per.year.TSF3plus$transJ[which(vr.per.year.TSF3plus$grazing == "low")] = predict(transJ.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")

vr.per.year.TSF3plus$transJ[which(vr.per.year.TSF3plus$grazing == "high")] = predict(transJ.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.HG$density), type = "response")


## 3.6. SR-LR transition ----
# ----------------------

vr.per.year.TSF3plus$transSR = NA

vr.per.year.TSF3plus$transSR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year != 2017)] = predict(transSR.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year[which(vr.per.year.TSF3plus$year != 2017)])), type = "response")
vr.per.year.TSF3plus$transSR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == 2017)] = 0

vr.per.year.TSF3plus$transSR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(transSR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")


## 3.7. LR-SR transition ----
# ----------------------

vr.per.year.TSF3plus$transLR = NA

vr.per.year.TSF3plus$transLR[which(vr.per.year.TSF3plus$grazing == "low")] = predict(transLR.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")

vr.per.year.TSF3plus$transLR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(transLR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.HG$density), type = "response")


## 3.8. SR probability of flowering ----
# ---------------------------------

vr.per.year.TSF3plus$floweringprobSR = NA

vr.per.year.TSF3plus$floweringprobSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year != 2016] = predict(floweringprobSR.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year[which(vr.per.year.TSF3plus$year != 2016)]), density = yearly.density.per.square.LG$density[which(yearly.density.per.square.LG$time != 2016)]), type = "response")
vr.per.year.TSF3plus$floweringprobSR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == 2016)] = mean(vr.per.year.TSF3plus$floweringprobSR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year != 2016)])

vr.per.year.TSF3plus$floweringprobSR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(floweringprobSR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)), type = "response")


## 3.9. LR probability of flowering ----
# ---------------------------------

vr.per.year.TSF3plus$floweringprobLR = NA

vr.per.year.TSF3plus$floweringprobLR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year != 2016)] = predict(floweringprobLR.LG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year[which(vr.per.year.TSF3plus$year != 2016)]), density = yearly.density.per.square.LG$density[which(yearly.density.per.square.LG$time != 2016)]), type = "response")
vr.per.year.TSF3plus$floweringprobLR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == 2016)] = mean(vr.per.year.TSF3plus$floweringprobLR[which(vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year != 2016)])

vr.per.year.TSF3plus$floweringprobLR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(floweringprobLR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.HG$density), type = "response")


## 3.10. SR number of flowering stalks ----
# ------------------------------------

vr.per.year.TSF3plus$nbfsSR = NA

vr.per.year.TSF3plus$nbfsSR[which(vr.per.year.TSF3plus$grazing == "low")] = nbfsSR.LG$fs[order(nbfsSR.LG$time)]

vr.per.year.TSF3plus$nbfsSR[which(vr.per.year.TSF3plus$grazing == "high")] = nbfsSR.HG$fs[which(nbfsSR.HG$time != 2019)]


## 3.11. LR number of flowering stalks ----
# ------------------------------------

vr.per.year.TSF3plus$nbfsLR = NA

vr.per.year.TSF3plus$nbfsLR[which(vr.per.year.TSF3plus$grazing == "low")] = nbfsLR.LG.pred$pred[order(nbfsLR.LG.pred$time)][- length(nbfsLR.LG.pred$pred)]

vr.per.year.TSF3plus$nbfsLR[which(vr.per.year.TSF3plus$grazing == "high")] = nbfsLR.HG$fs[which(nbfsLR.HG$time != 2019)]


## 3.12. SR number of flowers per stalk ----
# -------------------------------------

vr.per.year.TSF3plus$nbfpsSR = NA

vr.per.year.TSF3plus$nbfpsSR[which(vr.per.year.TSF3plus$grazing == "low")] = nbfpsSR.LG$fps[order(nbfpsSR.LG$time)]

vr.per.year.TSF3plus$nbfpsSR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(nbfpsSR.HG.quasiP, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year), density = yearly.density.per.square.HG$density))


## 3.13. LR number of flowers per stalk ----
# -------------------------------------

vr.per.year.TSF3plus$nbfpsLR = NA

vr.per.year.TSF3plus$nbfpsLR[which(vr.per.year.TSF3plus$grazing == "low")] = nbfpsLR.LG.pred$pred[order(nbfpsLR.LG.pred$time)][- length(nbfpsLR.LG.pred$time)]

vr.per.year.TSF3plus$nbfpsLR[which(vr.per.year.TSF3plus$grazing == "high")] = predict(nbfpsLR.HG, newdata = data.frame(time = unique(vr.per.year.TSF3plus$year)))




###########################################################################
#
# 4. Saving models and tables ----
#
###########################################################################

## 4.1. Saving vital rates table ----
# ------------------------------

write.csv(vr.per.year.TSF3plus, "DPVR_TSF3plus.csv", row.names = F)


## 4.2. Saving models ----
# -------------------

## 4.2.1. Seedling survival ----
# -----------------------

save(survSD.LG, file = "GLMM_survSD_LG.RData")

save(survSD.HG, file = "GLMM_survSD_HG.RData")


## 4.2.2. Juvenile survival ----
# -----------------------

save(survJ.LG, file = "GLMM_survJ_LG.RData")

save(survJ.HG, file = "GLMM_survJ_HG.RData")


## 4.2.3. Small reproductive survival ----
# ---------------------------------

save(survSR.LG, file = "GLMM_survSR_LG.RData")

save(survSR.HG, file = "GLMM_survSR_HG.RData")


## 4.2.4. Large reproductive survival ----
# ---------------------------------

save(survLR.LG, file = "GLMM_survLR_LG.RData")

save(survLR.HG, file = "GLMM_survLR_HG.RData")


## 4.2.5. J-LR transition ----
# ---------------------

save(transJ.LG, file = "GLMM_transJ_LG.RData")

save(transJ.HG, file = "GLMM_transJ_HG.RData")


## 4.2.6. SR-LR transition ----
# ----------------------

save(transSR.LG, file = "GLMM_transSR_LG.RData")

save(transSR.HG, file = "GLMM_transSR_LG.RData")


## 4.2.7. LR-SR transition ----
# ----------------------

save(transLR.LG, file = "GLMM_transLR_LG.RData")

save(transLR.HG, file = "GLMM_transLR_HG.RData")


## 4.2.8. SR probability of flowering ----
# ---------------------------------

save(floweringprobSR.LG, file = "GLMM_floweringprobSR_LG.RData")

save(floweringprobSR.HG, file = "GLMM_floweringprobSR_HG.RData")


## 4.2.9. LR probability of flowering ----
# ---------------------------------

save(floweringprobLR.LG, file = "GLMM_floweringprobLR_LG.RData")

save(floweringprobLR.HG, file = "GLMM_floweringprobLR_HG.RData")


## 4.2.10. LR number of flowering stalks ----
# --------------------------------------

save(nbfsLR.LG, file = "GLMM_nbfsLR_LG.RData")


## 4.2.11. SR number of flowers per stalk ----
# -------------------------------------

save(nbfpsSR.HG, file = "GLMM_nbfpsSR_HG.RData")


## 4.2.12. LR number of flowers per stalk ----
# -------------------------------------

save(nbfpsLR.LG, file = "GLMM_nbfpsLR_LG.RData")

save(nbfpsLR.HG, file = "GLMM_nbfpsLR_HG.RData")

