##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al. 2022, Ecology).
# 
# Eva Conquet, Arpat Ozgul, Daniel T. Blumstein, Kenneth B. Armitage,
# Madan K. Oli, Julien G. A. Martin, Tim H. Clutton-Brock, and Maria Paniw
# 
# This script uses the capture-recapture data of a marmot population in the upper East River Valley nearby 
# Gothic, Colorado, United States, collected between 1976 and 2016. 
# The aim of this script is to model the survival, transitions, and recruitment of four 
# life-history stages: juvenile, yearling, non-reproductive, and reproductive adults.
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
  library(lme4)
  library(MuMIn)
  library(bbmle)
  library(boot)
  library(ggplot2)
}

load.librairies()


## 1.3. Loading functions ----
# -----------------------

# Function to test for overdispersion in GLMMs.
# Taken from material from Prof. Ben Bolker (accessible at https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

overdisp_fun <- function(model){
  rdf <- df.residual(model)
  rp <- residuals(model,type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}


# Function to compute 95 % CIs for glmmPQL models
# Taken from Prof. Marc Girondot (accessible at https://biostatsr.blogspot.com/2016/02/predict-for-glm-and-glmm.html). 

easyPredCI <- function(model, newdata = NULL, alpha = 0.05){
  # Marc Girondot - 2016-01-09
  if (is.null(newdata)){
    
    if (any(class(model) == "glmerMod")) newdata <- model@frame
    if (any(class(model) == "glmmPQL") | any(class(model) == "glm")) newdata <- model$data
    if (any(class(model) == "glmmadmb")) newdata <- model$frame
  }
  
  ## baseline prediction, on the linear predictor scale:
  pred0 <- predict(model, re.form = NA, newdata = newdata)
  ## fixed-effects model matrix for new data
  if (any(class(model) == "glmmadmb")){
    X <- model.matrix(delete.response(model$terms), newdata)
  } 
  else{
    X <- model.matrix(formula(model, fixed.only = TRUE)[-2],
                      newdata)
  }
  
  if (any(class(model) == "glm")){
    # Marc Girondot - 2016-01-09
    # Note that beta is not used
    beta <- model$coefficients
  } 
  else{
    beta <- fixef(model) ## fixed-effects coefficients
  }
  
  V <- vcov(model)     ## variance-covariance matrix of beta
  
  # Marc Girondot - 2016-01-09
  if (any(!(colnames(V) %in% colnames(X)))) {
    dfi <- matrix(data = rep(0, dim(X)[1] * sum(!(colnames(V) %in% colnames(X)))), nrow = dim(X)[1])
    colnames(dfi) <- colnames(V)[!(colnames(V) %in% colnames(X))]
    X <- cbind(X, dfi)
  }
  
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  
  ## inverse-link function
  # Marc Girondot - 2016-01-09
  if (any(class(model) == "glmmPQL") | any(class(model) == "glm")) linkinv <- model$family$linkinv
  if (any(class(model) == "glmerMod")) linkinv <- model@resp$family$linkinv
  if (any(class(model) == "glmmadmb")) linkinv <- model$ilinkfun
  
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(lwr = pred0 - crit * pred.se,
                upr = pred0 + crit * pred.se))
}


## 1.4. Loading and preparing data ----
# --------------------------------

data.marmots = read.csv("MarmotsData.csv")
head(data.marmots)


# Defining year as a factor

data.marmots$year = factor(data.marmots$year)




###########################################################################
#
# 2. Fitting Generalized Linear Mixed Models (GLMMs) ----
# to estimate the vital rates
#
###########################################################################

## 2.1. Survival (surv) ----
# ---------------------

# 2.1.1. Juveniles (J) ----
# --------------------

# Only for the winter-summer transition, as juveniles don't "survive" in the summer-winter transition, there are only newborns from the reproductive adults in this transition.

plot(data.marmots$surv[data.marmots$stage_surv == "J"])


# Finding best random effect 

surv.J1 = glmer(surv ~ 1 + (1|year), data.marmots[data.marmots$stage_surv == "J" & data.marmots$season == "winter", ], family = binomial)


# Testing the part of variance explained by the RE structure

r.squaredGLMM(surv.J1, null = glmer(surv ~ 1 + (1|year), data.marmots[data.marmots$stage_surv == "J" & data.marmots$season == "winter", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.

summary(surv.J1)
survJ = surv.J1


# Plotting the observed data vs the predictions of the model

survJ.obs = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "J", ], FUN = mean)

survJ.pred = data.frame(levels(data.marmots$year), predict(survJ, newdata = expand.grid(year = levels(data.marmots$year)), type = "response"))
colnames(survJ.pred) = c("year", "pred")


# Obs vs. Pred

plot(survJ.obs$surv, survJ.pred$pred[which(survJ.pred$year %in% survJ.obs$year)], pch = 16, ylim = c(0, 1), xlab = "Observed survival probability",ylab="Predicted survival probability",main="Observed vs. Predicted values of juveniles (J) survival")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survJ.new.data = expand.grid(year = levels(data.marmots$year))


# 95% confidence intervals

bootstrap.survJ <- bootMer(survJ, FUN = function(x) predict(x, survJ.new.data),
                           nsim = 1000)

survJ.new.data$lwr = inv.logit(apply(bootstrap.survJ$t, 2, quantile, 0.025, na.rm = T))
survJ.new.data$upr = inv.logit(apply(bootstrap.survJ$t, 2, quantile, 0.975, na.rm = T))
survJ.new.data$pred = predict(survJ, newdata = survJ.new.data, type = "response")


# Plotting the predictions

plot.survJ = ggplot(survJ.new.data, aes(x = as.numeric(as.character(year)), y = pred)) + 
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Survival probability") +
  ggtitle("Juvenile survival") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), axis.text.x = element_text(size=40,colour="black", margin = margin(t=10,r=0,b=0,l=0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15,l = 0), face = "bold"), 
        legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), 
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



# 2.1.2. Yearlings (Y) ---- 
# --------------------

plot(data.marmots$surv[data.marmots$stage_surv == "Y"])


# Finding best random effect

surv.Y1 = glmer(surv ~ season + (1|year), data.marmots[data.marmots$stage_surv == "Y", ], family = binomial)
surv.Y2 = glmer(surv ~ season + (1 + season|year), data.marmots[data.marmots$stage_surv == "Y", ], family = binomial)


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.Y1, null = glmer(surv ~ 1 + (1|year), data.marmots[data.marmots$stage_surv == "Y", ], family = binomial)) 
r.squaredGLMM(surv.Y2, null = glmer(surv ~ 1 + (1 + season|year), data.marmots[data.marmots$stage_surv == "Y", ], family = binomial)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

surv.Y3 = glmer(surv ~ 1 + (1 + season|year), data.marmots[data.marmots$stage_surv == "Y", ], family = binomial)

AICctab(surv.Y2, surv.Y3) 

summary(surv.Y3)
survY = surv.Y3


# Plotting the observed data vs the predictions of the model

survY.obs.winter = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "Y" & data.marmots$season == "winter", ], FUN = mean)
survY.obs.summer = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "Y" & data.marmots$season == "summer", ], FUN = mean)


# Winter

survY.pred.winter = data.frame(year = levels(data.marmots$year), 
                               season = "winter", 
                               pred = predict(survY, newdata = data.frame(year = levels(data.marmots$year), season = "winter"), type = "response"))


# Obs vs. Pred

plot(survY.obs.winter$surv, survY.pred.winter$pred[which(survY.pred.winter$year %in% survY.obs.winter$year)], 
     pch = 16, ylim = c(0, 0.9), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of yearlings (Y) survival - Winter")
abline(a = 0, b = 1, col = "red")


# Summer

survY.pred.summer = data.frame(year = levels(data.marmots$year), 
                               season = "summer", 
                               pred = predict(survY, newdata = data.frame(year = levels(data.marmots$year), season = "summer"), type = "response"))


# Obs vs. Pred

plot(survY.obs.summer$surv, survY.pred.summer$pred[which(survY.pred.summer$year %in% survY.obs.summer$year)], 
     pch = 16, xlim = c(0, 1), ylim = c(0, 0.9), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of yearlings (Y) survival - Summer")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survY.new.data = expand.grid(year = levels(data.marmots$year),
                             season = unique(data.marmots$season))


# 95% confidence intervals

bootstrap.survY <- bootMer(survY, FUN = function(x) predict(x, survY.new.data),
                           nsim = 1000)
survY.new.data$lwr = inv.logit(apply(bootstrap.survY$t, 2, quantile, 0.025, na.rm = T))
survY.new.data$upr = inv.logit(apply(bootstrap.survY$t, 2, quantile, 0.975, na.rm = T))
survY.new.data$pred = predict(survY, newdata = survY.new.data, type = "response")

survY.new.data.mean = data.frame(season = c("winter", "summer"), 
                                 pred = c(mean(survY.new.data$pred[which(survY.new.data$season == "winter")]), mean(survY.new.data$pred[which(survY.new.data$season == "summer")])), 
                                 upr = c(mean(survY.new.data$upr[which(survY.new.data$season == "winter")]), mean(survY.new.data$upr[which(survY.new.data$season == "summer")])), 
                                 lwr = c(mean(survY.new.data$lwr[which(survY.new.data$season == "winter")]), mean(survY.new.data$lwr[which(survY.new.data$season == "summer")])))


# Plotting the predictions

plot.survY = ggplot(survY.new.data.mean, aes(x = season, y = pred)) + 
  geom_point(size = 10, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, position = position_dodge(width = 0.4)) +
  xlab("Year") +
  ylab("Survival probability") +
  ggtitle("Yearling survival") +
  scale_x_discrete(name = "Season",
                   labels = c("winter" = "Winter", "summer" = "Summer")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0 , l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "YearlingSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survY

dev.off()



# 2.1.3. Non Reproductive Adults (NRA) ----
# ------------------------------------

plot(data.marmots$surv[data.marmots$stage_surv == "NRA"])


# Finding best random effect

surv.NR1 = glmer(surv ~ season + (1|year), data.marmots[data.marmots$stage_surv == "NRA", ], family = binomial)
surv.NR2 = glmer(surv ~ season + (1 + season|year), data.marmots[data.marmots$stage_surv == "NRA", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.NR1, null = glmer(surv ~ 1 + (1|year), data.marmots[data.marmots$stage_surv == "NRA", ], family = binomial)) 
r.squaredGLMM(surv.NR2, null = glmer(surv ~ 1 + (1+season|year), data.marmots[data.marmots$stage_surv == "NRA", ], family = binomial)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

surv.NR3 = glmer(surv ~ 1 + (1 + season|year), data.marmots[data.marmots$stage_surv == "NRA", ], family = binomial)

AICctab(surv.NR2, surv.NR3, base = T) 

summary(surv.NR3)
survNR = surv.NR3


# Plotting the observed data vs the predictions of the model

survNR.obs.winter = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "NRA" & data.marmots$season == "winter", ], FUN = mean)
survNR.obs.summer = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "NRA" & data.marmots$season == "summer", ], FUN = mean)


# Winter

survNR.pred.winter = data.frame(year = levels(data.marmots$year), 
                                 season = "winter", 
                                 pred = predict(survNR, newdata = data.frame(year = levels(data.marmots$year), season = "winter"), type = "response"))


# Obs vs. Pred

plot(survNR.obs.winter$surv, survNR.pred.winter$pred, 
     pch = 16, ylim = c(0, 0.9), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of non reproductive adults (NRA) survival - Winter")
abline(a = 0, b = 1, col = "red")


# Summer

survNR.pred.summer = data.frame(year = levels(data.marmots$year), 
                                 season = "summer", 
                                 pred = predict(survNR, newdata = data.frame(year = levels(data.marmots$year), season = "summer"), type = "response"))


# Obs vs. Pred

plot(survNR.obs.summer$surv, survNR.pred.summer$pred[which(survNR.pred.summer$year %in% survNR.obs.summer$year)], 
     pch = 16, xlim = c(0, 1), ylim = c(0, 0.9), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of non reproductive adults (NRA) survival - Summer")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survNR.new.data = expand.grid(year = levels(data.marmots$year),
                               season = unique(data.marmots$season))


# 95% confidence intervals

bootstrap.survNR <- bootMer(survNR, FUN = function(x) predict(x, survNR.new.data,re.form = NA),
                             nsim = 1000)
survNR.new.data$lwr = inv.logit(apply(bootstrap.survNR$t, 2, quantile, 0.025, na.rm = T))
survNR.new.data$upr = inv.logit(apply(bootstrap.survNR$t, 2, quantile, 0.975, na.rm = T))
survNR.new.data$pred = predict(survNR, newdata = survNR.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survNR = ggplot(survNR.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 10, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, position = position_dodge(width = 0.4)) +
  xlab("Season") +
  ylab("Survival probability") +
  ggtitle("Non reproductive adult survival") +
  scale_x_discrete(name = "Season",
                   labels = c("winter" = "Winter", "summer" = "Summer")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"),
        legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"))

png(filename = "NonReproAdultSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survNR

dev.off()



# 2.1.4. Reproductive Adults (RA) ----
# -------------------------------

plot(data.marmots$surv[data.marmots$stage_surv=="RA"])


# Finding best random effect

surv.R1 = glmer(surv ~ season + (1|year), data = data.marmots[data.marmots$stage_surv == "RA", ], family = binomial)
surv.R2 = glmer(surv ~ season + (1 + season|year), data = data.marmots[data.marmots$stage_surv == "RA", ], family = binomial, control = glmerControl(optimizer = "bobyqa", calc.derivs = F))


# Testing the part of variance explained by each model (i.e. by each RE structure)

r.squaredGLMM(surv.R1, null = glmer(surv ~ 1 + (1|year), data.marmots[data.marmots$stage_surv == "RA", ], family = binomial)) 
r.squaredGLMM(surv.R2, null = glmer(surv ~ 1 + (1 + season|year), data.marmots[data.marmots$stage_surv == "RA", ], family = binomial)) # The variance explained by this model is slightly higher. We thus keep the random effect on the intercept and the slope.


# Finding best fixed effect

surv.R3 = glmer(surv ~ 1 + (1 + season|year), data.marmots[data.marmots$stage_surv == "RA", ], family = binomial)

AICctab(surv.R2, surv.R3)

summary(surv.R2)
survR = surv.R2


# Plotting the observed data vs the predictions of the model

survR.obs.winter = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "RA" & data.marmots$season == "winter", ], FUN = mean)
survR.obs.summer = aggregate(surv ~ year, data = data.marmots[data.marmots$stage_surv == "RA" & data.marmots$season == "summer", ], FUN = mean)


# Winter

survR.pred.winter = data.frame(year = levels(data.marmots$year), 
                                season = "winter", 
                                pred = predict(survR, newdata = data.frame(year = levels(data.marmots$year), season = "winter"), type = "response"))


# Obs vs. Pred

plot(survR.obs.winter$surv, survR.pred.winter$pred, 
     pch = 16, xlim = c(0, 1), ylim = c(0, 0.9), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of reproductive adults (RA) survival - Winter")
abline(a = 0, b = 1, col = "red")


# Summer

survR.pred.summer = data.frame(year = levels(data.marmots$year), 
                                season = "summer", 
                                pred = predict(survR, newdata = data.frame(year = levels(data.marmots$year), season = "summer"), type = "response"))


# Obs vs. Pred

plot(survR.obs.summer$surv, survR.pred.summer$pred,
     pch = 16, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "Observed survival probability", ylab = "Predicted survival probability", 
     main = "Observed vs. Predicted values of reproductive adults (RA) survival - Summer")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

survR.new.data = expand.grid(year = levels(data.marmots$year),
                              season = unique(data.marmots$season))


# 95% confidence intervals

bootstrap.survR <- bootMer(survR, FUN = function(x) predict(x, survR.new.data, re.form = NA),
                            nsim = 1000)
survR.new.data$lwr = inv.logit(apply(bootstrap.survR$t, 2, quantile, 0.025, na.rm = T))
survR.new.data$upr = inv.logit(apply(bootstrap.survR$t, 2, quantile, 0.975, na.rm = T))
survR.new.data$pred = predict(survR, newdata = survR.new.data, type = "response", re.form = NA)


# Plotting the predictions

plot.survR = ggplot(survR.new.data, aes(x = season, y = pred)) + 
  #facet_wrap(~year)+
  geom_point(size = 10, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, position = position_dodge(width = 0.4)) +
  xlab("Season") +
  ylab("Reproductive adult\nsurvival probability") +
  #ggtitle("Reproductive adult survival") +
  scale_x_discrete(name = "Season",
                   labels = c("winter" = "Winter", "summer" = "Summer")) +
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

png(filename = "ReproAdultSurvival.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.survR

dev.off()




## 2.2. Recruitment of reproductive adults (recruits, summer only)
# ----------------------------------------------------------------

hist(data.marmots$recruits[data.marmots$season == "summer"])


# Finding best random effect

recruit1 = glmer(recruits ~ 1 + (1|year), data = data.marmots[data.marmots$stage_surv == "RA" & data.marmots$season == "summer", ], family = poisson)


# Testing the part of variance explained by each the RE structure

r.squaredGLMM(recruit1, null = glmer(recruits ~ 1 + (1|year), data.marmots[data.marmots$stage_surv == "RA" & data.marmots$season == "summer", ], family = poisson)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.

summary(recruit1)
recruit = recruit1


# Checking for over/underdispersion and fitting a quasi-Poisson model

overdisp_fun(recruit)

recruit.quasiP = MASS::glmmPQL(recruits ~ 1, random =  ~ 1|year, family = poisson, data = data.marmots[data.marmots$stage_surv == "RA" & data.marmots$season == "summer", ])

summary(recruit.quasiP)

recruit = recruit.quasiP


# Plotting the observed data vs the predictions of the model

recruit.obs = aggregate(recruits ~ year, data = data.marmots[data.marmots$stage_surv == "RA", ], FUN = mean)

recruit.pred = data.frame(year = levels(data.marmots$year), 
                          pred = predict(recruit, newdata = data.frame(year = levels(data.marmots$year)), type = "response"))


# Obs vs. Pred

plot(recruit.obs$recruits, recruit.pred$pred, 
     ylim = c(0, 7), xlim = c(0, 7), 
     pch = 16, xlab = "Observed recruitment", ylab = "Predicted recruitment", 
     main = "Observed vs. Predicted values of recruitment")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

recruit.new.data = data.frame(year = levels(data.marmots$year))

recruit.new.data$lwr = easyPredCI(recruit, newdata = recruit.new.data)[, 1]
recruit.new.data$upr = easyPredCI(recruit, newdata = recruit.new.data)[, 2]
recruit.new.data$pred = predict(recruit, newdata = recruit.new.data, type = "response")


# Plotting the predictions

plot.recruit = ggplot(recruit.new.data, aes(x = seq(1976, 2016), y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Number of pups per female") +
  ggtitle("Recruitment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0 )), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"),
        legend.text = element_text(size = 40), 
        legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", 
        legend.key.size = unit(6,"lines"))


png(filename = "Recruitment.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.recruit

dev.off()




## 2.3. Transition to reproductive adults (transRA, winter only) ----
# --------------------------------------------------------------

plot(data.marmots$transRA[data.marmots$season=="winter"])


# 2.3.1. Yearlings ----
# ----------------

# Finding best random effect

transYR1 = glmer(transRA ~ 1 + (1|year), data = data.marmots[data.marmots$season == "winter" & data.marmots$stage_surv == "Y", ], family = binomial)


# Testing the part of variance explained by each the RE structure

r.squaredGLMM(transYR1, null = glmer(transRA ~ 1 + (1|year), data = data.marmots[data.marmots$season == "winter" & data.marmots$stage_surv == "Y", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.

summary(transYR1)
transYR = transYR1


# Plotting the observed data vs the predictions of the model

transYR.obs = aggregate(transRA ~ year, data = data.marmots[data.marmots$stage_surv == "Y", ], FUN = mean)


# Creating a dataframe with predictions of the model to have an estimate for every year from 1976 to 2015

coef(transYR)$year # There is no estimate for 1983
table(data.marmots$transRA[which(data.marmots$year == 1983)], data.marmots$stage_surv[which(data.marmots$year == 1983)]) # There are no transitions for the yearlings in 1983 so we set it to 0.
transYR.pred = c()

for(year in 1976:2015){
  
  if(year!=1983){
    
    transYR.pred = c(transYR.pred, predict(transYR, newdata = data.frame(year = year), type = "response"))
  }
  
  else{
    
    transYR.pred = c(transYR.pred, 0)
  }
}

transYR.pred = data.frame(year = seq(1976, 2015), pred = transYR.pred)


# Obs vs. Pred

plot(transYR.obs$transRA, transYR.pred$pred[which(transYR.pred$year %in% transYR.obs$year)], 
     pch = 16, xlab = "Observed transition probability", ylab = "Predicted transition probability", 
     main = "Observed vs. Predicted values of Y-RA transition probability")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transYR.new.data = expand.grid(year = c(seq(1976, 1982), seq(1984, 2015)))


# 95% confidence intervals

bootstrap.recruit <- bootMer(transYR, FUN = function(x) predict(x, transYR.new.data),
                             nsim = 1000)
transYR.new.data$lwr = inv.logit(apply(bootstrap.recruit$t, 2, quantile, 0.025, na.rm = T))
transYR.new.data$upr = inv.logit(apply(bootstrap.recruit$t, 2, quantile, 0.975, na.rm = T))
transYR.new.data$pred = predict(transYR,newdata = transYR.new.data, type = "response")


# Plotting the predictions

plot.transYR = ggplot(transYR.new.data, aes(x = year, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Transition probability") +
  ggtitle("Yearling - Reproductive adult transition") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 45, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 45, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 40, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = 40, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size = 50, hjust = 0.5, margin = margin(t = 0, r = 0, b = 15, l = 0), face = "bold"), 
        legend.text = element_text(size = 40), legend.title = element_text(face = "bold", size = 40), 
        legend.position = "right", legend.key.size = unit(6, "lines"))

png(filename = "Y-RATransition.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.transYR

dev.off()


# 2.3.2. Non-reproductive adults ----
# ------------------------------

# Finding best random effect

transNRR1 = glmer(transRA ~ 1 + (1|year), data = data.marmots[data.marmots$season == "winter" & data.marmots$stage_surv == "NRA", ], family = binomial)

# Testing the part of variance explained by each the RE structure

r.squaredGLMM(transNRR1, null = glmer(transRA ~ 1 + (1|year), data = data.marmots[data.marmots$season == "winter" & data.marmots$stage_surv == "NRA", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.

summary(transNRR1)
transNRR = transNRR1


# Plotting the observed data vs the predictions of the model

transNRR.obs = aggregate(transRA ~ year, data = data.marmots[data.marmots$stage_surv == "NRA", ], FUN = mean)


# Creating a dataframe with predictions of the model to have an estimate for every year from 1976 to 2015

coef(transNRR)$year # There is no estimate for 1977, 1981, 2013, and 2014
table(data.marmots$transRA[which(data.marmots$year %in% c(1977, 1981, 2013, 2014))], data.marmots$stage_surv[which(data.marmots$year %in% c(1977, 1981, 2013, 2014))], data.marmots$year[which(data.marmots$year %in% c(1977, 1981, 2013, 2014))]) # There are no transitions for the non-reproductive adults in 1977, 1981, 2013, and 2014 so we set it to 0.
transNRR.pred = c()

for(year in 1976:2015){

  if(!(year %in% c(1977, 1981, 2013, 2014))){
    
    transNRR.pred = c(transNRR.pred, predict(transNRR, newdata = data.frame(year = year), type = "response"))
  }
  
  else{
    transNRR.pred = c(transNRR.pred, 0)
  }
}

transNRR.pred = data.frame(year = seq(1976, 2015), pred = transNRR.pred)


# Obs vs. Pred

plot(transNRR.obs$transRA, transNRR.pred$pred[which(transNRR.pred$year %in% transNRR.obs$year)], ylim = c(0, 1), 
     pch = 16, xlab = "Observed transition probability", ylab = "Predicted transition probability", 
     main = "Observed vs. Predicted values of NRA-RA transition")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transNRR.new.data = expand.grid(year = c(1976, seq(1978, 1980), seq(1982, 2012), 2015))


# 95% confidence intervals

bootstrap.recruit <- bootMer(transNRR, FUN = function(x) predict(x, transNRR.new.data),
                             nsim = 1000)
transNRR.new.data$lwr = inv.logit(apply(bootstrap.recruit$t, 2, quantile, 0.025, na.rm = T))
transNRR.new.data$upr = inv.logit(apply(bootstrap.recruit$t, 2, quantile, 0.975, na.rm = T))
transNRR.new.data$pred = predict(transNRR, newdata = transNRR.new.data, type = "response")


# Plotting the predictions

plot.transNRR = ggplot(transNRR.new.data, aes(x = year, y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Transition probability") +
  ggtitle("Non reproductive - Reproductive adult transition") +
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

png(filename = "NRA-RATransition.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.transNRR

dev.off()



# 2.3.3. Reproductive adults ----
# --------------------------

# Finding best random effect

transRR1 = glmer(transRA ~ 1 + (1|year), data = data.marmots[data.marmots$season == "winter" & data.marmots$stage_surv == "RA", ], family = binomial)


# Testing the part of variance explained by each the RE structure

r.squaredGLMM(transRR1, null = glmer(transRA ~ 1 + (1|year), data = data.marmots[data.marmots$season == "winter" & data.marmots$stage_surv == "RA", ], family = binomial)) # The RE does explain part of the variance. We thus keep the random effect on the intercept only.

summary(transRR1)
transRR = transRR1


# Plotting the observed data vs the predictions of the model

transRR.obs = aggregate(transRA ~ year, data = data.marmots[data.marmots$stage_surv == "RA", ], FUN = mean)

transRR.pred = data.frame(year = levels(data.marmots$year)[1:40], 
                             pred = predict(transRR, newdata = data.frame(year = levels(data.marmots$year)[1:40], season = "winter"),type = "response"))


# Obs vs. Pred

plot(transRR.obs$transRA, transRR.pred$pred, ylim = c(0, 1), 
     pch = 16, xlab = "Observed transition probability", ylab = "Predicted transition probability", 
     main = "Observed vs. Predicted values of RA-RA transition")
abline(a = 0, b = 1, col = "red")


# Computing the predictions of the model and the 95 % CI

transRR.new.data = expand.grid(year = levels(data.marmots$year)[1:40])


# 95% confidence intervals

bootstrap.recruit <- bootMer(transRR, FUN = function(x) predict(x, transRR.new.data),
                             nsim = 1000)
transRR.new.data$lwr = inv.logit(apply(bootstrap.recruit$t, 2, quantile, 0.025, na.rm = T))
transRR.new.data$upr = inv.logit(apply(bootstrap.recruit$t, 2, quantile, 0.975, na.rm = T))
transRR.new.data$pred = predict(transRR, newdata = transRR.new.data, type = "response")


# Plotting the predictions

plot.transRR = ggplot(transRR.new.data, aes(x = seq(1976, 2015), y = pred)) + 
  #facet_wrap(~year)+
  geom_line(size = 3) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), linetype = 0, alpha = 0.2) +
  xlab("Year") +
  ylab("Stasis probability") +
  ggtitle("Reproductive adult stasis") +
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

png(filename = "RA-RATransition.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

plot.transRR

dev.off()




###########################################################################
#
# 3. Saving the vital rates ----
#
###########################################################################

vr.per.year = data.frame(year = unique(data.marmots$year))

# Survival of juveniles

surv.J.per.year = predict(survJ, newdata = data.frame(year = vr.per.year$year), type = "response")
vr.per.year$survJ = surv.J.per.year

# Survival of yearlings

surv.Y.per.year.winter = predict(survY, newdata = data.frame(year = vr.per.year$year, season = "winter"), type = "response")
surv.Y.per.year.summer = predict(survY, newdata = data.frame(year = vr.per.year$year, season = "summer"), type = "response")
vr.per.year$survY_winter = surv.Y.per.year.winter
vr.per.year$survY_summer = surv.Y.per.year.summer

# Survival of non reproductive adults

surv.NR.per.year.winter = predict(survNR, newdata = data.frame(year = vr.per.year$year, season = "winter"), type = "response")
surv.NR.per.year.summer = predict(survNR, newdata = data.frame(year = vr.per.year$year, season = "summer"), type = "response")
vr.per.year$survNR_winter = surv.NR.per.year.winter
vr.per.year$survNR_summer = surv.NR.per.year.summer

# Survival of reproductive adults

surv.R.per.year.winter = predict(survR, newdata = data.frame(year = vr.per.year$year, season = "winter"), type = "response")
surv.R.per.year.summer = predict(survR, newdata = data.frame(year = vr.per.year$year, season = "summer"), type = "response")
vr.per.year$survR_winter = surv.R.per.year.winter
vr.per.year$survR_summer = surv.R.per.year.summer

# Transition Y -> R

trans.Y.R.per.year = data.frame(year = rownames(coef(transYR)$year), 
                                transYR = predict(transYR, newdata = data.frame(year = rownames(coef(transYR)$year)), type = "response"))
trans.Y.R.per.year = rbind(trans.Y.R.per.year, data.frame(year = vr.per.year$year[which(!(vr.per.year$year %in% trans.Y.R.per.year$year))], 
                                                          transYR = rep(NA, length(which(!(vr.per.year$year %in% trans.Y.R.per.year$year))))))
trans.Y.R.per.year = trans.Y.R.per.year[order(trans.Y.R.per.year$year), ]
vr.per.year$transYR = trans.Y.R.per.year$transYR

# Transition NR -> R

trans.NR.R.per.year = data.frame(year = rownames(coef(transNRR)$year), 
                                 transNRR = predict(transNRR, newdata = data.frame(year = rownames(coef(transNRR)$year)), type = "response"))
trans.NR.R.per.year = rbind(trans.NR.R.per.year, data.frame(year = vr.per.year$year[which(!(vr.per.year$year %in% trans.NR.R.per.year$year))], 
                                                            transNRR = rep(NA, length(which(!(vr.per.year$year %in% trans.NR.R.per.year$year))))))
trans.NR.R.per.year = trans.NR.R.per.year[order(trans.NR.R.per.year$year), ]
vr.per.year$transNRR = trans.NR.R.per.year$transNRR

# Transition R -> R

trans.R.R.per.year = data.frame(year = rownames(coef(transRR)$year), 
                                transRR = predict(transRR, newdata = data.frame(year = rownames(coef(transRR)$year)), type = "response"))
trans.R.R.per.year = rbind(trans.R.R.per.year, data.frame(year = vr.per.year$year[which(!(vr.per.year$year %in% trans.R.R.per.year$year))], transRR = rep(NA, length(which(!(vr.per.year$year %in% trans.R.R.per.year$year))))))
trans.R.R.per.year = trans.R.R.per.year[order(trans.R.R.per.year$year), ]
vr.per.year$transRR = trans.R.R.per.year$transRR

# Recruitment of R individuals

recruit.per.year = predict(recruit, newdata = data.frame(year = vr.per.year$year), type = "response")
vr.per.year$recruit = recruit.per.year

# Turning NAs in transitions into 0s

vr.per.year[is.na(vr.per.year)] = 0




###########################################################################
#
# 4. Saving files and RData ----
#
###########################################################################

write.csv(vr.per.year, "Marmots_VitalRatesEstimates.csv", row.names = F)
save(survJ, file = "GLMM_survJ.RData")
save(survY, file = "GLMM_survY.RData")
save(survNR, file = "GLMM_survNR.RData")
save(survR, file = "GLMM_survR.RData")
save(recruit, file = "GLMM_recruit.RData")
save(transYR, file = "GLMM_transYR.RData")
save(transNRR, file = "GLMM_transNRR.RData")
save(transRR, file = "GLMM_transRR.RData")
