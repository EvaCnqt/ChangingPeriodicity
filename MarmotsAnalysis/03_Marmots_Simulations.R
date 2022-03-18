##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
#
# This script contains the function needed to project the marmot population dynamics under scenarios
# of higher or lower seasonality in vital rates. 
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
}

load.librairies()


## 1.3. Loading data ----
# ------------------

source("02_Marmots_MPMs.R")

marmots.vr = read.csv("Marmots_VitalRatesEstimates.csv")

load("GLMM_survY.RData")
load("GLMM_survNR.RData")
load("GLMM_survR.RData")




###########################################################################
#
# 2. Projection function ----
#
###########################################################################

stoch.sim.marmot <- function(nsimul = 500,
                             nyears = 100,
                             n0,
                             vr.model = NULL,
                             seas.threshold = 0.5,
                             seasonality = "control"){
  
  # (1) Control case = Simulating what happens in the normal conditions, with periodic vital rates and stochastic variation across the years - Random sampling of the year among all available ones (1976 to 2015).
  # (2) High seasonality case = Random sampling only among years with a significant difference between seasons in season-dependent vital rates (adults survival).
  # (3) Low seasonality case = Random sampling only among years with a non-significant difference between seasons in season-dependent vital rates (adults survival).
  
  
  # Error if seasonality is either "high" or "low" but no model has been given as a reference to find significant years
  
  if(is.null(vr.model) & any(seasonality == c("high", "low"))){
    
    stop("Please enter the reference vital rate model to determine significant/nonsignificant years for seasonality!")
  }
  
  
  # Storage arrays and vectors for log lambdas and population extinction
  
  log.lambda = array(NA, c(nsimul, nyears))
  ext = seq(0, 0, length.out = nsimul)
  
  
  # If specific seasonality significance, find related years: take significative slope at the response scale
  
  years = marmots.vr$year[1:40]
  
  # Compute absolute changes
  
  if(seasonality != "control"){
    
    abs.changes = abs(inv.logit(coef(vr.model)$year$`(Intercept)` + coef(vr.model)$year$seasonwinter) - inv.logit(coef(vr.model)$year$`(Intercept)`))[1:40]
    
  }
  
  # Get high- and low-seasonality years
  
  if(seasonality == "high"){
    
    high.seas.years = abs.changes >= quantile(abs.changes, probs = seas.threshold)
    high.seas.years = years[high.seas.years]
    
  }
  
  else if(seasonality == "low"){
    
    low.seas.years = abs.changes < quantile(abs.changes, probs = seas.threshold)
    low.seas.years = years[low.seas.years]
    
  }
  
  
  ## STARTING SIMULATIONS ##
  
  for(i in 1:nsimul){
    
    # Displaying simulations progress
    if(nsimul >= 20){
      
      if((i %% round(nsimul / 20)) == 0) cat(i, "simulations\n")
    }
    
    
    # Start every simulation with the same population vector
    n1 = n0
    
    for(j in 1:nyears){
      
      # Keep population vector of the beginning of the cycle
      pop.vec1 = n1
      
      # Random sampling among all years (Case (1))
      if(seasonality == "control"){
        
        year = sample(years, 1)
      }
      
      # Random sampling among high years only (Case (2))
      else if(seasonality == "high"){
        
        year = sample(high.seas.years, 1)
      }
      
      # Random sampling among season-non-significant years only (Case (3))
      else if(seasonality == "low"){
        
        year = sample(low.seas.years, 1)
      }
      
      # Building seasonal matrices
      summer.matrix = marmots.buildmat.summer(marmots.vr, year)
      winter.matrix = marmots.buildmat.winter(marmots.vr, year)
      
      # Computing the annual matrix
      annual.matrix = summer.matrix %*% winter.matrix
      
      # Computing the new population vector
      n1 = annual.matrix %*% n1
      
      # Keep population vector at the end of the cycle
      pop.vec2 = n1
      
      # Compute and save one-time log lambda (in the log.lambda array)
      log.lambda[i, j] = log(sum(pop.vec2) / sum(pop.vec1))
      
      # Record if extinction occurs in the simulation
      if(sum(n1) < 6 | n1[4] < 1){
        
        ext[i] = 1
        break
      }
    }
  }
  return(list(log.lambda, ext))
}




###########################################################################
#
# 3. Performing the simulations ----
#
###########################################################################

# We look at the effect of high vs low seasonality in vital rates for which the random effect structure is 1+season|year (i.e., random effect of the year on the intercept and the difference between seasons). These vital rates are: yearling, non reproductive, and reproductive adult survival. 

n0 = c(12, 11, 9, 7) # Initial population vector


## 3.1. Tests ----
# -----------

test1 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "control")
test2 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "high", vr.model = survY)
test3 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "low", vr.model = survY)
test4 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "high", vr.model = survR)
test5 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "low", vr.model = survR)
test6 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "high", vr.model = survNR)
test7 = stoch.sim.marmot(nsimul = 10, nyears = 10, n0 = n0, seasonality = "low", vr.model = survNR)


## 3.2.Control case ----
# -----------------

print("CONTROL CASE")
control.sim = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "control")

# Simulation-wise stochastic lambdas
control.sim.lambda = apply(control.sim[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
control.stoch.lambda = mean(control.sim.lambda, na.rm = T)

# Variance in log lambda across years
control.var.lambda = apply(control.sim[[1]], 1, var, na.rm = T)

# Mean extinction probability
control.ext.prob = mean(control.sim[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(control.sim[[1]]), type = "l", col = rainbow(nrow(control.sim[[1]])), main = "Simulations yearly log lambda - Control case", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(control.sim[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3. High seasonality case ----
# ---------------------------

## 3.3.1. Y survival as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - Y SURVIVAL")
high.seas.Y.surv = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "high", seas.threshold = 0.5, vr.model = survY)

# Simulation-wise stochastic lambdas
high.seas.Y.surv.sim.lambda = apply(high.seas.Y.surv[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.Y.surv.stoch.lambda = mean(high.seas.Y.surv.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.Y.surv.var.lambda = apply(high.seas.Y.surv[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.Y.surv.ext.prob = mean(high.seas.Y.surv[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.Y.surv[[1]]), type = "l", col = rainbow(nrow(high.seas.Y.surv[[1]])), main = "Simulations yearly log lambda - High seasonality case (Y survival)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.Y.surv[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.2. NR survival as reference ----
# ---------------------------------

print("HIGH SEASONALITY CASE - NR SURVIVAL")
high.seas.NR.surv = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "high", seas.threshold = 0.5, vr.model = survNR)

# Simulation-wise stochastic lambdas
high.seas.NR.surv.sim.lambda = apply(high.seas.NR.surv[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.NR.surv.stoch.lambda = mean(high.seas.NR.surv.sim.lambda, na.rm = T)

# Variance in log lambda across years 
high.seas.NR.surv.var.lambda = apply(high.seas.NR.surv[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.NR.surv.ext.prob = mean(high.seas.NR.surv[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.NR.surv[[1]]), type = "l", col = rainbow(nrow(high.seas.NR.surv[[1]])), main = "Simulations yearly log lambda - High seasonality case (NR survival)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.NR.surv[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.3. R survival as reference ----
# --------------------------------

print("HIGH SEASONALITY CASE - R SURVIVAL")
high.seas.R.surv = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "high", seas.threshold = 0.5, vr.model = survR)

# Simulation-wise stochastic lambdas
high.seas.R.surv.sim.lambda = apply(high.seas.R.surv[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.R.surv.stoch.lambda = mean(high.seas.R.surv.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.R.surv.var.lambda = apply(high.seas.R.surv[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.R.surv.ext.prob = mean(high.seas.R.surv[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.R.surv[[1]]), type = "l", col = rainbow(nrow(high.seas.R.surv[[1]])), main = "Simulations yearly log lambda - High seasonality case (R survival)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.R.surv[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.4. Low seasonality case ----
# --------------------------

## 3.4.1. Y survival as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - Y SURVIVAL")
low.seas.Y.surv = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "low", seas.threshold = 0.5, vr.model = survY)

# Simulation-wise stochastic lambdas
low.seas.Y.surv.sim.lambda = apply(low.seas.Y.surv[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.Y.surv.stoch.lambda = mean(low.seas.Y.surv.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.Y.surv.var.lambda = apply(low.seas.Y.surv[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.Y.surv.ext.prob = mean(low.seas.Y.surv[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.Y.surv[[1]]), type = "l", col = rainbow(nrow(low.seas.Y.surv[[1]])), main = "Simulations yearly log lambda - Low seasonality case (Y survival)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.Y.surv[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.4.2. NR survival as reference ----
# --------------------------------

print("LOW SEASONALITY CASE - NR SURVIVAL")
low.seas.NR.surv = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "low", seas.threshold = 0.5, vr.model = survNR)

# Simulation-wise stochastic lambdas
low.seas.NR.surv.sim.lambda = apply(low.seas.NR.surv[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.NR.surv.stoch.lambda = mean(low.seas.NR.surv.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.NR.surv.var.lambda = apply(low.seas.NR.surv[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.NR.surv.ext.prob = mean(low.seas.NR.surv[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.NR.surv[[1]]), type = "l", col = rainbow(nrow(low.seas.NR.surv[[1]])), main = "Simulations yearly log lambda - Low seasonality case (NR survival)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.NR.surv[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.4.3. R survival as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - R SURVIVAL")
low.seas.R.surv = stoch.sim.marmot(nsimul = 500, nyears = 100, n0 = n0, seasonality = "low", seas.threshold = 0.5, vr.model = survR)

# Simulation-wise stochastic lambdas
low.seas.R.surv.sim.lambda = apply(low.seas.R.surv[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.R.surv.stoch.lambda = mean(low.seas.R.surv.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.R.surv.var.lambda = apply(low.seas.R.surv[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.R.surv.ext.prob = mean(low.seas.R.surv[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.R.surv[[1]]), type = "l", col = rainbow(nrow(low.seas.R.surv[[1]])), main= "Simulations yearly log lambda - Low seasonality case (R survival)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.R.surv[[1]], 2, mean, na.rm = T), lwd = 2)




###########################################################################
#
# 4. Saving results ----
#
###########################################################################

save.image("MarmotsProjectionsResults.RData")

