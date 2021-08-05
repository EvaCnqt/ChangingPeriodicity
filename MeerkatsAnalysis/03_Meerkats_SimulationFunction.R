##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script contains the function needed to project the meerkat population dynamics under scenarios
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
  library(lubridate)
}

load.librairies()


## 1.3. Loading data ----
# ------------------

source("02_Meerkats_MPMs.R")

data.meerkats = read.csv("MeerkatsData.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

meerkats.vr = read.csv("Meerkats_VitalRatesEstimates.csv")

load("GLMM_poprange.RData")




###########################################################################
#
# 2. Projection function
#
###########################################################################

stoch.sim.meerkat <- function(nsimul = 500,
                              nyears = 100,
                              nseas = 2,
                              n0,
                              dens0,
                              vr.model = NULL,
                              seas.threshold = 0.5,
                              density = T,
                              seasonality = "control"){
  
  # (1) Control case = Simulating what happens in the normal conditions, with periodic vital rates and stochastic variation across the years - Random sampling of the year among all available ones (1998 to 2016).
  # (2) High seasonality case = Random sampling only among years with a significant difference between seasons in season-dependent vital rates (adults survival).
  # (3) Low seasonality case = Random sampling only among years with a non-significant difference between seasons in season-dependent vital rates (adults survival).
  
  # Error if seasonality is either yes or no but no model has been given as a reference to find significant years
  
  if(is.null(vr.model)&seasonality != "control"){
    
    stop("Please enter the reference vital rate model to determine significant/nonsignificant years for seasonality!")
  
  }
  
  # Storage arrays for log lambdas, population extinction and density
  
  log.lambda = array(NA, c(nsimul, nyears))
  ext = seq(0, 0, length.out = nsimul)
  dens = array(NA, c(nsimul, nyears * nseas))
  
  # If specific seasonality significance, find related years: take significative slope at the response scale, not link scale
  
  years = meerkats.vr$year
  
  if(!is.null(vr.model)){
  
    # Compute absolute changes
  
      # Recruitment
  
    if(formula(vr.model) == formula(recruitment.D) | formula(vr.model) == formula(recruitment.H)){
    
      abs.changes = abs(exp(coef(vr.model)$`(Intercept)` + coef(vr.model)$seasonrain) - exp(coef(vr.model)$`(Intercept)`))
    
    }
  
      # Survival or transition 
  
    else{
    
      abs.changes = abs(inv.logit(coef(vr.model)$year$`(Intercept)` + coef(vr.model)$year$seasonrain) - inv.logit(coef(vr.model)$year$`(Intercept)`))
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
  }
  
  ## STARTING SIMULATIONS ##
  
  for(i in 1:nsimul){
    
    # Displaying simulations progress
    if(nsimul >= 20){
      
      if((i %% round(nsimul / 20)) == 0) cat(i, "simulations\n")
    }
    
    time.step = 1
    
    # Start every simulation with the same population vector and density
    n1 = n0
    
    # If density dependence is not required, set density to the average density, which is the optimum density of the population
    
    if(density == F) dens1 = mean(unique(data.meerkats$density))

    else dens1 = dens0 
    
    for(j in 1:nyears){
      
      # Keep population vector of the beginning of the cycle
      pop.vec1 = n1
      
      # Random sampling among all years (Cases (1))
      
      if(seasonality == "control") year = sample(years, 1)
      
      # Random sampling among season-significant years only (Case (2))
      else if(seasonality == "high") year = sample(high.seas.years, 1)
      
      # Random sampling among season-non-significant years only (Case (3))
      else if(seasonality == "low") year = sample(low.seas.years, 1)
      
      
      for(season in c("dry", "rain")){
        
        # Build seasonal MPM and obtain the new population vector n1
        
        if(season == "dry"){
          
          dry.matrix = meerkats.buildmat(year, season = season, dens1)
          n1 = dry.matrix %*% as.numeric(n1)
        }
        
        else{
          
          rain.matrix = meerkats.buildmat(year, season = season, dens1)
          n1 = rain.matrix %*% as.numeric(n1)
        }
      
        
        # If density dependence is required, predict new population range from year and number of dominants (we use exp() because the response variable of the model was on the log scale)
        
        if(density == T){
          
          new.pop.range = as.numeric(exp(predict(pop.range, newdata = data.frame(year = year, domNum = n1[4]))))
          
          # Computing the new density value. We use the size of the population times 2 because we want the size of the whole population, not only the females.
          
          dens1 = (2 * sum(n1)) / new.pop.range
          
          # if(density == F & is.infinite(dens1)){
          #   
          #   break
          # }
          
          dens[i, time.step] = dens1
        }
        
        time.step = time.step + 1
      }
      
      # Keep population vector of the end of the cycle
      pop.vec2 = n1
      
      # Calculate and save one-time log lambda (in the log.lambda array)
      log.lambda[i, j] = log(sum(pop.vec2)/sum(pop.vec1))
      
      # Record if extinction occurs:
      if(sum(n1) < 20 | n1[4] < 5){
        
        ext[i] = 1
        break
      }
    }
  }
  
  return(list(log.lambda, ext, dens))
}




###########################################################################
#
# 3. Performing the simulations ----
#
###########################################################################

# We look at the effect of high vs low seasonality in vital rates for which the random effect structure is 1+season|year (i.e., random effect of the year on the intercept and the difference between seasons). These vital rates are: subadult, helper, and dominant survival, helper emigration, transition from helper to dominant, and helper and dominant recruitment.

# We start with intial population vector and density picked from 2016 dry season
nbJ.init = length(data.meerkats$stage[data.meerkats$stage == "J" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])
nbS.init = length(data.meerkats$stage[data.meerkats$stage == "S" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])
nbH.init = length(data.meerkats$stage[data.meerkats$stage == "H" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])
nbD.init = length(data.meerkats$stage[data.meerkats$stage == "D" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])

n0 = c(nbJ.init, nbS.init, nbH.init, nbD.init)

dens0 = unique(data.meerkats$density[data.meerkats$year == "2016" & data.meerkats$season == "dry"])


## 3.1. Tests ----
# -----------

test1 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "control")
test2 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = survH)
test2a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", density = F, vr.model = survH)

test3 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = survH)
test3a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", density = F, vr.model = survH)

test4 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = emig)
test4a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", density = F, vr.model = emig)

test5 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = emig)
test5a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", density = F, vr.model = emig)

test6 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = recruitment.H)
test6a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", density = F, vr.model = recruitment.H)

test7 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = recruitment.H)
test7a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", density = F, vr.model = recruitment.H)

test8 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = recruitment.D)
test8a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", density = F, vr.model = recruitment.D)

test9 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = recruitment.D)
test9a = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", density = F, vr.model = recruitment.D)


## 3.2. Control case ----
# ------------------

## 3.2.1. With density dependence ----
# -------------------------------

print("CONTROL CASE - WITH DENSITY DEPENDENCE")
control.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "control")

# Simulation-wise stochastic lambdas
control.dens.dep.sim.lambda = apply(control.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
control.dens.dep.stoch.lambda = mean(control.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
control.dens.dep.var.lambda = apply(control.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
control.dens.dep.ext.prob = mean(control.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(control.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(control.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(control.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(control.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(control.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(control.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.2.2. With optimum density ----
# ----------------------------

print("CONTROL CASE - WITH OPTIMUM DENSITY")
control.sim.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "control")

# Simulation-wise stochastic lambdas
control.opt.dens.sim.lambda = apply(control.sim.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
control.opt.dens.stoch.lambda = mean(control.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
control.opt.dens.var.lambda = apply(control.sim.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
control.opt.dens.ext.prob = mean(control.sim.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(control.sim.opt.dens[[1]]), type = "l", col = rainbow(nrow(control.sim.opt.dens[[1]])), main = "Simulations yearly log lambda - Control scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(control.sim.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3. High seasonality case ----
# --------------------------- 

## 3.3.1. S survival as reference ----
# -------------------------------

## 3.3.1.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - S SURVIVAL - WITH DENSITY DEPENDENCE")
high.seas.S.surv.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = survS)

# Simulation-wise stochastic lambdas
high.seas.S.surv.dens.dep.sim.lambda = apply(high.seas.S.surv.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.S.surv.dens.dep.stoch.lambda = mean(high.seas.S.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.S.surv.dens.dep.var.lambda = apply(high.seas.S.surv.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.S.surv.dens.dep.ext.prob = mean(high.seas.S.surv.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.S.surv.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.S.surv.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in subadults survival scenario - Density dependence", xlab  = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.S.surv.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.S.surv.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.S.surv.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in subadults survival scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.S.surv.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.1.2. With optimum density ----
# ------------------------------

print("HIGH SEASONALITY CASE - S SURVIVAL - WITH OPTIMUM DENSITY")
high.seas.S.surv.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5, vr.model = survS)

# Simulation-wise stochastic lambdas
high.seas.S.surv.opt.dens.sim.lambda = apply(high.seas.S.surv.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.S.surv.opt.dens.stoch.lambda = mean(high.seas.S.surv.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.S.surv.opt.dens.var.lambda = apply(high.seas.S.surv.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.S.surv.opt.dens.ext.prob = mean(high.seas.S.surv.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.S.surv.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.S.surv.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in subadults survival scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.S.surv.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)



## 3.3.2. H survival as reference ----
# -------------------------------

## 3.3.2.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - H SURVIVAL - WITH DENSITY DEPENDENCE")
high.seas.H.surv.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = survH)

# Simulation-wise stochastic lambdas
high.seas.H.surv.dens.dep.sim.lambda = apply(high.seas.H.surv.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.H.surv.dens.dep.stoch.lambda = mean(high.seas.H.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.H.surv.dens.dep.var.lambda = apply(high.seas.H.surv.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.H.surv.dens.dep.ext.prob = mean(high.seas.H.surv.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.H.surv.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.H.surv.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers survival scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.H.surv.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.H.surv.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.H.surv.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in helpers survival scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.H.surv.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.2.2. With optimum density ----
# ------------------------------

print("HIGH SEASONALITY CASE - H SURVIVAL - WITH OPTIMUM DENSITY")
high.seas.H.surv.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5, vr.model = survH)

# Simulation-wise stochastic lambdas
high.seas.H.surv.opt.dens.sim.lambda = apply(high.seas.H.surv.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.H.surv.opt.dens.stoch.lambda = mean(high.seas.H.surv.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.H.surv.opt.dens.var.lambda = apply(high.seas.H.surv.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.H.surv.opt.dens.ext.prob = mean(high.seas.H.surv.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.H.surv.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.H.surv.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers survival scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.H.surv.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.3. D survival as reference ----
# -------------------------------

## 3.3.3.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - D SURVIVAL - WITH DENSITY DEPENDENCE")
high.seas.D.surv.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = survD)

# Simulation-wise stochastic lambdas
high.seas.D.surv.dens.dep.sim.lambda = apply(high.seas.D.surv.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.D.surv.dens.dep.stoch.lambda = mean(high.seas.D.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.D.surv.dens.dep.var.lambda = apply(high.seas.D.surv.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.D.surv.dens.dep.ext.prob = mean(high.seas.D.surv.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.D.surv.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.D.surv.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in dominants survival scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.D.surv.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.D.surv.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.D.surv.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in dominants survival scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.D.surv.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.3.2. With optimum density ----
# ------------------------------

print("HIGH SEASONALITY CASE - D SURVIVAL - WITH OPTIMUM DENSITY")
high.seas.D.surv.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5,vr.model = survD)

# Simulation-wise stochastic lambdas
high.seas.D.surv.opt.dens.sim.lambda = apply(high.seas.D.surv.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.D.surv.opt.dens.stoch.lambda=mean(high.seas.D.surv.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.D.surv.opt.dens.var.lambda = apply(high.seas.D.surv.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.D.surv.opt.dens.ext.prob = mean(high.seas.D.surv.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.D.surv.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.D.surv.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in dominants survival scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.D.surv.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.4. Emigration as reference ----
# -------------------------------

## 3.3.4.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - EMIGRATION - WITH DENSITY DEPENDENCE")
high.seas.emig.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
high.seas.emig.dens.dep.sim.lambda = apply(high.seas.emig.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.emig.dens.dep.stoch.lambda = mean(high.seas.emig.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.emig.dens.dep.var.lambda = apply(high.seas.emig.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.emig.dens.dep.ext.prob = mean(high.seas.emig.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.emig.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.emig.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers emigration scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.emig.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.emig.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.emig.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in helpers emigration scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.emig.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.4.2. With optimum density ----
# ------------------------------

print("HIGH SEASONALITY CASE - EMIGRATION - WITH OPTIMUM DENSITY")
high.seas.emig.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
high.seas.emig.opt.dens.sim.lambda = apply(high.seas.emig.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.emig.opt.dens.stoch.lambda = mean(high.seas.emig.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.emig.opt.dens.var.lambda = apply(high.seas.emig.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.emig.opt.dens.ext.prob = mean(high.seas.emig.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.emig.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.emig.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers emigration scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.emig.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.5. Transition as reference ----
# -------------------------------

## 3.3.5.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - TRANSITION - WITH DENSITY DEPENDENCE")
high.seas.trans.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
high.seas.trans.dens.dep.sim.lambda = apply(high.seas.trans.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.trans.dens.dep.stoch.lambda = mean(high.seas.trans.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.trans.dens.dep.var.lambda = apply(high.seas.trans.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.trans.dens.dep.ext.prob = mean(high.seas.trans.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.trans.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.trans.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers transition scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.trans.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.trans.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.trans.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in helpers transition scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.trans.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.5.2. With optimum density ----
# ------------------------------

print("HIGH SEASONALITY CASE - TRANSITION - WITH OPTIMUM DENSITY")
high.seas.trans.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
high.seas.trans.opt.dens.sim.lambda = apply(high.seas.trans.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.trans.opt.dens.stoch.lambda = mean(high.seas.trans.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.trans.opt.dens.var.lambda = apply(high.seas.trans.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.trans.opt.dens.ext.prob = mean(high.seas.trans.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.trans.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.trans.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers transition scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.trans.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.6. H recruitment as reference ----
# ----------------------------------

## 3.3.6.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - H RECRUITMENT - WITH DENSITY DEPENDENCE")
high.seas.h.recruit.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = recruitment.H)

# Simulation-wise stochastic lambdas
high.seas.h.recruit.dens.dep.sim.lambda = apply(high.seas.h.recruit.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.h.recruit.dens.dep.stoch.lambda = mean(high.seas.h.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.h.recruit.dens.dep.var.lambda = apply(high.seas.h.recruit.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.h.recruit.dens.dep.ext.prob = mean(high.seas.h.recruit.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.h.recruit.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.h.recruit.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers recruitment scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.h.recruit.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.h.recruit.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.h.recruit.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in helpers recruitment scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.h.recruit.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.6.2. With optimum density ----
# ------------------------------

print("HIGH SEASONALITY CASE - H RECRUITMENT - WITH OPTIMUM DENSITY")
high.seas.h.recruit.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5, vr.model = recruitment.H)

# Simulation-wise stochastic lambdas
high.seas.h.recruit.opt.dens.sim.lambda = apply(high.seas.h.recruit.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.h.recruit.opt.dens.stoch.lambda = mean(high.seas.h.recruit.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.h.recruit.opt.dens.var.lambda = apply(high.seas.h.recruit.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.h.recruit.opt.dens.ext.prob = mean(high.seas.h.recruit.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.h.recruit.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.h.recruit.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in helpers recruitment scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.h.recruit.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)


## 3.3.7. D recruitment as reference ----
# ----------------------------------

## 3.3.7.1. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - D RECRUITMENT - WITH DENSITY DEPENDENCE")
high.seas.d.recruit.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", seas.threshold = 0.5, vr.model = recruitment.D)

# Simulation-wise stochastic lambdas
high.seas.d.recruit.dens.dep.sim.lambda = apply(high.seas.d.recruit.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.d.recruit.dens.dep.stoch.lambda = mean(high.seas.d.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.d.recruit.dens.dep.var.lambda = apply(high.seas.d.recruit.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.d.recruit.dens.dep.ext.prob = mean(high.seas.d.recruit.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.d.recruit.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.d.recruit.dens.dep[[1]])), main = "Simulations yearly log lambda - High seasonality in dominants recruitment scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.d.recruit.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.d.recruit.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.d.recruit.dens.dep[[3]])), main = "Density over seasonal time steps - High seasonality in dominants recruitment scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.d.recruit.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.7.2. With density dependence ----
# ---------------------------------

print("HIGH SEASONALITY CASE - D RECRUITMENT - WITH OPTIMUM DENSITY")
high.seas.d.recruit.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "high", seas.threshold = 0.5, vr.model = recruitment.D)

# Simulation-wise stochastic lambdas
high.seas.d.recruit.opt.dens.sim.lambda = apply(high.seas.d.recruit.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.d.recruit.opt.dens.stoch.lambda = mean(high.seas.d.recruit.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.d.recruit.opt.dens.var.lambda = apply(high.seas.d.recruit.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.d.recruit.opt.dens.ext.prob = mean(high.seas.d.recruit.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.d.recruit.opt.dens[[1]]), type = "l", col = rainbow(nrow(high.seas.d.recruit.opt.dens[[1]])), main = "Simulations yearly log lambda - High seasonality in dominants recruitment scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.d.recruit.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)




## 3.4. Low seasonality case ----
# --------------------------

## 3.4.1. S survival as reference ----
# -------------------------------

## 3.4.1.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - S SURVIVAL - WITH DENSITY DEPENDENCE")
low.seas.S.surv.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = survS)

# Simulation-wise stochastic lambdas
low.seas.S.surv.dens.dep.sim.lambda = apply(low.seas.S.surv.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.S.surv.dens.dep.stoch.lambda = mean(low.seas.S.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.S.surv.dens.dep.var.lambda = apply(low.seas.S.surv.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.S.surv.dens.dep.ext.prob = mean(low.seas.S.surv.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.S.surv.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.S.surv.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in subadults survival scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.S.surv.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.S.surv.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.S.surv.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in subadults survival scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.S.surv.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.1.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - S SURVIVAL - WITH OPTIMUM DENSITY")
low.seas.S.surv.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = survS)

# Simulation-wise stochastic lambdas
low.seas.S.surv.opt.dens.sim.lambda = apply(low.seas.S.surv.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.S.surv.opt.dens.stoch.lambda = mean(low.seas.S.surv.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.S.surv.opt.dens.var.lambda = apply(low.seas.S.surv.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.S.surv.opt.dens.ext.prob = mean(low.seas.S.surv.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.S.surv.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.S.surv.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in subadults survival scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.S.surv.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)



## 3.4.2. H survival as reference ----
# -------------------------------

## 3.4.2.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - H SURVIVAL - WITH DENSITY DEPENDENCE")
low.seas.H.surv.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = survH)

# Simulation-wise stochastic lambdas
low.seas.H.surv.dens.dep.sim.lambda = apply(low.seas.H.surv.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.H.surv.dens.dep.stoch.lambda = mean(low.seas.H.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.H.surv.dens.dep.var.lambda = apply(low.seas.H.surv.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.H.surv.dens.dep.ext.prob = mean(low.seas.H.surv.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.H.surv.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.H.surv.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers survival scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.H.surv.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.H.surv.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.H.surv.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in helpers survival scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.H.surv.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.2.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - H SURVIVAL - WITH OPTIMUM DENSITY")
low.seas.H.surv.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = survH)

# Simulation-wise stochastic lambdas
low.seas.H.surv.opt.dens.sim.lambda = apply(low.seas.H.surv.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.H.surv.opt.dens.stoch.lambda = mean(low.seas.H.surv.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.H.surv.opt.dens.var.lambda = apply(low.seas.H.surv.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.H.surv.opt.dens.ext.prob = mean(low.seas.H.surv.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.H.surv.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.H.surv.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers survival scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.H.surv.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)



## 3.4.3. D survival as reference ----
# -------------------------------

## 3.4.3.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - D SURVIVAL - WITH DENSITY DEPENDENCE")
low.seas.D.surv.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = survD)

# Simulation-wise stochastic lambdas
low.seas.D.surv.dens.dep.sim.lambda = apply(low.seas.D.surv.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.D.surv.dens.dep.stoch.lambda = mean(low.seas.D.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.D.surv.dens.dep.var.lambda = apply(low.seas.D.surv.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.D.surv.dens.dep.ext.prob = mean(low.seas.D.surv.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.D.surv.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.D.surv.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in dominants survival scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.D.surv.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.D.surv.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.D.surv.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in dominants survival scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.D.surv.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.3.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - D SURVIVAL - WITH OPTIMUM DENSITY")
low.seas.D.surv.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = survD)

# Simulation-wise stochastic lambdas
low.seas.D.surv.opt.dens.sim.lambda = apply(low.seas.D.surv.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.D.surv.opt.dens.stoch.lambda = mean(low.seas.D.surv.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.D.surv.opt.dens.var.lambda = apply(low.seas.D.surv.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.D.surv.opt.dens.ext.prob = mean(low.seas.D.surv.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.D.surv.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.D.surv.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in dominants survival scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.D.surv.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)



## 3.4.4. Emigration as reference ----
# -------------------------------

## 3.4.4.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - EMIGRATION - WITH DENSITY DEPENDENCE")
low.seas.emig.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
low.seas.emig.dens.dep.sim.lambda = apply(low.seas.emig.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.emig.dens.dep.stoch.lambda = mean(low.seas.emig.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.emig.dens.dep.var.lambda = apply(low.seas.emig.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.emig.dens.dep.ext.prob = mean(low.seas.emig.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.emig.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.emig.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers emigration scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.emig.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.emig.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.emig.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in helpers emigration scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.emig.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.4.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - EMIGRATION - WITH OPTIMUM DENSITY")
low.seas.emig.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
low.seas.emig.opt.dens.sim.lambda = apply(low.seas.emig.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.emig.opt.dens.stoch.lambda = mean(low.seas.emig.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.emig.opt.dens.var.lambda = apply(low.seas.emig.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.emig.opt.dens.ext.prob = mean(low.seas.emig.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.emig.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.emig.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers emigration scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.emig.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)




## 3.4.5. Transition as reference ----
# -------------------------------

## 3.4.5.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - TRANSITION - WITH DENSITY DEPENDENCE")
low.seas.trans.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
low.seas.trans.dens.dep.sim.lambda = apply(low.seas.trans.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.trans.dens.dep.stoch.lambda = mean(low.seas.trans.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.trans.dens.dep.var.lambda = apply(low.seas.trans.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.trans.dens.dep.ext.prob = mean(low.seas.trans.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.trans.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.trans.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers transition scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.trans.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.trans.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.trans.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in helpers transition scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.trans.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.5.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - TRANSITION - WITH OPTIMUM DENSITY")
low.seas.trans.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = emig)

# Simulation-wise stochastic lambdas
low.seas.trans.opt.dens.sim.lambda = apply(low.seas.trans.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.trans.opt.dens.stoch.lambda = mean(low.seas.trans.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.trans.opt.dens.var.lambda = apply(low.seas.trans.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.trans.opt.dens.ext.prob = mean(low.seas.trans.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.trans.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.trans.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers transition scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.trans.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)



## 3.4.6. H recruitment as reference ----
# ----------------------------------

## 3.4.6.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - H RECRUITMENT - WITH DENSITY DEPENDENCE")
low.seas.h.recruit.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = recruitment.H)

# Simulation-wise stochastic lambdas
low.seas.h.recruit.dens.dep.sim.lambda = apply(low.seas.h.recruit.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.h.recruit.dens.dep.stoch.lambda = mean(low.seas.h.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.h.recruit.dens.dep.var.lambda = apply(low.seas.h.recruit.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.h.recruit.dens.dep.ext.prob = mean(low.seas.h.recruit.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.h.recruit.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.h.recruit.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers recruitment scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.h.recruit.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.h.recruit.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.h.recruit.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in helpers recruitment scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.h.recruit.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.6.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - H RECRUITMENT - WITH OPTIMUM DENSITY")
low.seas.h.recruit.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = recruitment.H)

# Simulation-wise stochastic lambdas
low.seas.h.recruit.opt.dens.sim.lambda = apply(low.seas.h.recruit.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.h.recruit.opt.dens.stoch.lambda = mean(low.seas.h.recruit.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.h.recruit.opt.dens.var.lambda = apply(low.seas.h.recruit.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.h.recruit.opt.dens.ext.prob = mean(low.seas.h.recruit.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.h.recruit.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.h.recruit.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in helpers recruitment scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.h.recruit.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)




## 3.4.7. D recruitment as reference ----
# ----------------------------------

## 3.4.7.1. With density dependence ----
# ---------------------------------

print("LOW SEASONALITY CASE - D RECRUITMENT - WITH DENSITY DEPENDENCE")
low.seas.d.recruit.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", seas.threshold = 0.5, vr.model = recruitment.D)

# Simulation-wise stochastic lambdas
low.seas.d.recruit.dens.dep.sim.lambda = apply(low.seas.d.recruit.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.d.recruit.dens.dep.stoch.lambda = mean(low.seas.d.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.d.recruit.dens.dep.var.lambda = apply(low.seas.d.recruit.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.d.recruit.dens.dep.ext.prob = mean(low.seas.d.recruit.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.d.recruit.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.d.recruit.dens.dep[[1]])), main = "Simulations yearly log lambda - Low seasonality in dominants recruitment scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.d.recruit.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.d.recruit.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.d.recruit.dens.dep[[3]])), main = "Density over seasonal time steps - Low seasonality in dominants recruitment scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.d.recruit.dens.dep[[3]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.4.7.2. With optimum density ----
# ------------------------------

print("LOW SEASONALITY CASE - D RECRUITMENT - WITH OPTIMUM DENSITY")
low.seas.d.recruit.opt.dens = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, density = F, seasonality = "low", seas.threshold = 0.5, vr.model = recruitment.D)

# Simulation-wise stochastic lambdas
low.seas.d.recruit.opt.dens.sim.lambda = apply(low.seas.d.recruit.opt.dens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.d.recruit.opt.dens.stoch.lambda = mean(low.seas.d.recruit.opt.dens.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.d.recruit.opt.dens.var.lambda = apply(low.seas.d.recruit.opt.dens[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.d.recruit.opt.dens.ext.prob = mean(low.seas.d.recruit.opt.dens[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.d.recruit.opt.dens[[1]]), type = "l", col = rainbow(nrow(low.seas.d.recruit.opt.dens[[1]])), main = "Simulations yearly log lambda - Low seasonality in dominants recruitment scenario - Optimum density", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.d.recruit.opt.dens[[1]], 2, mean, na.rm = T), lwd = 2)




###########################################################################
#
# 4. Saving files and data
#
###########################################################################

save.image("MeerkatsProjectionsResults.RData")
