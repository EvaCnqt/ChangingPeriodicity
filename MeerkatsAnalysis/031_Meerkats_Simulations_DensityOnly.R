##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al. 2022, Ecology).
# 
# Eva Conquet, Arpat Ozgul, Daniel T. Blumstein, Kenneth B. Armitage,
# Madan K. Oli, Julien G. A. Martin, Tim H. Clutton-Brock, and Maria Paniw
#
# This script contains the function needed to project the meerkat population dynamics under scenarios
# of higher or lower seasonality in vital rates. 
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

source("021_Meerkats_MPMs_DensityOnly.R")

data.meerkats = read.csv("MeerkatsData.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

meerkats.vr = read.csv("Meerkats_VitalRatesEstimates_DensityOnly.csv")

load("GLMM_poprange_DensityOnly.RData")




###########################################################################
#
# 2. Projection function ----
#
###########################################################################

stoch.sim.meerkat <- function(nsimul = 500,
                              nyears = 100,
                              nseas = 2,
                              n0,
                              dens0,
                              vr.model = NULL,
                              seas.threshold = 0.5,
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
    
    if(formula(vr.model) == formula(recruitment.D_WithSeason) | formula(vr.model) == formula(recruitment.H_WithSeason)){
      
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
    
    dens1 = dens0
    
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
          
          dry.matrix = meerkats.buildmat(year, density = dens1)
          n1 = dry.matrix %*% as.numeric(n1)
        }
        
        else{
          
          rain.matrix = meerkats.buildmat(year, density = dens1)
          n1 = rain.matrix %*% as.numeric(n1)
        }
        
        # If density dependence is required, predict new population range from year and number of dominants (we use exp() because the response variable of the model was on the log scale)
        
        new.pop.range = as.numeric(exp(predict(pop.range, newdata = data.frame(year = year, domNum = n1[4]))))
        
        # Computing the new density value. We use the size of the population times 2 because we want the size of the whole population, not only the females.
        
        dens1 = (2 * sum(n1)) / new.pop.range
        
        # if(density == F & is.infinite(dens1)){
        #   
        #   break
        # }
        
        dens[i, time.step] = dens1
        
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

test2 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = survH_WithSeason)

test3 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = survH_WithSeason)

test4 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = emig_WithSeason)

test5 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = emig_WithSeason)

test6 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = recruitment.H_WithSeason)

test7 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = recruitment.H_WithSeason)

test8 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "high", vr.model = recruitment.D_WithSeason)

test9 = stoch.sim.meerkat(nsimul = 10, nyears = 10, n0 = n0, dens0 = dens0, seasonality = "low", vr.model = recruitment.D_WithSeason)


## 3.2. Control case ----
# ------------------

print("CONTROL CASE - WITH DENSITY DEPENDENCE")
control.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0)

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


## 3.3. High seasonality case ----
# ---------------------------

## 3.3.1. S survival as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - S SURVIVAL")
high.seas.S.surv.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = survS_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.S.surv.dens.dep.sim.lambda = apply(high.seas.S.surv.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.S.surv.dens.dep.stoch.lambda = mean(high.seas.S.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.S.surv.dens.dep.var.lambda = apply(high.seas.S.surv.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.S.surv.dens.dep.ext.prob = mean(high.seas.S.surv.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.S.surv.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.S.surv.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.S.surv.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.S.surv.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.S.surv.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.S.surv.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.2. H survival as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - H SURVIVAL")
high.seas.H.surv.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = survH_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.H.surv.dens.dep.sim.lambda = apply(high.seas.H.surv.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.H.surv.dens.dep.stoch.lambda = mean(high.seas.H.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.H.surv.dens.dep.var.lambda = apply(high.seas.H.surv.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.H.surv.dens.dep.ext.prob = mean(high.seas.H.surv.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.H.surv.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.H.surv.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.H.surv.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.H.surv.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.H.surv.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.H.surv.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step



## 3.3.3. D survival as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - D SURVIVAL")
high.seas.D.surv.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = survD_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.D.surv.dens.dep.sim.lambda = apply(high.seas.D.surv.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.D.surv.dens.dep.stoch.lambda = mean(high.seas.D.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.D.surv.dens.dep.var.lambda = apply(high.seas.D.surv.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.D.surv.dens.dep.ext.prob = mean(high.seas.D.surv.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.D.surv.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.D.surv.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.D.surv.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.D.surv.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.D.surv.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.D.surv.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.4. Emigration as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - EMIGRATION")
high.seas.emig.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = emig_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.emig.dens.dep.sim.lambda = apply(high.seas.emig.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.emig.dens.dep.stoch.lambda = mean(high.seas.emig.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.emig.dens.dep.var.lambda = apply(high.seas.emig.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.emig.dens.dep.ext.prob = mean(high.seas.emig.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.emig.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.emig.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.emig.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.emig.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.emig.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.emig.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.5. Transition as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - TRANSITION")
high.seas.trans.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = transition_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.trans.dens.dep.sim.lambda = apply(high.seas.trans.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.trans.dens.dep.stoch.lambda = mean(high.seas.trans.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.trans.dens.dep.var.lambda = apply(high.seas.trans.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.trans.dens.dep.ext.prob = mean(high.seas.trans.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.trans.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.trans.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.trans.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.trans.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.trans.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.trans.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.6. H recruitment as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - H RECRUITMENT")
high.seas.h.recruit.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = recruitment.H_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.h.recruit.dens.dep.sim.lambda = apply(high.seas.h.recruit.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.h.recruit.dens.dep.stoch.lambda = mean(high.seas.h.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.h.recruit.dens.dep.var.lambda = apply(high.seas.h.recruit.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.h.recruit.dens.dep.ext.prob = mean(high.seas.h.recruit.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.h.recruit.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.h.recruit.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.h.recruit.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.h.recruit.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.h.recruit.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.h.recruit.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.3.7. D recruitment as reference ----
# -------------------------------

print("HIGH SEASONALITY CASE - D RECRUITMENT")
high.seas.d.recruit.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "high", vr = recruitment.D_WithSeason)

# Simulation-wise stochastic lambdas
high.seas.d.recruit.dens.dep.sim.lambda = apply(high.seas.d.recruit.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
high.seas.d.recruit.dens.dep.stoch.lambda = mean(high.seas.d.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
high.seas.d.recruit.dens.dep.var.lambda = apply(high.seas.d.recruit.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
high.seas.d.recruit.dens.dep.ext.prob = mean(high.seas.d.recruit.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(high.seas.d.recruit.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(high.seas.d.recruit.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(high.seas.d.recruit.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(high.seas.d.recruit.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(high.seas.d.recruit.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(high.seas.d.recruit.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step






## 3.4. Low seasonality case ----
# --------------------------

## 3.4.1. S survival as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - S SURVIVAL")
low.seas.S.surv.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = survS_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.S.surv.dens.dep.sim.lambda = apply(low.seas.S.surv.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.S.surv.dens.dep.stoch.lambda = mean(low.seas.S.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.S.surv.dens.dep.var.lambda = apply(low.seas.S.surv.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.S.surv.dens.dep.ext.prob = mean(low.seas.S.surv.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.S.surv.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.S.surv.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.S.surv.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.S.surv.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.S.surv.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.S.surv.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.4.2. H survival as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - H SURVIVAL")
low.seas.H.surv.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = survH_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.H.surv.dens.dep.sim.lambda = apply(low.seas.H.surv.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.H.surv.dens.dep.stoch.lambda = mean(low.seas.H.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.H.surv.dens.dep.var.lambda = apply(low.seas.H.surv.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.H.surv.dens.dep.ext.prob = mean(low.seas.H.surv.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.H.surv.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.H.surv.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.H.surv.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.H.surv.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.H.surv.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.H.surv.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step



## 3.4.3. D survival as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - D SURVIVAL")
low.seas.D.surv.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = survD_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.D.surv.dens.dep.sim.lambda = apply(low.seas.D.surv.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.D.surv.dens.dep.stoch.lambda = mean(low.seas.D.surv.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.D.surv.dens.dep.var.lambda = apply(low.seas.D.surv.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.D.surv.dens.dep.ext.prob = mean(low.seas.D.surv.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.D.surv.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.D.surv.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.D.surv.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.D.surv.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.D.surv.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.D.surv.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.4.4. Emigration as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - EMIGRATION")
low.seas.emig.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = emig_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.emig.dens.dep.sim.lambda = apply(low.seas.emig.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.emig.dens.dep.stoch.lambda = mean(low.seas.emig.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.emig.dens.dep.var.lambda = apply(low.seas.emig.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.emig.dens.dep.ext.prob = mean(low.seas.emig.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.emig.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.emig.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.emig.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.emig.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.emig.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.emig.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.4.5. Transition as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - TRANSITION")
low.seas.trans.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = transition_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.trans.dens.dep.sim.lambda = apply(low.seas.trans.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.trans.dens.dep.stoch.lambda = mean(low.seas.trans.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.trans.dens.dep.var.lambda = apply(low.seas.trans.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.trans.dens.dep.ext.prob = mean(low.seas.trans.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.trans.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.trans.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.trans.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.trans.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.trans.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.trans.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.4.6. H recruitment as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - H RECRUITMENT")
low.seas.h.recruit.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = recruitment.H_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.h.recruit.dens.dep.sim.lambda = apply(low.seas.h.recruit.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.h.recruit.dens.dep.stoch.lambda = mean(low.seas.h.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.h.recruit.dens.dep.var.lambda = apply(low.seas.h.recruit.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.h.recruit.dens.dep.ext.prob = mean(low.seas.h.recruit.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.h.recruit.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.h.recruit.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.h.recruit.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.h.recruit.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.h.recruit.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.h.recruit.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.4.7. D recruitment as reference ----
# -------------------------------

print("LOW SEASONALITY CASE - D RECRUITMENT")
low.seas.d.recruit.sim.dens.dep = stoch.sim.meerkat(nsimul = 500, nyears = 100, n0 = n0, dens0 = dens0, seasonality = "low", vr = recruitment.D_WithSeason)

# Simulation-wise stochastic lambdas
low.seas.d.recruit.dens.dep.sim.lambda = apply(low.seas.d.recruit.sim.dens.dep[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
low.seas.d.recruit.dens.dep.stoch.lambda = mean(low.seas.d.recruit.dens.dep.sim.lambda, na.rm = T)

# Variance in log lambda across years
low.seas.d.recruit.dens.dep.var.lambda = apply(low.seas.d.recruit.sim.dens.dep[[1]], 1, var, na.rm = T)

# Mean extinction probability
low.seas.d.recruit.dens.dep.ext.prob = mean(low.seas.d.recruit.sim.dens.dep[[2]])

# Yearly log lambdas for each simulation, and mean
matplot(t(low.seas.d.recruit.sim.dens.dep[[1]]), type = "l", col = rainbow(nrow(low.seas.d.recruit.sim.dens.dep[[1]])), main = "Simulations yearly log lambda - Control scenario - Density dependence", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(low.seas.d.recruit.sim.dens.dep[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(low.seas.d.recruit.sim.dens.dep[[3]]), type = "l", col = rainbow(nrow(low.seas.d.recruit.sim.dens.dep[[3]])), main = "Density over seasonal time steps - Control scenario - Density dependence", xlab = "Time steps (seasons)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(low.seas.d.recruit.sim.dens.dep[[3]], 2, mean), lwd = 2) # Average density accross simulations for each time step




###########################################################################
#
# 4. Saving files and data ----
#
###########################################################################

save.image("MeerkatsProjectionsResults_DensityOnly.RData")
