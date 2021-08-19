############################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script uses the seasonal vital rate-estimates for the marmots population. 
# The aim of this script is to calculate (1) the elasticity of the population growth rate to changes in the mean or standard deviation of each vital rate, and (2) the relative effect of variability (Morris et al. 2008), i.e. the proportion of the stochastic elasticity E(S) attributed to changes in the variability of a given vital rate category. 
#
# Author: Eva Conquet
#
# Morris, W. F., et al. 2008. Longevity can buffer plant and animal populations 
# against changing climatic variability. Ecology 89: 19–25.
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
  library(ggplot2)
}

load.librairies()


## 1.3. Loading and preparing data ----
# --------------------------------

source("02_Marmots_MPMs.R")
marmots.vr = read.csv("Marmots_VitalRatesEstimates.csv")




###########################################################################
#
# 2. Preparing the simulations ----
#
###########################################################################

## 2.1. Getting mean vital rates ----
# ------------------------------

mean.surv.j.winter    = mean(marmots.vr$survJ)
mean.surv.y.winter    = mean(marmots.vr$survY_winter)
mean.surv.y.summer    = mean(marmots.vr$survY_summer)
mean.surv.nr.winter   = mean(marmots.vr$survNR_winter)
mean.surv.nr.summer   = mean(marmots.vr$survNR_summer)
mean.surv.r.winter    = mean(marmots.vr$survR_winter)
mean.surv.r.summer    = mean(marmots.vr$survR_summer)
mean.trans.yr.winter  = mean(marmots.vr$transYR)
mean.trans.nrr.winter = mean(marmots.vr$transNRR)
mean.trans.rr.winter  = mean(marmots.vr$transRR)
mean.recruit.summer   = mean(marmots.vr$recruit)


## 2.2. Simulations parameters and initial values ----
# -----------------------------------------------

# Define population vector
n0 = c(12, 11, 9, 7)


# 100 simulations of 1000 years with a cutoff at 500 years to look at asymptotic dynamics
n.sim = 100
sim.n.years.total = 1000
sim.n.years.discard = 500
sim.n.years.asymptotic = sim.n.years.total - sim.n.years.discard
t2 <- sim.n.years.asymptotic + 2 * sim.n.years.discard # Time vector for the backward iteration to calculate the left eigenvectors


# Elasticities storage

elasticity.mean.j.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.y.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.y.surv.summer    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.nr.surv.winter   = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.nr.surv.summer   = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.r.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.r.surv.summer    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.yr.trans.winter  = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.nrr.trans.winter = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.rr.trans.winter  = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.recruit.summer   = array(0, c(length(n0), length(n0), n.sim))

elasticity.sd.j.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.y.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.y.surv.summer    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.nr.surv.winter   = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.nr.surv.summer   = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.r.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.r.surv.summer    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.yr.trans.winter  = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.nrr.trans.winter = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.rr.trans.winter  = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.recruit.summer   = array(0, c(length(n0), length(n0), n.sim))

elasticity.stochastic.j.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.y.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.y.surv.summer    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.nr.surv.winter   = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.nr.surv.summer   = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.r.surv.winter    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.r.surv.summer    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.yr.trans.winter  = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.nrr.trans.winter = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.rr.trans.winter  = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.recruit.summer   = array(0, c(length(n0), length(n0), n.sim))




###########################################################################
#
# 3. Simulating dynamics to calculate elasticities -----
#
###########################################################################

for(x in 1:n.sim){ 
  
  # Initialize list of C matrices (with difference between perturbed and unperturbed matrices)
  
  mean.C.mat.juv.surv = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.yearling.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.yearling.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.nr.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.nr.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.r.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.r.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.yr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.nrr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.rr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.recruitment = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  
  sd.C.mat.juv.surv = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.yearling.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.yearling.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.nr.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.nr.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.r.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.r.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.yr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.nrr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.rr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.recruitment = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  
  vr.value.C.mat.juv.surv = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.yearling.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.yearling.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.nr.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.nr.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.r.surv.summer = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.r.surv.winter = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.yr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.nrr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.rr.trans = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.recruitment = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  
  # Initialize population vectors (n0)
  
  vec1 = n0
  vec1 <- vec1 / sum(vec1)
  vec1 <- t(vec1) 
  vec2 <- vec1
  
  uvecs <- array(0, c(length(n0), sim.n.years.asymptotic)) # Right eigenvectors
  vvecs <- array(0, c(length(n0), sim.n.years.asymptotic)) # Left eigenvectors
  growth <- array(0, sim.n.years.asymptotic)
  
  # Simulate environmental states, i.e., years over trun simulations
  
  year.sim = sample(marmots.vr$year, t2 + 1, replace = T)
  
  for(i in 1:sim.n.years.total){
    
    year = year.sim[i]
    
    # Create annual matrix
    annual.mat = marmots.buildmat.summer(marmots.vr, year) %*% marmots.buildmat.winter(marmots.vr, year)
    
    # Create a perturbed data frame for each vital rate for the current year and an annual matrix for each perturbed vital rate, for a perturbation with the mean vital rate, with the standard deviation, and the value of the vital rate.
    
      # Juvenile survival
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.juv.surv = marmots.vr
    mean.perturbed.marmots.vr.juv.surv$survJ[which(mean.perturbed.marmots.vr.juv.surv$year == year)] = mean.perturbed.marmots.vr.juv.surv$survJ[which(mean.perturbed.marmots.vr.juv.surv$year == year)] + mean.surv.j.winter
    
    mean.perturbed.annual.mat.juv.surv = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.juv.surv, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.juv.surv, year = year)
    
    mean.C.mat.juv.surv[, , i] = mean.perturbed.annual.mat.juv.surv - annual.mat # C matrix with perturbed elements
     
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.juv.surv = marmots.vr
    sd.perturbed.marmots.vr.juv.surv$survJ[which(sd.perturbed.marmots.vr.juv.surv$year == year)] = sd.perturbed.marmots.vr.juv.surv$survJ[which(sd.perturbed.marmots.vr.juv.surv$year == year)] + (sd.perturbed.marmots.vr.juv.surv$survJ[which(sd.perturbed.marmots.vr.juv.surv$year == year)] - mean.surv.j.winter)
    
    sd.perturbed.annual.mat.juv.surv = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.juv.surv, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.juv.surv, year = year)
    
    sd.C.mat.juv.surv[, , i] = sd.perturbed.annual.mat.juv.surv - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.juv.surv = marmots.vr
    vr.value.perturbed.marmots.vr.juv.surv$survJ[which(vr.value.perturbed.marmots.vr.juv.surv$year == year)] = vr.value.perturbed.marmots.vr.juv.surv$survJ[which(vr.value.perturbed.marmots.vr.juv.surv$year == year)] + vr.value.perturbed.marmots.vr.juv.surv$survJ[which(vr.value.perturbed.marmots.vr.juv.surv$year == year)]
    
    vr.value.perturbed.annual.mat.juv.surv = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.juv.surv, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.juv.surv, year = year)
    
    vr.value.C.mat.juv.surv[, , i] = vr.value.perturbed.annual.mat.juv.surv - annual.mat # C matrix with perturbed elements    
    
    ################################
    
      # Yearling winter survival
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.yearling.surv.winter = marmots.vr
    mean.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(mean.perturbed.marmots.vr.yearling.surv.winter$year == year)] = mean.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(mean.perturbed.marmots.vr.yearling.surv.winter$year == year)] + mean.surv.y.winter
    
    mean.perturbed.annual.mat.yearling.surv.winter = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.yearling.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.yearling.surv.winter, year = year)
    
    mean.C.mat.yearling.surv.winter[, , i] = mean.perturbed.annual.mat.yearling.surv.winter - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.yearling.surv.winter = marmots.vr
    sd.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(sd.perturbed.marmots.vr.yearling.surv.winter$year == year)] = sd.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(sd.perturbed.marmots.vr.yearling.surv.winter$year == year)] + (sd.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(sd.perturbed.marmots.vr.yearling.surv.winter$year == year)] - mean.surv.y.winter)
    
    sd.perturbed.annual.mat.yearling.surv.winter = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.yearling.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.yearling.surv.winter, year = year)
    
    sd.C.mat.yearling.surv.winter[, , i] = sd.perturbed.annual.mat.yearling.surv.winter - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.yearling.surv.winter = marmots.vr
    vr.value.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(vr.value.perturbed.marmots.vr.yearling.surv.winter$year == year)] = vr.value.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(vr.value.perturbed.marmots.vr.yearling.surv.winter$year == year)] + vr.value.perturbed.marmots.vr.yearling.surv.winter$survY_winter[which(vr.value.perturbed.marmots.vr.yearling.surv.winter$year == year)]
    
    vr.value.perturbed.annual.mat.yearling.surv.winter = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.yearling.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.yearling.surv.winter, year = year)
    
    vr.value.C.mat.yearling.surv.winter[, , i] = vr.value.perturbed.annual.mat.yearling.surv.winter - annual.mat # C matrix with perturbed elements 
    
    ################################  
    
      # Yearling summer survival
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.yearling.surv.summer = marmots.vr
    mean.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(mean.perturbed.marmots.vr.yearling.surv.summer$year == year)] = mean.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(mean.perturbed.marmots.vr.yearling.surv.summer$year == year)] + mean.surv.y.summer
    
    mean.perturbed.annual.mat.yearling.surv.summer = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.yearling.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.yearling.surv.summer, year = year)
    
    mean.C.mat.yearling.surv.summer[, , i] = mean.perturbed.annual.mat.yearling.surv.summer - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.yearling.surv.summer = marmots.vr
    sd.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(sd.perturbed.marmots.vr.yearling.surv.summer$year == year)] = sd.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(sd.perturbed.marmots.vr.yearling.surv.summer$year == year)] + (sd.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(sd.perturbed.marmots.vr.yearling.surv.summer$year == year)] - mean.surv.y.summer)
    
    sd.perturbed.annual.mat.yearling.surv.summer = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.yearling.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.yearling.surv.summer, year = year)
    
    sd.C.mat.yearling.surv.summer[, , i] = sd.perturbed.annual.mat.yearling.surv.summer - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.yearling.surv.summer = marmots.vr
    vr.value.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(vr.value.perturbed.marmots.vr.yearling.surv.summer$year == year)] = vr.value.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(vr.value.perturbed.marmots.vr.yearling.surv.summer$year == year)] + vr.value.perturbed.marmots.vr.yearling.surv.summer$survY_summer[which(vr.value.perturbed.marmots.vr.yearling.surv.summer$year == year)]
    
    vr.value.perturbed.annual.mat.yearling.surv.summer = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.yearling.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.yearling.surv.summer, year = year)
    
    vr.value.C.mat.yearling.surv.summer[, , i] = vr.value.perturbed.annual.mat.yearling.surv.summer - annual.mat # C matrix with perturbed elements 
    
    
    ################################
    
      # Non reproductive winter survival
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.nr.surv.winter = marmots.vr
    mean.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(mean.perturbed.marmots.vr.nr.surv.winter$year == year)] = mean.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(mean.perturbed.marmots.vr.nr.surv.winter$year == year)] + mean.surv.nr.winter
    
    mean.perturbed.annual.mat.nr.surv.winter = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.nr.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.nr.surv.winter, year = year)
    
    mean.C.mat.nr.surv.winter[, , i] = mean.perturbed.annual.mat.nr.surv.winter - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.nr.surv.winter = marmots.vr
    sd.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(sd.perturbed.marmots.vr.nr.surv.winter$year == year)] = sd.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(sd.perturbed.marmots.vr.nr.surv.winter$year == year)] + (sd.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(sd.perturbed.marmots.vr.nr.surv.winter$year == year)] - mean.surv.nr.winter)
    
    sd.perturbed.annual.mat.nr.surv.winter = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.nr.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.nr.surv.winter, year = year)
    
    sd.C.mat.nr.surv.winter[, , i] = sd.perturbed.annual.mat.nr.surv.winter - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.nr.surv.winter = marmots.vr
    vr.value.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(vr.value.perturbed.marmots.vr.nr.surv.winter$year == year)] = vr.value.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(vr.value.perturbed.marmots.vr.nr.surv.winter$year == year)] + vr.value.perturbed.marmots.vr.nr.surv.winter$survNR_winter[which(vr.value.perturbed.marmots.vr.nr.surv.winter$year == year)]
    
    vr.value.perturbed.annual.mat.nr.surv.winter = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.nr.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.nr.surv.winter, year = year)
    
    vr.value.C.mat.nr.surv.winter[, , i] = vr.value.perturbed.annual.mat.nr.surv.winter - annual.mat # C matrix with perturbed elements 
    
    ################################
    
      # Non reproductive summer survival
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.nr.surv.summer = marmots.vr
    mean.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(mean.perturbed.marmots.vr.nr.surv.summer$year == year)] = mean.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(mean.perturbed.marmots.vr.nr.surv.summer$year == year)] + mean.surv.nr.summer
    
    mean.perturbed.annual.mat.nr.surv.summer = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.nr.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.nr.surv.summer, year = year)
    
    mean.C.mat.nr.surv.summer[, , i] = mean.perturbed.annual.mat.nr.surv.summer - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.nr.surv.summer = marmots.vr
    sd.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(sd.perturbed.marmots.vr.nr.surv.summer$year == year)] = sd.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(sd.perturbed.marmots.vr.nr.surv.summer$year == year)] + (sd.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(sd.perturbed.marmots.vr.nr.surv.summer$year == year)] - mean.surv.nr.summer)
    
    sd.perturbed.annual.mat.nr.surv.summer = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.nr.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.nr.surv.summer, year = year)
    
    sd.C.mat.nr.surv.summer[, , i] = sd.perturbed.annual.mat.nr.surv.summer - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.nr.surv.summer = marmots.vr
    vr.value.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(vr.value.perturbed.marmots.vr.nr.surv.summer$year == year)] = vr.value.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(vr.value.perturbed.marmots.vr.nr.surv.summer$year == year)] + vr.value.perturbed.marmots.vr.nr.surv.summer$survNR_summer[which(vr.value.perturbed.marmots.vr.nr.surv.summer$year == year)]
    
    vr.value.perturbed.annual.mat.nr.surv.summer = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.nr.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.nr.surv.summer, year = year)
    
    vr.value.C.mat.nr.surv.summer[, , i] = vr.value.perturbed.annual.mat.nr.surv.summer - annual.mat # C matrix with perturbed elements
    
    ################################
    
      # Reproductive winter survival

    # Perturbation with the mean
    mean.perturbed.marmots.vr.r.surv.winter = marmots.vr
    mean.perturbed.marmots.vr.r.surv.winter$survR_winter[which(mean.perturbed.marmots.vr.r.surv.winter$year == year)] = mean.perturbed.marmots.vr.r.surv.winter$survR_winter[which(mean.perturbed.marmots.vr.r.surv.winter$year == year)] + mean.surv.r.winter
    
    mean.perturbed.annual.mat.r.surv.winter = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.r.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.r.surv.winter, year = year)
    
    mean.C.mat.r.surv.winter[, , i] = mean.perturbed.annual.mat.r.surv.winter - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.r.surv.winter = marmots.vr
    sd.perturbed.marmots.vr.r.surv.winter$survR_winter[which(sd.perturbed.marmots.vr.r.surv.winter$year == year)] = sd.perturbed.marmots.vr.r.surv.winter$survR_winter[which(sd.perturbed.marmots.vr.r.surv.winter$year == year)] + (sd.perturbed.marmots.vr.r.surv.winter$survR_winter[which(sd.perturbed.marmots.vr.r.surv.winter$year == year)] - mean.surv.r.winter)
    
    sd.perturbed.annual.mat.r.surv.winter = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.r.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.r.surv.winter, year = year)
    
    sd.C.mat.r.surv.winter[, , i] = sd.perturbed.annual.mat.r.surv.winter - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.r.surv.winter = marmots.vr
    vr.value.perturbed.marmots.vr.r.surv.winter$survR_winter[which(vr.value.perturbed.marmots.vr.r.surv.winter$year == year)] = vr.value.perturbed.marmots.vr.r.surv.winter$survR_winter[which(vr.value.perturbed.marmots.vr.r.surv.winter$year == year)] + vr.value.perturbed.marmots.vr.r.surv.winter$survR_winter[which(vr.value.perturbed.marmots.vr.r.surv.winter$year == year)]
    
    vr.value.perturbed.annual.mat.r.surv.winter = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.r.surv.winter, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.r.surv.winter, year = year)
    
    vr.value.C.mat.r.surv.winter[, , i] = vr.value.perturbed.annual.mat.r.surv.winter - annual.mat # C matrix with perturbed elements 
    
    
    ################################
    
      # Reproductive summer survival 
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.r.surv.summer = marmots.vr
    mean.perturbed.marmots.vr.r.surv.summer$survR_summer[which(mean.perturbed.marmots.vr.r.surv.summer$year == year)] = mean.perturbed.marmots.vr.r.surv.summer$survR_summer[which(mean.perturbed.marmots.vr.r.surv.summer$year == year)] + mean.surv.r.summer
    
    mean.perturbed.annual.mat.r.surv.summer = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.r.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.r.surv.summer, year = year)
    
    mean.C.mat.r.surv.summer[, , i] = mean.perturbed.annual.mat.r.surv.summer - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.r.surv.summer = marmots.vr
    sd.perturbed.marmots.vr.r.surv.summer$survR_summer[which(sd.perturbed.marmots.vr.r.surv.summer$year == year)] = sd.perturbed.marmots.vr.r.surv.summer$survR_summer[which(sd.perturbed.marmots.vr.r.surv.summer$year == year)] + (sd.perturbed.marmots.vr.r.surv.summer$survR_summer[which(sd.perturbed.marmots.vr.r.surv.summer$year == year)] - mean.surv.r.summer)
    
    sd.perturbed.annual.mat.r.surv.summer = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.r.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.r.surv.summer, year = year)
    
    sd.C.mat.r.surv.summer[, , i] = sd.perturbed.annual.mat.r.surv.summer - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.r.surv.summer = marmots.vr
    vr.value.perturbed.marmots.vr.r.surv.summer$survR_summer[which(vr.value.perturbed.marmots.vr.r.surv.summer$year == year)] = vr.value.perturbed.marmots.vr.r.surv.summer$survR_summer[which(vr.value.perturbed.marmots.vr.r.surv.summer$year == year)] + vr.value.perturbed.marmots.vr.r.surv.summer$survR_summer[which(vr.value.perturbed.marmots.vr.r.surv.summer$year == year)]
    
    vr.value.perturbed.annual.mat.r.surv.summer = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.r.surv.summer, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.r.surv.summer, year = year)
    
    vr.value.C.mat.r.surv.summer[, , i] = vr.value.perturbed.annual.mat.r.surv.summer - annual.mat # C matrix with perturbed elements
    
    ################################
    
      # Yearling-Reproductive transition

    # Perturbation with the mean
    mean.perturbed.marmots.vr.yr.trans = marmots.vr
    mean.perturbed.marmots.vr.yr.trans$transYR[which(mean.perturbed.marmots.vr.yr.trans$year == year)] = mean.perturbed.marmots.vr.yr.trans$transYR[which(mean.perturbed.marmots.vr.yr.trans$year == year)] + mean.trans.yr.winter
    
    mean.perturbed.annual.mat.yr.trans = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.yr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.yr.trans, year = year)
    
    mean.C.mat.yr.trans[, , i] = mean.perturbed.annual.mat.yr.trans - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.yr.trans = marmots.vr
    sd.perturbed.marmots.vr.yr.trans$transYR[which(sd.perturbed.marmots.vr.yr.trans$year == year)] = sd.perturbed.marmots.vr.yr.trans$transYR[which(sd.perturbed.marmots.vr.yr.trans$year == year)] + (sd.perturbed.marmots.vr.yr.trans$transYR[which(sd.perturbed.marmots.vr.yr.trans$year == year)] - mean.trans.yr.winter)
    
    sd.perturbed.annual.mat.yr.trans = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.yr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.yr.trans, year = year)
    
    sd.C.mat.yr.trans[, , i] = sd.perturbed.annual.mat.yr.trans - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.yr.trans = marmots.vr
    vr.value.perturbed.marmots.vr.yr.trans$transYR[which(vr.value.perturbed.marmots.vr.yr.trans$year == year)] = vr.value.perturbed.marmots.vr.yr.trans$transYR[which(vr.value.perturbed.marmots.vr.yr.trans$year == year)] + vr.value.perturbed.marmots.vr.yr.trans$transYR[which(vr.value.perturbed.marmots.vr.yr.trans$year == year)]
    
    vr.value.perturbed.annual.mat.yr.trans = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.yr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.yr.trans, year = year)
    
    vr.value.C.mat.yr.trans[, , i] = vr.value.perturbed.annual.mat.yr.trans - annual.mat # C matrix with perturbed elements
     
    ################################
    
      # Non reproductive-Reproductive transition
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.nrr.trans = marmots.vr
    mean.perturbed.marmots.vr.nrr.trans$transNRR[which(mean.perturbed.marmots.vr.nrr.trans$year == year)] = mean.perturbed.marmots.vr.nrr.trans$transNRR[which(mean.perturbed.marmots.vr.nrr.trans$year == year)] + mean.trans.nrr.winter
    
    mean.perturbed.annual.mat.nrr.trans = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.nrr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.nrr.trans, year = year)
    
    mean.C.mat.nrr.trans[, , i] = mean.perturbed.annual.mat.nrr.trans - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.nrr.trans = marmots.vr
    sd.perturbed.marmots.vr.nrr.trans$transNRR[which(sd.perturbed.marmots.vr.nrr.trans$year == year)] = sd.perturbed.marmots.vr.nrr.trans$transNRR[which(sd.perturbed.marmots.vr.nrr.trans$year == year)] + (sd.perturbed.marmots.vr.nrr.trans$transNRR[which(sd.perturbed.marmots.vr.nrr.trans$year == year)] - mean.trans.nrr.winter)
    
    sd.perturbed.annual.mat.nrr.trans = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.nrr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.nrr.trans, year = year)
    
    sd.C.mat.nrr.trans[, , i] = sd.perturbed.annual.mat.nrr.trans - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.nrr.trans = marmots.vr
    vr.value.perturbed.marmots.vr.nrr.trans$transNRR[which(vr.value.perturbed.marmots.vr.nrr.trans$year == year)] = vr.value.perturbed.marmots.vr.nrr.trans$transNRR[which(vr.value.perturbed.marmots.vr.nrr.trans$year == year)] + vr.value.perturbed.marmots.vr.nrr.trans$transNRR[which(vr.value.perturbed.marmots.vr.nrr.trans$year == year)]
    
    vr.value.perturbed.annual.mat.nrr.trans = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.nrr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.nrr.trans, year = year)
    
    vr.value.C.mat.nrr.trans[, , i] = vr.value.perturbed.annual.mat.nrr.trans - annual.mat # C matrix with perturbed elements
    
    ################################
    
      # Reproductive-Reproductive stasis
    
    # Perturbation with the mean
    mean.perturbed.marmots.vr.rr.trans = marmots.vr
    mean.perturbed.marmots.vr.rr.trans$transRR[which(mean.perturbed.marmots.vr.rr.trans$year == year)] = mean.perturbed.marmots.vr.rr.trans$transRR[which(mean.perturbed.marmots.vr.rr.trans$year == year)] + mean.trans.rr.winter
    
    mean.perturbed.annual.mat.rr.trans = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.rr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.rr.trans, year = year)
    
    mean.C.mat.rr.trans[, , i] = mean.perturbed.annual.mat.rr.trans - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.rr.trans = marmots.vr
    sd.perturbed.marmots.vr.rr.trans$transRR[which(sd.perturbed.marmots.vr.rr.trans$year == year)] = sd.perturbed.marmots.vr.rr.trans$transRR[which(sd.perturbed.marmots.vr.rr.trans$year == year)] + (sd.perturbed.marmots.vr.rr.trans$transRR[which(sd.perturbed.marmots.vr.rr.trans$year == year)] - mean.trans.rr.winter)
    
    sd.perturbed.annual.mat.rr.trans = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.rr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.rr.trans, year = year)
    
    sd.C.mat.rr.trans[, , i] = sd.perturbed.annual.mat.rr.trans - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.rr.trans = marmots.vr
    vr.value.perturbed.marmots.vr.rr.trans$transRR[which(vr.value.perturbed.marmots.vr.rr.trans$year == year)] = vr.value.perturbed.marmots.vr.rr.trans$transRR[which(vr.value.perturbed.marmots.vr.rr.trans$year == year)] + vr.value.perturbed.marmots.vr.rr.trans$transRR[which(vr.value.perturbed.marmots.vr.rr.trans$year == year)]
    
    vr.value.perturbed.annual.mat.rr.trans = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.rr.trans, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.rr.trans, year = year)
    
    vr.value.C.mat.rr.trans[, , i] = vr.value.perturbed.annual.mat.rr.trans - annual.mat # C matrix with perturbed elements
    
    ################################
    
      # Recruitment
  
    # Perturbation with the mean
    mean.perturbed.marmots.vr.recruitment = marmots.vr
    mean.perturbed.marmots.vr.recruitment$recruit[which(mean.perturbed.marmots.vr.recruitment$year == year)] = mean.perturbed.marmots.vr.recruitment$recruit[which(mean.perturbed.marmots.vr.recruitment$year == year)] + mean.recruit.summer
    
    mean.perturbed.annual.mat.recruitment = marmots.buildmat.summer(marmots.vr = mean.perturbed.marmots.vr.recruitment, year = year) %*% marmots.buildmat.winter(marmots.vr = mean.perturbed.marmots.vr.recruitment, year = year)
    
    mean.C.mat.recruitment[, , i] = mean.perturbed.annual.mat.recruitment - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.marmots.vr.recruitment = marmots.vr
    sd.perturbed.marmots.vr.recruitment$recruit[which(sd.perturbed.marmots.vr.recruitment$year == year)] = sd.perturbed.marmots.vr.recruitment$recruit[which(sd.perturbed.marmots.vr.recruitment$year == year)] + (sd.perturbed.marmots.vr.recruitment$recruit[which(sd.perturbed.marmots.vr.recruitment$year == year)] - mean.recruit.summer)
    
    sd.perturbed.annual.mat.recruitment = marmots.buildmat.summer(marmots.vr = sd.perturbed.marmots.vr.recruitment, year = year) %*% marmots.buildmat.winter(marmots.vr = sd.perturbed.marmots.vr.recruitment, year = year)
    
    sd.C.mat.recruitment[, , i] = sd.perturbed.annual.mat.recruitment - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.marmots.vr.recruitment = marmots.vr
    vr.value.perturbed.marmots.vr.recruitment$recruit[which(vr.value.perturbed.marmots.vr.recruitment$year == year)] = vr.value.perturbed.marmots.vr.recruitment$recruit[which(vr.value.perturbed.marmots.vr.recruitment$year == year)] + vr.value.perturbed.marmots.vr.recruitment$recruit[which(vr.value.perturbed.marmots.vr.recruitment$year == year)]
    
    vr.value.perturbed.annual.mat.recruitment = marmots.buildmat.summer(marmots.vr = vr.value.perturbed.marmots.vr.recruitment, year = year) %*% marmots.buildmat.winter(marmots.vr = vr.value.perturbed.marmots.vr.recruitment, year = year)
    
    vr.value.C.mat.recruitment[, , i] = vr.value.perturbed.annual.mat.recruitment - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Forward and backward iterations to get v and w vectors 
    
    vec1 <- annual.mat %*% as.numeric(vec1)
    growth1 <- sum(vec1) # population growth at one time step 
    vec1 <- vec1 / growth1
    
      # Forward iteration
    if(i > sim.n.years.discard){      
      i1 <- i - sim.n.years.discard
      uvecs[, i1] <- vec1
      growth[i1] <- growth1
    }
    
      # Backward iteration
    j <- (t2 - i + 1)
    
    year = year.sim[j]
    
    annual.mat = marmots.buildmat.summer(marmots.vr, year) %*% marmots.buildmat.winter(marmots.vr, year) # Annual matrix

    vec2 <- t(as.numeric(vec2) %*% annual.mat)
    vec2 <- vec2 / (sum(vec2))
    
    if (i > sim.n.years.discard){     
      vvecs[, j - sim.n.years.discard] <- vec2
    }
  }
  
  
  for (z in (sim.n.years.discard + 1):(sim.n.years.discard + sim.n.years.asymptotic - 1)){ # Start after a time lag (sim.n.years.discard)
    
    itime <- z + 1 - sim.n.years.discard
    i1 <- z - sim.n.years.discard
    
    # Calculating stochastic elasticities (Translation Tuljapurkar et al. 2003)
    
    # Scalar eigenvector product (denominator)
    scale1 <- (t(vvecs[, itime])) %*% (uvecs[, itime])
    a = as.numeric((1 / (scale1 * growth[itime])))
    
    # Elasticities juvenile survival
    mean.pert.juv.surv <- a * (diag(vvecs[, itime]) %*% mean.C.mat.juv.surv[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.juv.surv <- a * (diag(vvecs[, itime]) %*% sd.C.mat.juv.surv[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.juv.surv <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.juv.surv[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.j.surv.winter[, , x] <- elasticity.mean.j.surv.winter[, , x] + mean.pert.juv.surv
    elasticity.sd.j.surv.winter[, , x] <- elasticity.sd.j.surv.winter[, , x] + sd.pert.juv.surv
    elasticity.stochastic.j.surv.winter[, , x] <- elasticity.stochastic.j.surv.winter[, , x] + vr.value.pert.juv.surv
    
    # Elasticities yearling winter survival
    mean.pert.yearling.surv.winter <- a * (diag(vvecs[, itime]) %*% mean.C.mat.yearling.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.yearling.surv.winter <- a * (diag(vvecs[, itime]) %*% sd.C.mat.yearling.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.yearling.surv.winter <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.yearling.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.y.surv.winter[, , x] <- elasticity.mean.y.surv.winter[, , x] + mean.pert.yearling.surv.winter
    elasticity.sd.y.surv.winter[, , x] <- elasticity.sd.y.surv.winter[, , x] + sd.pert.yearling.surv.winter
    elasticity.stochastic.y.surv.winter[, , x] <- elasticity.stochastic.y.surv.winter[, , x] + vr.value.pert.yearling.surv.winter
    
    # Elasticities yearling summer survival
    mean.pert.yearling.surv.summer <- a * (diag(vvecs[, itime]) %*% mean.C.mat.yearling.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.yearling.surv.summer <- a * (diag(vvecs[, itime]) %*% sd.C.mat.yearling.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.yearling.surv.summer <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.yearling.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.y.surv.summer[, , x] <- elasticity.mean.y.surv.summer[,,x] + mean.pert.yearling.surv.summer
    elasticity.sd.y.surv.summer[, , x] <- elasticity.sd.y.surv.summer[,,x] + sd.pert.yearling.surv.summer
    elasticity.stochastic.y.surv.summer[, , x] <- elasticity.stochastic.y.surv.summer[,,x] + vr.value.pert.yearling.surv.summer
    
    # Elasticities non reproductive winter survival
    mean.pert.nr.surv.winter <- a * (diag(vvecs[, itime]) %*% mean.C.mat.nr.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.nr.surv.winter <- a * (diag(vvecs[, itime]) %*% sd.C.mat.nr.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.nr.surv.winter <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.nr.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.nr.surv.winter[, , x] <- elasticity.mean.nr.surv.winter[, , x] + mean.pert.nr.surv.winter
    elasticity.sd.nr.surv.winter[, , x] <- elasticity.sd.nr.surv.winter[, , x] + sd.pert.nr.surv.winter
    elasticity.stochastic.nr.surv.winter[, , x] <- elasticity.stochastic.nr.surv.winter[, , x] + vr.value.pert.nr.surv.winter
    
    # Elasticities non reproductive summer survival
    mean.pert.nr.surv.summer <- a * (diag(vvecs[, itime]) %*% mean.C.mat.nr.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.nr.surv.summer <- a * (diag(vvecs[, itime]) %*% sd.C.mat.nr.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.nr.surv.summer <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.nr.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.nr.surv.summer[, , x] <- elasticity.mean.nr.surv.summer[, , x] + mean.pert.nr.surv.summer
    elasticity.sd.nr.surv.summer[, , x] <- elasticity.sd.nr.surv.summer[, , x] + sd.pert.nr.surv.summer
    elasticity.stochastic.nr.surv.summer[, , x] <- elasticity.stochastic.nr.surv.summer[, , x] + vr.value.pert.nr.surv.summer
    
    # Elasticities reproductive winter survival
    mean.pert.r.surv.winter <- a * (diag(vvecs[, itime]) %*% mean.C.mat.r.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.r.surv.winter <- a * (diag(vvecs[, itime]) %*% sd.C.mat.r.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.r.surv.winter <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.r.surv.winter[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.r.surv.winter[, , x] <- elasticity.mean.r.surv.winter[, , x] + mean.pert.r.surv.winter
    elasticity.sd.r.surv.winter[, , x] <- elasticity.sd.r.surv.winter[, , x] + sd.pert.r.surv.winter
    elasticity.stochastic.r.surv.winter[, , x] <- elasticity.stochastic.r.surv.winter[, , x] + vr.value.pert.r.surv.winter
    
    # Elasticities reproductive summer survival
    mean.pert.r.surv.summer <- a * (diag(vvecs[, itime]) %*% mean.C.mat.r.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.r.surv.summer <- a * (diag(vvecs[, itime]) %*% sd.C.mat.r.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.r.surv.summer <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.r.surv.summer[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.r.surv.summer[, , x] <- elasticity.mean.r.surv.summer[, , x] + mean.pert.r.surv.summer
    elasticity.sd.r.surv.summer[, , x] <- elasticity.sd.r.surv.summer[, , x] + sd.pert.r.surv.summer
    elasticity.stochastic.r.surv.summer[, , x] <- elasticity.stochastic.r.surv.summer[, , x] + vr.value.pert.r.surv.summer
    
    # Elasticities yearling-reproductive transition
    mean.pert.yr.trans <- a * (diag(vvecs[, itime]) %*% mean.C.mat.yr.trans[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.yr.trans <- a * (diag(vvecs[, itime]) %*% sd.C.mat.yr.trans[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.yr.trans <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.yr.trans[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.yr.trans.winter[, , x] <- elasticity.mean.yr.trans.winter[, , x] + mean.pert.yr.trans
    elasticity.sd.yr.trans.winter[, , x] <- elasticity.sd.yr.trans.winter[, , x] + sd.pert.yr.trans
    elasticity.stochastic.yr.trans.winter[, , x] <- elasticity.stochastic.yr.trans.winter[, , x] + vr.value.pert.yr.trans
    
    # Elasticities non reproductive-reproductive transition
    mean.pert.nrr.trans <- a * (diag(vvecs[, itime]) %*% mean.C.mat.nrr.trans[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.nrr.trans <- a * (diag(vvecs[, itime]) %*% sd.C.mat.nrr.trans[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.nrr.trans <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.nrr.trans[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.nrr.trans.winter[, , x] <- elasticity.mean.nrr.trans.winter[, , x] + mean.pert.nrr.trans
    elasticity.sd.nrr.trans.winter[, , x] <- elasticity.sd.nrr.trans.winter[, , x] + sd.pert.nrr.trans
    elasticity.stochastic.nrr.trans.winter[, , x] <- elasticity.stochastic.nrr.trans.winter[, , x] + vr.value.pert.nrr.trans
    
    # Elasticities reproductive-reproductive transition
    mean.pert.rr.trans <- a * (diag(vvecs[, itime]) %*% mean.C.mat.rr.trans[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.rr.trans <- a * (diag(vvecs[, itime]) %*% sd.C.mat.rr.trans[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.rr.trans <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.rr.trans[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.rr.trans.winter[, , x] <- elasticity.mean.rr.trans.winter[, , x] + mean.pert.rr.trans
    elasticity.sd.rr.trans.winter[, , x] <- elasticity.sd.rr.trans.winter[, , x] + sd.pert.rr.trans
    elasticity.stochastic.rr.trans.winter[, , x] <- elasticity.stochastic.rr.trans.winter[, , x] + vr.value.pert.rr.trans
    
  }
  
  elasticity.mean.j.surv.winter[, , x] = elasticity.mean.j.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.j.surv.winter[, , x] = elasticity.sd.j.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.j.surv.winter[, , x] = elasticity.stochastic.j.surv.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.y.surv.winter[, , x] = elasticity.mean.y.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.y.surv.winter[, , x] = elasticity.sd.y.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.y.surv.winter[, , x] = elasticity.stochastic.y.surv.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.y.surv.summer[, , x] = elasticity.mean.y.surv.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.y.surv.summer[, , x] = elasticity.sd.y.surv.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.y.surv.summer[, , x] = elasticity.stochastic.y.surv.summer[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.nr.surv.winter[, , x] = elasticity.mean.nr.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.nr.surv.winter[, , x] = elasticity.sd.nr.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.nr.surv.winter[, , x] = elasticity.stochastic.nr.surv.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.nr.surv.summer[, , x] = elasticity.mean.nr.surv.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.nr.surv.summer[, , x] = elasticity.sd.nr.surv.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.nr.surv.summer[, , x] = elasticity.stochastic.nr.surv.summer[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.r.surv.winter[, , x] = elasticity.mean.r.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.r.surv.winter[, , x] = elasticity.sd.r.surv.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.r.surv.winter[, , x] = elasticity.stochastic.r.surv.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.r.surv.summer[, , x] = elasticity.mean.r.surv.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.r.surv.summer[, , x] = elasticity.sd.r.surv.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.r.surv.summer[, , x] = elasticity.stochastic.r.surv.summer[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.yr.trans.winter[, , x] = elasticity.mean.yr.trans.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.yr.trans.winter[, , x] = elasticity.sd.yr.trans.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.yr.trans.winter[, , x] = elasticity.stochastic.yr.trans.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.nrr.trans.winter[, , x] = elasticity.mean.nrr.trans.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.nrr.trans.winter[, , x] = elasticity.sd.nrr.trans.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.nrr.trans.winter[, , x] = elasticity.stochastic.nrr.trans.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.rr.trans.winter[, , x] = elasticity.mean.rr.trans.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.rr.trans.winter[, , x] = elasticity.sd.rr.trans.winter[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.rr.trans.winter[, , x] = elasticity.stochastic.rr.trans.winter[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.recruit.summer[, , x] = elasticity.mean.recruit.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.recruit.summer[, , x] = elasticity.sd.recruit.summer[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.recruit.summer[, , x] = elasticity.stochastic.recruit.summer[, , x]/(sim.n.years.asymptotic)
  
}




###########################################################################
#
# 4. Processing the results ----
#
###########################################################################

## 4.1. Check that the elasticities sum to 0 ---- 
# ------------------------------------------

elast.sums.vec = array(0, c(ncol(marmots.vr) - 1, n.sim))

for(x in 1:n.sim){
  
  elast.sums.vec[1, x] = sum(round(elasticity.mean.j.surv.winter[,, x] + elasticity.sd.j.surv.winter[,, x] - elasticity.stochastic.j.surv.winter[,, x], 10))
  
  elast.sums.vec[2, x] = sum(round(elasticity.mean.y.surv.winter[,, x] + elasticity.sd.y.surv.winter[,, x] - elasticity.stochastic.y.surv.winter[,, x], 10))
  elast.sums.vec[3, x] = sum(round(elasticity.mean.y.surv.summer[,, x] + elasticity.sd.y.surv.summer[,, x] - elasticity.stochastic.y.surv.summer[,, x], 10))
  
  elast.sums.vec[4, x] = sum(round(elasticity.mean.nr.surv.winter[,, x] + elasticity.sd.nr.surv.winter[,, x] - elasticity.stochastic.nr.surv.winter[,, x], 10))
  elast.sums.vec[5, x] = sum(round(elasticity.mean.nr.surv.summer[,, x] + elasticity.sd.nr.surv.summer[,, x] - elasticity.stochastic.nr.surv.summer[,, x], 10))
  
  elast.sums.vec[6, x] = sum(round(elasticity.mean.r.surv.winter[,, x] + elasticity.sd.r.surv.winter[,, x] - elasticity.stochastic.r.surv.winter[,, x], 10))
  elast.sums.vec[7, x] = sum(round(elasticity.mean.r.surv.summer[,, x] + elasticity.sd.r.surv.summer[,, x] - elasticity.stochastic.r.surv.summer[,, x], 10))
  
  elast.sums.vec[8, x] = sum(round(elasticity.mean.yr.trans.winter[,, x] + elasticity.sd.yr.trans.winter[,, x] - elasticity.stochastic.yr.trans.winter[,, x], 10))
  
  elast.sums.vec[9, x] = sum(round(elasticity.mean.nrr.trans.winter[,, x] + elasticity.sd.nrr.trans.winter[,, x] - elasticity.stochastic.nrr.trans.winter[,, x], 10))
  
  elast.sums.vec[10, x] = sum(round(elasticity.mean.rr.trans.winter[,, x] + elasticity.sd.rr.trans.winter[,, x] - elasticity.stochastic.rr.trans.winter[,, x], 10))
  
  elast.sums.vec[11, x] = sum(round(elasticity.mean.recruit.summer[,, x] + elasticity.sd.recruit.summer[,, x] - elasticity.stochastic.recruit.summer[,, x], 10))

}


## 4.2. Processing elasticities results ---- 
# -------------------------------------

## 4.2.1. Creating data frame with elasticities ----
# ---------------------------------------------

elast.df = expand.grid(vr = c("J survival", "Y survival (winter)", "Y survival (summer)", "NR survival (winter)", "NR survival (summer)", "R survival (winter)", "R survival (summer)", "Y-R transition", "NR-R transition", "R-R stasis", "Recruitment"), 
                       elast.type = c("Mean", "SD", "Stochastic"), 
                       elast.mean = NA, 
                       elast.low = NA, 
                       elast.up = NA)


## 4.2.2. Juvenile survival ----
# -------------------------

elast.df$elast.mean[which(elast.df$vr == "J survival")] = c(mean(abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum))),
                     mean(abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum))),
                     mean(abs(apply(elasticity.stochastic.j.surv.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "J survival")] = c(quantile(abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum)), probs = 0.025),
                     quantile(abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum)), probs = 0.025),
                     quantile(abs(apply(elasticity.stochastic.j.surv.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "J survival")] = c(quantile(abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                               quantile(abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                               quantile(abs(apply(elasticity.stochastic.j.surv.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.3. Yearling winter survival ----
# --------------------------------

elast.df$elast.mean[which(elast.df$vr == "Y survival (winter)")] = c(mean(abs(apply(elasticity.mean.y.surv.winter, 3, FUN = sum))),
                     mean(abs(apply(elasticity.sd.y.surv.winter, 3, FUN = sum))),
                     mean(abs(apply(elasticity.stochastic.y.surv.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "Y survival (winter)")] = c(quantile(abs(apply(elasticity.mean.y.surv.winter, 3, FUN = sum)), probs = 0.025),
                      quantile(abs(apply(elasticity.sd.y.surv.winter, 3, FUN = sum)), probs = 0.025),
                      quantile(abs(apply(elasticity.stochastic.y.surv.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "Y survival (winter)")] = c(quantile(abs(apply(elasticity.mean.y.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.sd.y.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.stochastic.y.surv.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.4. Yearling summer survival ----
# --------------------------------

elast.df$elast.mean[which(elast.df$vr == "Y survival (summer)")] = c(mean(abs(apply(elasticity.mean.y.surv.summer, 3, FUN = sum))),
                                                                                                         mean(abs(apply(elasticity.sd.y.surv.summer, 3, FUN = sum))),
                                                                                                         mean(abs(apply(elasticity.stochastic.y.surv.summer, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "Y survival (summer)")] = c(quantile(abs(apply(elasticity.mean.y.surv.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                        quantile(abs(apply(elasticity.sd.y.surv.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                        quantile(abs(apply(elasticity.stochastic.y.surv.summer, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "Y survival (summer)")] = c(quantile(abs(apply(elasticity.mean.y.surv.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.sd.y.surv.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.stochastic.y.surv.summer, 3, FUN = sum)), probs = 0.975))


## 4.2.5. Non reproductive winter survival ----
# ----------------------------------------

elast.df$elast.mean[which(elast.df$vr == "NR survival (winter)")] = c(mean(abs(apply(elasticity.mean.nr.surv.winter, 3, FUN = sum))),
                     mean(abs(apply(elasticity.sd.nr.surv.winter, 3, FUN = sum))),
                     mean(abs(apply(elasticity.stochastic.nr.surv.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "NR survival (winter)")] = c(quantile(abs(apply(elasticity.mean.nr.surv.winter, 3, FUN = sum)), probs = 0.025),
                      quantile(abs(apply(elasticity.sd.nr.surv.winter, 3, FUN = sum)), probs = 0.025),
                      quantile(abs(apply(elasticity.stochastic.nr.surv.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "NR survival (winter)")] = c(quantile(abs(apply(elasticity.mean.nr.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.sd.nr.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.stochastic.nr.surv.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.6. Non reproductive summer survival ----
# ----------------------------------------

elast.df$elast.mean[which(elast.df$vr == "NR survival (summer)")] = c(mean(abs(apply(elasticity.mean.nr.surv.summer, 3, FUN = sum))),
                                                                                                          mean(abs(apply(elasticity.sd.nr.surv.summer, 3, FUN = sum))),
                                                                                                          mean(abs(apply(elasticity.stochastic.nr.surv.summer, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "NR survival (summer)")] = c(quantile(abs(apply(elasticity.mean.nr.surv.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                         quantile(abs(apply(elasticity.sd.nr.surv.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                         quantile(abs(apply(elasticity.stochastic.nr.surv.summer, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "NR survival (summer)")] = c(quantile(abs(apply(elasticity.mean.nr.surv.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.sd.nr.surv.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.stochastic.nr.surv.summer, 3, FUN = sum)), probs = 0.975))


## 4.2.7. Reproductive winter survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "R survival (winter)")] = c(mean(abs(apply(elasticity.mean.r.surv.winter, 3, FUN = sum))),
                                                                                                          mean(abs(apply(elasticity.sd.r.surv.winter, 3, FUN = sum))),
                                                                                                          mean(abs(apply(elasticity.stochastic.r.surv.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "R survival (winter)")] = c(quantile(abs(apply(elasticity.mean.r.surv.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                         quantile(abs(apply(elasticity.sd.r.surv.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                         quantile(abs(apply(elasticity.stochastic.r.surv.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "R survival (winter)")] = c(quantile(abs(apply(elasticity.mean.r.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.sd.r.surv.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.stochastic.r.surv.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.8. Reproductive summer survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "R survival (summer)")] = c(mean(abs(apply(elasticity.mean.r.surv.summer, 3, FUN = sum))),
                                                                                                          mean(abs(apply(elasticity.sd.r.surv.summer, 3, FUN = sum))),
                                                                                                          mean(abs(apply(elasticity.stochastic.r.surv.summer, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "R survival (summer)")] = c(quantile(abs(apply(elasticity.mean.r.surv.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                         quantile(abs(apply(elasticity.sd.r.surv.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                         quantile(abs(apply(elasticity.stochastic.r.surv.summer, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "R survival (summer)")] = c(quantile(abs(apply(elasticity.mean.r.surv.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.sd.r.surv.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                        quantile(abs(apply(elasticity.stochastic.r.surv.summer, 3, FUN = sum)), probs = 0.975))


## 4.2.9. Y-R transition ----
# ----------------------

elast.df$elast.mean[which(elast.df$vr == "Y-R transition")] = c(mean(abs(apply(elasticity.mean.yr.trans.winter, 3, FUN = sum))),
                                                                                                         mean(abs(apply(elasticity.sd.yr.trans.winter, 3, FUN = sum))),
                                                                                                         mean(abs(apply(elasticity.stochastic.yr.trans.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "Y-R transition")] = c(quantile(abs(apply(elasticity.mean.yr.trans.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                        quantile(abs(apply(elasticity.sd.yr.trans.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                        quantile(abs(apply(elasticity.stochastic.yr.trans.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "Y-R transition")] = c(quantile(abs(apply(elasticity.mean.yr.trans.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                       quantile(abs(apply(elasticity.sd.yr.trans.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                       quantile(abs(apply(elasticity.stochastic.yr.trans.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.10. NR-R transition ----
# ------------------------

elast.df$elast.mean[which(elast.df$vr == "NR-R transition")] = c(mean(abs(apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum))),
                                                                                                    mean(abs(apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum))),
                                                                                                    mean(abs(apply(elasticity.stochastic.nrr.trans.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "NR-R transition")] = c(quantile(abs(apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                   quantile(abs(apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                   quantile(abs(apply(elasticity.stochastic.nrr.trans.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "NR-R transition")] = c(quantile(abs(apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                  quantile(abs(apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                  quantile(abs(apply(elasticity.stochastic.nrr.trans.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.11. R-R stasis ----
# -------------------

elast.df$elast.mean[which(elast.df$vr == "R-R stasis")] = c(mean(abs(apply(elasticity.mean.rr.trans.winter, 3, FUN = sum))),
                                                                                                    mean(abs(apply(elasticity.sd.rr.trans.winter, 3, FUN = sum))),
                                                                                                    mean(abs(apply(elasticity.stochastic.rr.trans.winter, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "R-R stasis")] = c(quantile(abs(apply(elasticity.mean.rr.trans.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                   quantile(abs(apply(elasticity.sd.rr.trans.winter, 3, FUN = sum)), probs = 0.025),
                                                                                                   quantile(abs(apply(elasticity.stochastic.rr.trans.winter, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "R-R stasis")] = c(quantile(abs(apply(elasticity.mean.rr.trans.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                  quantile(abs(apply(elasticity.sd.rr.trans.winter, 3, FUN = sum)), probs = 0.975),
                                                                                                  quantile(abs(apply(elasticity.stochastic.rr.trans.winter, 3, FUN = sum)), probs = 0.975))


## 4.2.12. Recruitment ----
# --------------------

elast.df$elast.mean[which(elast.df$vr == "Recruitment")] = c(mean(abs(apply(elasticity.mean.recruit.summer, 3, FUN = sum))),
                                                                                                    mean(abs(apply(elasticity.sd.recruit.summer, 3, FUN = sum))),
                                                                                                    mean(abs(apply(elasticity.stochastic.recruit.summer, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "Recruitment")] = c(quantile(abs(apply(elasticity.mean.recruit.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                   quantile(abs(apply(elasticity.sd.recruit.summer, 3, FUN = sum)), probs = 0.025),
                                                                                                   quantile(abs(apply(elasticity.stochastic.recruit.summer, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "Recruitment")] = c(quantile(abs(apply(elasticity.mean.recruit.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                  quantile(abs(apply(elasticity.sd.recruit.summer, 3, FUN = sum)), probs = 0.975),
                                                                                                  quantile(abs(apply(elasticity.stochastic.recruit.summer, 3, FUN = sum)), probs = 0.975))


## 4.3. Processing relative effect of variability results ---- 
# -------------------------------------------------------

# Following Morris et al. (2008), we compute the ratio: sum(elasticity to SD of all vital rates in one of the survival/transition/reproductive rates categories) / sum(elasticity to the mean + elasticity to SD of all vital rates)

## 4.3.1. Survival rates ----
# ----------------------

relative.elast.surv = abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum) +
                          apply(elasticity.sd.y.surv.summer, 3, FUN = sum) + 
                          apply(elasticity.sd.y.surv.winter, 3, FUN = sum) +
                          apply(elasticity.sd.nr.surv.summer, 3, FUN = sum) + 
                          apply(elasticity.sd.nr.surv.winter, 3, FUN = sum) +
                          apply(elasticity.sd.r.surv.summer, 3, FUN = sum) +
                          apply(elasticity.sd.r.surv.winter, 3, FUN = sum)) / (abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.y.surv.summer, 3, FUN = sum) + 
                                                                                   apply(elasticity.sd.y.surv.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.nr.surv.summer, 3, FUN = sum) + 
                                                                                   apply(elasticity.sd.nr.surv.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.r.surv.summer, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.r.surv.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.yr.trans.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.rr.trans.winter, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.recruit.summer, 3, FUN = sum)) + abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.y.surv.summer, 3, FUN = sum) + 
                                                                                                                                            apply(elasticity.mean.y.surv.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.nr.surv.summer, 3, FUN = sum) + 
                                                                                                                                            apply(elasticity.mean.nr.surv.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.r.surv.summer, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.r.surv.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.yr.trans.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.rr.trans.winter, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.recruit.summer, 3, FUN = sum)))


## 4.3.2. Reproductive rates ----
# --------------------------

relative.elast.repro = abs(apply(elasticity.sd.recruit.summer, 3, FUN = sum)) / (abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.y.surv.summer, 3, FUN = sum) + 
                                                                                     apply(elasticity.sd.y.surv.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.nr.surv.summer, 3, FUN = sum) + 
                                                                                     apply(elasticity.sd.nr.surv.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.r.surv.summer, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.r.surv.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.yr.trans.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.rr.trans.winter, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.recruit.summer, 3, FUN = sum)) + abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.y.surv.summer, 3, FUN = sum) + 
                                                                                                                                              apply(elasticity.mean.y.surv.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.nr.surv.summer, 3, FUN = sum) + 
                                                                                                                                              apply(elasticity.mean.nr.surv.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.r.surv.summer, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.r.surv.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.yr.trans.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.rr.trans.winter, 3, FUN = sum) +
                                                                                                                                              apply(elasticity.mean.recruit.summer, 3, FUN = sum)))


## 4.3.3. Transition rates ----
# ------------------------

relative.elast.trans = abs(apply(elasticity.sd.yr.trans.winter, 3, FUN = sum) +
                           apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum) +
                           apply(elasticity.sd.rr.trans.winter, 3, FUN = sum)) / (abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.y.surv.summer, 3, FUN = sum) + 
                                                                                      apply(elasticity.sd.y.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.nr.surv.summer, 3, FUN = sum) + 
                                                                                      apply(elasticity.sd.nr.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.r.surv.summer, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.r.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.yr.trans.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.rr.trans.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.recruit.summer, 3, FUN = sum)) + abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.y.surv.summer, 3, FUN = sum) + 
                                                                                                                                               apply(elasticity.mean.y.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.nr.surv.summer, 3, FUN = sum) + 
                                                                                                                                               apply(elasticity.mean.nr.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.r.surv.summer, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.r.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.yr.trans.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.rr.trans.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.recruit.summer, 3, FUN = sum)))


## 4.3.4. All vital rates ----
# -----------------------

relative.elast.all.vr = abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum) +
                            apply(elasticity.sd.y.surv.summer, 3, FUN = sum) + 
                            apply(elasticity.sd.y.surv.winter, 3, FUN = sum) +
                            apply(elasticity.sd.nr.surv.summer, 3, FUN = sum) + 
                            apply(elasticity.sd.nr.surv.winter, 3, FUN = sum) +
                            apply(elasticity.sd.r.surv.summer, 3, FUN = sum) +
                            apply(elasticity.sd.r.surv.winter, 3, FUN = sum) +
                            apply(elasticity.sd.yr.trans.winter, 3, FUN = sum) +
                            apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum) +
                            apply(elasticity.sd.rr.trans.winter, 3, FUN = sum) +
                            apply(elasticity.sd.recruit.summer, 3, FUN = sum)) / (abs(apply(elasticity.sd.j.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.y.surv.summer, 3, FUN = sum) + 
                                                                                      apply(elasticity.sd.y.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.nr.surv.summer, 3, FUN = sum) + 
                                                                                      apply(elasticity.sd.nr.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.r.surv.summer, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.r.surv.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.yr.trans.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.nrr.trans.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.rr.trans.winter, 3, FUN = sum) +
                                                                                      apply(elasticity.sd.recruit.summer, 3, FUN = sum)) + abs(apply(elasticity.mean.j.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.y.surv.summer, 3, FUN = sum) + 
                                                                                                                                               apply(elasticity.mean.y.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.nr.surv.summer, 3, FUN = sum) + 
                                                                                                                                               apply(elasticity.mean.nr.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.r.surv.summer, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.r.surv.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.yr.trans.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.nrr.trans.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.rr.trans.winter, 3, FUN = sum) +
                                                                                                                                               apply(elasticity.mean.recruit.summer, 3, FUN = sum)))

relative.elast.df = data.frame(vr = c("Survival", "Reproduction", "Transition", "All vital rates"),
                               relative.elast.mean = c(mean(relative.elast.surv), mean(relative.elast.repro), mean(relative.elast.trans), mean(relative.elast.all.vr)), 
                               relative.elast.low = c(quantile(relative.elast.surv, probs = 0.025), quantile(relative.elast.repro, probs = 0.025), quantile(relative.elast.trans, probs = 0.025), quantile(relative.elast.all.vr, probs = 0.025)), 
                               relative.elast.up = c(quantile(relative.elast.surv, probs = 0.975), quantile(relative.elast.repro, probs = 0.975), quantile(relative.elast.trans, probs = 0.975), quantile(relative.elast.all.vr, probs = 0.975)))




###########################################################################
#
# 5. Plotting the results -----
#
###########################################################################

## 5.1. Preparing plotting ----
# ------------------------

elast.df$elast.type = factor(elast.df$elast.type, levels = c("Stochastic", "Mean", "SD"))
lbs = setNames(c(expression(bold("Stochastic"*"("*E^S*")")), expression(bold("Stochastic"*"("*E^S*")")), expression(bold("Stochastic"*"("*E^S*")"))), c("Stochastic", "Mean", "SD"))[levels(elast.df$elast.type)]
labels.elast.type = c(expression(bold("Stochastic ("*E^S*")")), expression(bold("Stochastic ("*E^S*")")), expression(bold("Stochastic ("*E^S*")")))

elast.df$elast.type.parsed = c(rep("S\U003BC", 11), rep("S\u03C3", 11), rep("S", 11))


## 5.2. Elasticities ----
# ------------------

png(filename = "Marmots_ElastAll.png",
    width = 3000,
    height = 2500,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

ggplot(elast.df, aes(x = vr, y = elast.mean, shape = elast.type.parsed)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = elast.low, ymax = elast.up, width = 0.25)) +
  labs(x = "Vital rate", y="Absolute elasticity") +
  scale_shape_manual(values = c(16, 17, 15),
                     labels = c(expression(paste(E^"S\U003BC")), expression(paste(E^"S\u03C3")), expression(paste(E^"S  "))), 
                     name = "Elasticity type") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "right", legend.key.size = unit(2, "lines"), 
        strip.text.x = element_text(size = 20, face = "bold"))

dev.off()


## 5.3. Relative effect of variability ----
# ------------------------------------

png(filename = "Marmots_RelativeEffectVariability.png",
    width = 4000,
    height = 3000,
    units = "px",
    bg = "white",
    res = 300,
    type = "cairo")

ggplot(relative.elast.df, aes(x = vr, y = relative.elast.mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = relative.elast.low, ymax = relative.elast.up, width = 0.25)) +
  labs(x = "Vital rate", y="Relative effect of variability") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 25, colour = "black", margin = margin(t = 18, r = 0, b = 0, l = 0), face = "bold"), 
        axis.title.y = element_text(size = 25, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"), 
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), hjust = 0.5), 
        axis.text.y = element_text(size = 22, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)), legend.text = element_text(size = 22), 
        legend.title = element_text(face = "bold", size = 25), 
        legend.position = "right", 
        legend.key.size = unit(6, "lines"),
        strip.text.x = element_text(size = 20, face = "bold"))

dev.off()




###########################################################################
#
# 6. Saving the results -----
#
###########################################################################

write.csv(elast.df, "Marmots_Elasticities.csv", row.names = F)
write.csv(relative.elast.df, "Marmots_RelativeElasticities.csv", row.names = F)
