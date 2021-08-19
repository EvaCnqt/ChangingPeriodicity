############################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script uses the seasonal vital rate-estimates for the meerkat population. 
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
  library(viridis)
  library(bbmle)
  library(lubridate)
  library(lme4)
  library(ggplot2)
}

load.librairies()


## 1.3. Loading and preparing data ----
# --------------------------------  
data.meerkats = read.csv("RScripts_ForSubmission/MeerkatsData.csv", stringsAsFactors = F, na.strings = c("", "NA"))
head(data.meerkats)

load("GLMM_survJ.RData")
load("GLMM_survS.RData")
load("GLMM_survH.RData")
load("GLMM_survD.RData")

load("GLMM_transHD.RData")
load("GLMM_emig.RData")

load("GLMM_recruitH.RData")
load("GLMM_recruitD.RData")

load("GLMM_poprange.RData")




###########################################################################
#
# 2. Preparing the simulations ----
#
###########################################################################

## 2.1. Creating functions to create the MPMs ----
# -------------------------------------------

# Wet season

buildmat <- function(survJ, survS, survH, survD, transHD, emigH, recruitH, recruitD){
  
  # Creating the matrix
  mat = matrix(0, nrow = 4, ncol = 4)
  colnames(mat) = c("J", "S", "H", "D")
  rownames(mat) = colnames(mat)
  
  # Juveniles
  mat["J", "H"] = recruitH
  mat["J", "D"] = recruitD
  
  # Subadults
  mat["S", "J"] = survJ
  
  # Helpers
  mat["H", "S"] = survS
  mat["H", "H"] = survH * (1 - emigH) * (1 - transHD)
  
  # Dominants 
  mat["D", "H"] = survH * (1 - emigH) * transHD
  mat["D", "D"] = survD
  
  return(mat)
}


## 2.2. Getting mean vital rates ----
# ------------------------------

data.meerkats$year = as.numeric(as.character(data.meerkats$year))
data.meerkats$season = as.character(data.meerkats$season)
year.sim = as.numeric(as.character(unique(data.meerkats$year)))
n.years = length(year.sim)

# Annual vital rates

annual.surv.j.dry    = rep(0, n.years)
annual.surv.j.wet    = rep(0, n.years)
annual.surv.sa.dry   = rep(0, n.years)
annual.surv.sa.wet   = rep(0, n.years)
annual.surv.h.dry    = rep(0, n.years)
annual.surv.h.wet    = rep(0, n.years)
annual.surv.d.dry    = rep(0, n.years)
annual.surv.d.wet    = rep(0, n.years)
annual.trans.hd.dry  = rep(0, n.years)
annual.trans.hd.wet  = rep(0, n.years)
annual.emig.h.dry    = rep(0, n.years)
annual.emig.h.wet    = rep(0, n.years)
annual.recruit.h.dry = rep(0, n.years)
annual.recruit.h.wet = rep(0, n.years)
annual.recruit.d.dry = rep(0, n.years)
annual.recruit.d.wet = rep(0, n.years)

for(y in 1:n.years){
  
  year = year.sim[y]
  
  pop.dens.dry = unique(data.meerkats$density[which(data.meerkats$season == "dry" & data.meerkats$year == year)])
  pop.dens.rain = unique(data.meerkats$density[which(data.meerkats$season == "rain" & data.meerkats$year == year)])
  
  annual.surv.j.dry[y] = predict(survJ, newdata = data.frame(season = "dry", density = pop.dens.dry, density2 = pop.dens.dry^2, year = year), type = "response")
  annual.surv.j.wet[y] = predict(survJ, newdata = data.frame(season = "rain", density = pop.dens.rain, density2 = pop.dens.rain^2, year = year), type = "response")
  
  annual.surv.sa.dry[y] = predict(survS, newdata = data.frame(season = "dry", density = pop.dens.dry, year = year), type = "response")
  annual.surv.sa.wet[y] = predict(survS, newdata = data.frame(season = "rain", density = pop.dens.rain, year = year), type = "response")
  
  annual.surv.h.dry[y] = predict(survH, newdata = data.frame(season = "dry", year = year), type = "response")
  annual.surv.h.wet[y] = predict(survH, newdata = data.frame(season = "rain", year = year), type = "response")
  
  annual.surv.d.dry[y] = predict(survD, newdata = data.frame(season = "dry", year = year), type = "response")
  annual.surv.d.wet[y] = predict(survD, newdata = data.frame(season = "rain", year = year), type = "response")
  
  annual.trans.hd.dry[y] = predict(transition, newdata = data.frame(season = "dry", density = pop.dens.dry, year = year), type = "response")
  annual.trans.hd.wet[y] = predict(transition, newdata = data.frame(season = "rain", density = pop.dens.rain, year = year), type = "response")
  
  annual.emig.h.dry[y] = predict(emig, newdata = data.frame(season = "dry", density = pop.dens.dry, year = year), type = "response")
  annual.emig.h.wet[y] = predict(emig, newdata = data.frame(season = "rain", density = pop.dens.rain, year = year), type = "response")
  
  annual.recruit.h.dry[y] = predict(recruitment.H, newdata = data.frame(season = c("dry", "rain"), density = pop.dens.dry, year = year), type = "response")[1]
  annual.recruit.h.wet[y] = predict(recruitment.H, newdata = data.frame(season = c("dry", "rain"), density = pop.dens.rain, year = year), type = "response")[2]
  
  annual.recruit.d.dry[y] = predict(recruitment.D, newdata = data.frame(season = c("dry", "rain"), density = pop.dens.dry, density2 = pop.dens.dry^2, year = year), type = "response")[1]
  annual.recruit.d.wet[y] = predict(recruitment.D, newdata = data.frame(season = c("dry", "rain"), density = pop.dens.rain, density2 = pop.dens.rain^2, year = year), type = "response")[2]
  
}

mean.surv.j.dry    = mean(annual.surv.j.dry)
mean.surv.j.wet    = mean(annual.surv.j.wet)
mean.surv.sa.dry   = mean(annual.surv.sa.dry)
mean.surv.sa.wet   = mean(annual.surv.sa.wet)
mean.surv.h.dry    = mean(annual.surv.h.dry)
mean.surv.h.wet    = mean(annual.surv.h.wet)
mean.surv.d.dry    = mean(annual.surv.d.dry)
mean.surv.d.wet    = mean(annual.surv.d.wet)
mean.trans.hd.dry  = mean(annual.trans.hd.dry)
mean.trans.hd.wet  = mean(annual.trans.hd.wet)
mean.emig.h.dry    = mean(annual.emig.h.dry)
mean.emig.h.wet    = mean(annual.emig.h.wet)
mean.recruit.h.dry = mean(annual.recruit.h.dry)
mean.recruit.h.wet = mean(annual.recruit.h.wet)
mean.recruit.d.dry = mean(annual.recruit.d.dry)
mean.recruit.d.wet = mean(annual.recruit.d.wet)


## 2.3. Creating yearly seasonal MPMs ----
# -----------------------------------

dry.matrices = array(0, c(length(unique(data.meerkats$stage)), length(unique(data.meerkats$stage)), n.years))
wet.matrices = array(0, c(length(unique(data.meerkats$stage)), length(unique(data.meerkats$stage)), n.years))

for(y in 1:n.years){
  
  dry.matrices[, , y] = buildmat(survJ = annual.surv.j.dry[y], 
                                 survS = annual.surv.sa.dry[y], 
                                 survH = annual.surv.h.dry[y], 
                                 survD = annual.surv.d.dry[y],
                                 transHD = annual.trans.hd.dry[y],
                                 emigH = annual.emig.h.dry[y],
                                 recruitH = annual.recruit.h.dry[y],
                                 recruitD = annual.recruit.d.dry[y])
  
  wet.matrices[, , y] = buildmat(survJ = annual.surv.j.wet[y], 
                                 survS = annual.surv.sa.wet[y], 
                                 survH = annual.surv.h.wet[y], 
                                 survD = annual.surv.d.wet[y],
                                 transHD = annual.trans.hd.wet[y],
                                 emigH = annual.emig.h.wet[y],
                                 recruitH = annual.recruit.h.wet[y],
                                 recruitD = annual.recruit.d.wet[y])
  
}




###########################################################################
#
# 3. Simulating dynamics to calculate elasticities ----
#
###########################################################################

## 3.1. Prepare simulations ----
# -------------------------

# Define population vector. We start with intial population vector and density picked from 2016 dry season
nbJ.init = length(data.meerkats$stage[data.meerkats$stage == "J" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])
nbS.init = length(data.meerkats$stage[data.meerkats$stage == "S" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])
nbH.init = length(data.meerkats$stage[data.meerkats$stage == "H" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])
nbD.init = length(data.meerkats$stage[data.meerkats$stage == "D" & data.meerkats$year == "2016" & data.meerkats$season == "dry"])

n0 = c(nbJ.init, nbS.init, nbH.init, nbD.init)


# Simulations parameters

n.sim = 100
sim.n.years.total = 1000
sim.n.years.discard = 500
sim.n.years.asymptotic = sim.n.years.total - sim.n.years.discard
t2 <- sim.n.years.asymptotic + 2 * sim.n.years.discard # Time vector for the backward iteration to calculate the left eigenvectors

# Elasticities storage

elasticity.mean.surv.j.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.j.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.sa.dry   = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.sa.wet   = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.h.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.h.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.d.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.surv.d.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.trans.hd.dry  = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.trans.hd.wet  = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.emig.h.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.emig.h.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.recruit.h.dry = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.recruit.h.wet = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.recruit.d.dry = array(0, c(length(n0), length(n0), n.sim))
elasticity.mean.recruit.d.wet = array(0, c(length(n0), length(n0), n.sim))

elasticity.sd.surv.j.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.j.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.sa.dry   = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.sa.wet   = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.h.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.h.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.d.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.surv.d.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.trans.hd.dry  = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.trans.hd.wet  = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.emig.h.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.emig.h.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.recruit.h.dry = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.recruit.h.wet = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.recruit.d.dry = array(0, c(length(n0), length(n0), n.sim))
elasticity.sd.recruit.d.wet = array(0, c(length(n0), length(n0), n.sim))

elasticity.stochastic.surv.j.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.j.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.sa.dry   = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.sa.wet   = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.h.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.h.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.d.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.surv.d.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.trans.hd.dry  = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.trans.hd.wet  = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.emig.h.dry    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.emig.h.wet    = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.recruit.h.dry = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.recruit.h.wet = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.recruit.d.dry = array(0, c(length(n0), length(n0), n.sim))
elasticity.stochastic.recruit.d.wet = array(0, c(length(n0), length(n0), n.sim))


## 3.2. Simulate dynamics ----
# -----------------------

for(x in 1:n.sim){ 
  
  # Initialize list of C matrices (with difference between perturbed and unperturbed matrices)
  
  mean.C.mat.surv.j.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.j.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.sa.dry   = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.sa.wet   = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.h.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.h.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.d.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.surv.d.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.trans.hd.dry  = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.trans.hd.wet  = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.emig.h.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.emig.h.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.recruit.h.dry = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.recruit.h.wet = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.recruit.d.dry = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  mean.C.mat.recruit.d.wet = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  
  sd.C.mat.surv.j.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.j.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.sa.dry   = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.sa.wet   = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.h.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.h.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.d.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.surv.d.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.trans.hd.dry  = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.trans.hd.wet  = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.emig.h.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.emig.h.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.recruit.h.dry = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.recruit.h.wet = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.recruit.d.dry = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  sd.C.mat.recruit.d.wet = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  
  vr.value.C.mat.surv.j.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.j.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.sa.dry   = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.sa.wet   = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.h.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.h.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.d.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.surv.d.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.trans.hd.dry  = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.trans.hd.wet  = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.emig.h.dry    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.emig.h.wet    = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.recruit.h.dry = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.recruit.h.wet = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.recruit.d.dry = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  vr.value.C.mat.recruit.d.wet = array(0, c(length(n0), c(length(n0)), sim.n.years.total))
  
  # Initialize population vectors (n0)
  
  vec1 = n0
  vec1 <- vec1 / sum(vec1)
  vec1 <- t(vec1) 
  vec2 <- vec1
  
  uvecs <- array(0, c(length(n0), sim.n.years.asymptotic)) # Right eigenvectors
  vvecs <- array(0, c(length(n0), sim.n.years.asymptotic)) # Left eigenvectors
  growth <- array(0, sim.n.years.asymptotic)
  
  # Simulate environmental states, i.e., years over trun simulations
  
  year.sim = sample(unique(data.meerkats$year), t2 + 1, replace = T)
  
  for(i in 1:sim.n.years.total){
    
    year = year.sim[i]
    year.index = which(unique(data.meerkats$year) == year)
    
    # Create annual matrix
    annual.mat = dry.matrices[, , year.index] %*% wet.matrices[, , year.index]
    
    # Perturb each vital rate with the mean vital rate, with the standard deviation, and the value of the vital rate.

    ################################

    # Juvenile dry-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.j.dry = annual.surv.j.dry[year.index] + mean.surv.j.dry
    
    mean.perturbed.annual.mat.surv.j.dry = buildmat(survJ = mean.perturbed.surv.j.dry, 
                                                    survS = annual.surv.sa.dry[year.index],
                                                    survH = annual.surv.h.dry[year.index],
                                                    survD = annual.surv.d.dry[year.index],
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = annual.emig.h.dry[year.index],
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = annual.surv.d.wet[year.index],
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index],
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.j.dry[, , i] = mean.perturbed.annual.mat.surv.j.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.j.dry = annual.surv.j.dry[year.index] + (annual.surv.j.dry[year.index] - mean.surv.j.dry)
    
    sd.perturbed.annual.mat.surv.j.dry = buildmat(survJ = sd.perturbed.surv.j.dry,
                                                  survS = annual.surv.sa.dry[year.index],
                                                  survH = annual.surv.h.dry[year.index],
                                                  survD = annual.surv.d.dry[year.index],
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = annual.emig.h.dry[year.index],
                                                  recruitH = annual.recruit.h.dry[year.index],
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = annual.surv.h.wet[year.index],
                                                                                                            survD = annual.surv.d.wet[year.index],
                                                                                                            transHD = annual.trans.hd.wet[year.index],
                                                                                                            emigH = annual.emig.h.wet[year.index],
                                                                                                            recruitH = annual.recruit.h.wet[year.index],
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.j.dry[, , i] = sd.perturbed.annual.mat.surv.j.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.j.dry = annual.surv.j.dry[year.index] + annual.surv.j.dry[year.index]
    
    vr.value.perturbed.annual.mat.surv.j.dry = buildmat(survJ = vr.value.perturbed.surv.j.dry,
                                                        survS = annual.surv.sa.dry[year.index],
                                                        survH = annual.surv.h.dry[year.index],
                                                        survD = annual.surv.d.dry[year.index],
                                                        transHD = annual.trans.hd.dry[year.index],
                                                        emigH = annual.emig.h.dry[year.index],
                                                        recruitH = annual.recruit.h.dry[year.index],
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                  survS = annual.surv.sa.wet[year.index],
                                                                                                                  survH = annual.surv.h.wet[year.index],
                                                                                                                  survD = annual.surv.d.wet[year.index],
                                                                                                                  transHD = annual.trans.hd.wet[year.index],
                                                                                                                  emigH = annual.emig.h.wet[year.index],
                                                                                                                  recruitH = annual.recruit.h.wet[year.index],
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.j.dry[, , i] = vr.value.perturbed.annual.mat.surv.j.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Juvenile wet-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.j.wet = annual.surv.j.wet[year.index] + mean.surv.j.wet
    
    mean.perturbed.annual.mat.surv.j.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                    survS = annual.surv.sa.dry[year.index],
                                                    survH = annual.surv.h.dry[year.index],
                                                    survD = annual.surv.d.dry[year.index],
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = annual.emig.h.dry[year.index],
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = mean.perturbed.surv.j.wet,
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = annual.surv.d.wet[year.index],
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index],
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.j.wet[, , i] = mean.perturbed.annual.mat.surv.j.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.j.wet = annual.surv.j.wet[year.index] + (annual.surv.j.wet[year.index] - mean.surv.j.wet)
    
    sd.perturbed.annual.mat.surv.j.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                  survS = annual.surv.sa.dry[year.index],
                                                  survH = annual.surv.h.dry[year.index],
                                                  survD = annual.surv.d.dry[year.index],
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = annual.emig.h.dry[year.index],
                                                  recruitH = annual.recruit.h.dry[year.index],
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = sd.perturbed.surv.j.wet,
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = annual.surv.h.wet[year.index],
                                                                                                            survD = annual.surv.d.wet[year.index],
                                                                                                            transHD = annual.trans.hd.wet[year.index],
                                                                                                            emigH = annual.emig.h.wet[year.index],
                                                                                                            recruitH = annual.recruit.h.wet[year.index],
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.j.wet[, , i] = sd.perturbed.annual.mat.surv.j.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.j.wet = annual.surv.j.wet[year.index] + annual.surv.j.wet[year.index]
    
    vr.value.perturbed.annual.mat.surv.j.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                        survS = annual.surv.sa.dry[year.index],
                                                        survH = annual.surv.h.dry[year.index],
                                                        survD = annual.surv.d.dry[year.index],
                                                        transHD = annual.trans.hd.dry[year.index],
                                                        emigH = annual.emig.h.dry[year.index],
                                                        recruitH = annual.recruit.h.dry[year.index],
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = vr.value.perturbed.surv.j.wet,
                                                                                                                  survS = annual.surv.sa.wet[year.index],
                                                                                                                  survH = annual.surv.h.wet[year.index],
                                                                                                                  survD = annual.surv.d.wet[year.index],
                                                                                                                  transHD = annual.trans.hd.wet[year.index],
                                                                                                                  emigH = annual.emig.h.wet[year.index],
                                                                                                                  recruitH = annual.recruit.h.wet[year.index],
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.j.wet[, , i] = vr.value.perturbed.annual.mat.surv.j.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Subadult dry-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.sa.dry = annual.surv.sa.dry[year.index] + mean.surv.sa.dry
    
    mean.perturbed.annual.mat.surv.sa.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                     survS = mean.perturbed.surv.sa.dry,
                                                     survH = annual.surv.h.dry[year.index],
                                                     survD = annual.surv.d.dry[year.index],
                                                     transHD = annual.trans.hd.dry[year.index],
                                                     emigH = annual.emig.h.dry[year.index],
                                                     recruitH = annual.recruit.h.dry[year.index],
                                                     recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                               survS = annual.surv.sa.wet[year.index],
                                                                                                               survH = annual.surv.h.wet[year.index],
                                                                                                               survD = annual.surv.d.wet[year.index],
                                                                                                               transHD = annual.trans.hd.wet[year.index],
                                                                                                               emigH = annual.emig.h.wet[year.index],
                                                                                                               recruitH = annual.recruit.h.wet[year.index],
                                                                                                               recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.sa.dry[, , i] = mean.perturbed.annual.mat.surv.sa.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.sa.dry = annual.surv.sa.dry[year.index] + (annual.surv.sa.dry[year.index] - mean.surv.sa.dry)
    
    sd.perturbed.annual.mat.surv.sa.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                   survS = sd.perturbed.surv.sa.dry,
                                                   survH = annual.surv.h.dry[year.index],
                                                   survD = annual.surv.d.dry[year.index],
                                                   transHD = annual.trans.hd.dry[year.index],
                                                   emigH = annual.emig.h.dry[year.index],
                                                   recruitH = annual.recruit.h.dry[year.index],
                                                   recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                             survS = annual.surv.sa.wet[year.index],
                                                                                                             survH = annual.surv.h.wet[year.index],
                                                                                                             survD = annual.surv.d.wet[year.index],
                                                                                                             transHD = annual.trans.hd.wet[year.index],
                                                                                                             emigH = annual.emig.h.wet[year.index],
                                                                                                             recruitH = annual.recruit.h.wet[year.index],
                                                                                                             recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.sa.dry[, , i] = sd.perturbed.annual.mat.surv.sa.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.sa.dry = annual.surv.sa.dry[year.index] + annual.surv.sa.dry[year.index]
    
    vr.value.perturbed.annual.mat.surv.sa.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                         survS = vr.value.perturbed.surv.sa.dry,
                                                         survH = annual.surv.h.dry[year.index],
                                                         survD = annual.surv.d.dry[year.index],
                                                         transHD = annual.trans.hd.dry[year.index],
                                                         emigH = annual.emig.h.dry[year.index],
                                                         recruitH = annual.recruit.h.dry[year.index],
                                                         recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                   survS = annual.surv.sa.wet[year.index],
                                                                                                                   survH = annual.surv.h.wet[year.index],
                                                                                                                   survD = annual.surv.d.wet[year.index],
                                                                                                                   transHD = annual.trans.hd.wet[year.index],
                                                                                                                   emigH = annual.emig.h.wet[year.index],
                                                                                                                   recruitH = annual.recruit.h.wet[year.index],
                                                                                                                   recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.sa.dry[, , i] = vr.value.perturbed.annual.mat.surv.sa.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Subadult wet-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.sa.wet = annual.surv.sa.wet[year.index] + mean.surv.sa.wet
    
    mean.perturbed.annual.mat.surv.sa.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                     survS = annual.surv.sa.dry[year.index],
                                                     survH = annual.surv.h.dry[year.index],
                                                     survD = annual.surv.d.dry[year.index],
                                                     transHD = annual.trans.hd.dry[year.index],
                                                     emigH = annual.emig.h.dry[year.index],
                                                     recruitH = annual.recruit.h.dry[year.index],
                                                     recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                               survS = mean.perturbed.surv.sa.wet,
                                                                                                               survH = annual.surv.h.wet[year.index],
                                                                                                               survD = annual.surv.d.wet[year.index],
                                                                                                               transHD = annual.trans.hd.wet[year.index],
                                                                                                               emigH = annual.emig.h.wet[year.index],
                                                                                                               recruitH = annual.recruit.h.wet[year.index],
                                                                                                               recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.sa.wet[, , i] = mean.perturbed.annual.mat.surv.sa.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.sa.wet = annual.surv.sa.wet[year.index] + (annual.surv.sa.wet[year.index] - mean.surv.sa.wet)
    
    sd.perturbed.annual.mat.surv.sa.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                   survS = annual.surv.sa.dry[year.index],
                                                   survH = annual.surv.h.dry[year.index],
                                                   survD = annual.surv.d.dry[year.index],
                                                   transHD = annual.trans.hd.dry[year.index],
                                                   emigH = annual.emig.h.dry[year.index],
                                                   recruitH = annual.recruit.h.dry[year.index],
                                                   recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                             survS = sd.perturbed.surv.sa.wet,
                                                                                                             survH = annual.surv.h.wet[year.index],
                                                                                                             survD = annual.surv.d.wet[year.index],
                                                                                                             transHD = annual.trans.hd.wet[year.index],
                                                                                                             emigH = annual.emig.h.wet[year.index],
                                                                                                             recruitH = annual.recruit.h.wet[year.index],
                                                                                                             recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.sa.wet[, , i] = sd.perturbed.annual.mat.surv.sa.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.sa.wet = annual.surv.sa.wet[year.index] + annual.surv.sa.wet[year.index]
    
    vr.value.perturbed.annual.mat.surv.sa.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                         survS = annual.surv.sa.dry[year.index],
                                                         survH = annual.surv.h.dry[year.index],
                                                         survD = annual.surv.d.dry[year.index],
                                                         transHD = annual.trans.hd.dry[year.index],
                                                         emigH = annual.emig.h.dry[year.index],
                                                         recruitH = annual.recruit.h.dry[year.index],
                                                         recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                   survS = vr.value.perturbed.surv.sa.wet, 
                                                                                                                   survH = annual.surv.h.wet[year.index],
                                                                                                                   survD = annual.surv.d.wet[year.index],
                                                                                                                   transHD = annual.trans.hd.wet[year.index],
                                                                                                                   emigH = annual.emig.h.wet[year.index],
                                                                                                                   recruitH = annual.recruit.h.wet[year.index],
                                                                                                                   recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.sa.wet[, , i] = vr.value.perturbed.annual.mat.surv.sa.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper dry-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.h.dry = annual.surv.h.dry[year.index] + mean.surv.h.dry
    
    mean.perturbed.annual.mat.surv.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                    survS = annual.surv.sa.dry[year.index],
                                                    survH = mean.perturbed.surv.h.dry,
                                                    survD = annual.surv.d.dry[year.index],
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = annual.emig.h.dry[year.index],
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = annual.surv.d.wet[year.index],
                                                                                                              transHD = annual.trans.hd.wet[year.index], 
                                                                                                              emigH = annual.emig.h.wet[year.index],
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.h.dry[, , i] = mean.perturbed.annual.mat.surv.h.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.h.dry = annual.surv.h.dry[year.index] + (annual.surv.h.dry[year.index] - mean.surv.h.dry)
    
    sd.perturbed.annual.mat.surv.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                  survS = annual.surv.sa.dry[year.index],
                                                  survH = sd.perturbed.surv.h.dry, 
                                                  survD = annual.surv.d.dry[year.index],
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = annual.emig.h.dry[year.index], 
                                                  recruitH = annual.recruit.h.dry[year.index], 
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = annual.surv.h.wet[year.index], 
                                                                                                            survD = annual.surv.d.wet[year.index],
                                                                                                            transHD = annual.trans.hd.wet[year.index], 
                                                                                                            emigH = annual.emig.h.wet[year.index],
                                                                                                            recruitH = annual.recruit.h.wet[year.index],
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.h.dry[, , i] = sd.perturbed.annual.mat.surv.h.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.h.dry = annual.surv.h.dry[year.index] + annual.surv.h.dry[year.index]
    
    vr.value.perturbed.annual.mat.surv.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                        survS = annual.surv.sa.dry[year.index],
                                                        survH = vr.value.perturbed.surv.h.dry,
                                                        survD = annual.surv.d.dry[year.index],
                                                        transHD = annual.trans.hd.dry[year.index],
                                                        emigH = annual.emig.h.dry[year.index],
                                                        recruitH = annual.recruit.h.dry[year.index],
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                  survS = annual.surv.sa.wet[year.index], 
                                                                                                                  survH = annual.surv.h.wet[year.index], 
                                                                                                                  survD = annual.surv.d.wet[year.index],
                                                                                                                  transHD = annual.trans.hd.wet[year.index],
                                                                                                                  emigH = annual.emig.h.wet[year.index], 
                                                                                                                  recruitH = annual.recruit.h.wet[year.index],
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.h.dry[, , i] = vr.value.perturbed.annual.mat.surv.h.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper wet-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.h.wet = annual.surv.h.wet[year.index] + mean.surv.h.wet
    
    mean.perturbed.annual.mat.surv.h.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                    survS = annual.surv.sa.dry[year.index],
                                                    survH = annual.surv.h.dry[year.index],
                                                    survD = annual.surv.d.dry[year.index],
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = annual.emig.h.dry[year.index], 
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = mean.perturbed.surv.h.wet, 
                                                                                                              survD = annual.surv.d.wet[year.index], 
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index], 
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.h.wet[, , i] = mean.perturbed.annual.mat.surv.h.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.h.wet = annual.surv.h.wet[year.index] + (annual.surv.h.wet[year.index] - mean.surv.h.wet)
    
    sd.perturbed.annual.mat.surv.h.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                  survS = annual.surv.sa.dry[year.index],
                                                  survH = annual.surv.h.dry[year.index],
                                                  survD = annual.surv.d.dry[year.index],
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = annual.emig.h.dry[year.index],
                                                  recruitH = annual.recruit.h.dry[year.index],
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = sd.perturbed.surv.h.wet, 
                                                                                                            survD = annual.surv.d.wet[year.index],
                                                                                                            transHD = annual.trans.hd.wet[year.index], 
                                                                                                            emigH = annual.emig.h.wet[year.index], 
                                                                                                            recruitH = annual.recruit.h.wet[year.index],
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.h.wet[, , i] = sd.perturbed.annual.mat.surv.h.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.h.wet = annual.surv.h.wet[year.index] + annual.surv.h.wet[year.index]
    
    vr.value.perturbed.annual.mat.surv.h.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                        survS = annual.surv.sa.dry[year.index],
                                                        survH = annual.surv.h.dry[year.index],
                                                        survD = annual.surv.d.dry[year.index],
                                                        transHD = annual.trans.hd.dry[year.index], 
                                                        emigH = annual.emig.h.dry[year.index], 
                                                        recruitH = annual.recruit.h.dry[year.index], 
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                  survS = annual.surv.sa.wet[year.index],
                                                                                                                  survH = vr.value.perturbed.surv.h.wet,
                                                                                                                  survD = annual.surv.d.wet[year.index], 
                                                                                                                  transHD = annual.trans.hd.wet[year.index],
                                                                                                                  emigH = annual.emig.h.wet[year.index],
                                                                                                                  recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.h.wet[, , i] = vr.value.perturbed.annual.mat.surv.h.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Dominant dry-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.d.dry = annual.surv.d.dry[year.index] + mean.surv.d.dry
    
    mean.perturbed.annual.mat.surv.d.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                    survS = annual.surv.sa.dry[year.index],
                                                    survH = annual.surv.h.dry[year.index], 
                                                    survD = mean.perturbed.surv.d.dry,
                                                    transHD = annual.trans.hd.dry[year.index], 
                                                    emigH = annual.emig.h.dry[year.index], 
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = annual.surv.d.wet[year.index], 
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index], 
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.d.dry[, , i] = mean.perturbed.annual.mat.surv.d.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.d.dry = annual.surv.d.dry[year.index] + (annual.surv.d.dry[year.index] - mean.surv.d.dry)
    
    sd.perturbed.annual.mat.surv.d.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                  survS = annual.surv.sa.dry[year.index],
                                                  survH = annual.surv.h.dry[year.index],
                                                  survD = sd.perturbed.surv.d.dry,
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = annual.emig.h.dry[year.index], 
                                                  recruitH = annual.recruit.h.dry[year.index],
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = annual.surv.h.wet[year.index], 
                                                                                                            survD = annual.surv.d.wet[year.index], 
                                                                                                            transHD = annual.trans.hd.wet[year.index],
                                                                                                            emigH = annual.emig.h.wet[year.index], 
                                                                                                            recruitH = annual.recruit.h.wet[year.index],
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.d.dry[, , i] = sd.perturbed.annual.mat.surv.d.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.d.dry = annual.surv.d.dry[year.index] + annual.surv.d.dry[year.index]
    
    vr.value.perturbed.annual.mat.surv.d.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                        survS = annual.surv.sa.dry[year.index],
                                                        survH = annual.surv.h.dry[year.index],
                                                        survD = vr.value.perturbed.surv.d.dry,
                                                        transHD = annual.trans.hd.dry[year.index], 
                                                        emigH = annual.emig.h.dry[year.index],
                                                        recruitH = annual.recruit.h.dry[year.index],
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                  survS = annual.surv.sa.wet[year.index], 
                                                                                                                  survH = annual.surv.h.wet[year.index],
                                                                                                                  survD = annual.surv.d.wet[year.index],
                                                                                                                  transHD = annual.trans.hd.wet[year.index], 
                                                                                                                  emigH = annual.emig.h.wet[year.index],
                                                                                                                  recruitH = annual.recruit.h.wet[year.index],
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.d.dry[, , i] = vr.value.perturbed.annual.mat.surv.d.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Dominant wet-season survival
    
    # Perturbation with the mean
    mean.perturbed.surv.d.wet = annual.surv.d.wet[year.index] + mean.surv.d.wet
    
    mean.perturbed.annual.mat.surv.d.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                    survS = annual.surv.sa.dry[year.index], 
                                                    survH = annual.surv.h.dry[year.index], 
                                                    survD = annual.surv.d.dry[year.index], 
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = annual.emig.h.dry[year.index],
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                              survS = annual.surv.sa.wet[year.index], 
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = mean.perturbed.surv.d.wet, 
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index], 
                                                                                                              recruitH = annual.recruit.h.wet[year.index], 
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.surv.d.wet[, , i] = mean.perturbed.annual.mat.surv.d.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.surv.d.wet = annual.surv.d.wet[year.index] + (annual.surv.d.wet[year.index] - mean.surv.d.wet)
    
    sd.perturbed.annual.mat.surv.d.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                  survS = annual.surv.sa.dry[year.index], 
                                                  survH = annual.surv.h.dry[year.index], 
                                                  survD = annual.surv.d.dry[year.index],
                                                  transHD = annual.trans.hd.dry[year.index], 
                                                  emigH = annual.emig.h.dry[year.index], 
                                                  recruitH = annual.recruit.h.dry[year.index], 
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                            survS = annual.surv.sa.wet[year.index], 
                                                                                                            survH = annual.surv.h.wet[year.index],
                                                                                                            survD = sd.perturbed.surv.d.wet, 
                                                                                                            transHD = annual.trans.hd.wet[year.index], 
                                                                                                            emigH = annual.emig.h.wet[year.index], 
                                                                                                            recruitH = annual.recruit.h.wet[year.index], 
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.surv.d.wet[,  , i] = sd.perturbed.annual.mat.surv.d.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.surv.d.wet = annual.surv.d.wet[year.index] + annual.surv.d.wet[year.index]
    
    vr.value.perturbed.annual.mat.surv.d.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                        survS = annual.surv.sa.dry[year.index], 
                                                        survH = annual.surv.h.dry[year.index],
                                                        survD = annual.surv.d.dry[year.index],
                                                        transHD = annual.trans.hd.dry[year.index], 
                                                        emigH = annual.emig.h.dry[year.index], 
                                                        recruitH = annual.recruit.h.dry[year.index],
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                  survS = annual.surv.sa.wet[year.index], 
                                                                                                                  survH = annual.surv.h.wet[year.index], 
                                                                                                                  survD = vr.value.perturbed.surv.d.wet, 
                                                                                                                  transHD = annual.trans.hd.wet[year.index],
                                                                                                                  emigH = annual.emig.h.wet[year.index], 
                                                                                                                  recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.surv.d.wet[, , i] = vr.value.perturbed.annual.mat.surv.d.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper-Dominant dry-season transition
    
    # Perturbation with the mean
    mean.perturbed.trans.hd.dry = annual.trans.hd.dry[year.index] + mean.trans.hd.dry
    
    mean.perturbed.annual.mat.trans.hd.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                      survS = annual.surv.sa.dry[year.index], 
                                                      survH = annual.surv.h.dry[year.index], 
                                                      survD = annual.surv.d.dry[year.index], 
                                                      transHD = mean.perturbed.trans.hd.dry, 
                                                      emigH = annual.emig.h.dry[year.index], 
                                                      recruitH = annual.recruit.h.dry[year.index], 
                                                      recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                survS = annual.surv.sa.wet[year.index], 
                                                                                                                survH = annual.surv.h.wet[year.index], 
                                                                                                                survD = annual.surv.d.wet[year.index], 
                                                                                                                transHD = annual.trans.hd.wet[year.index], 
                                                                                                                emigH = annual.emig.h.wet[year.index], 
                                                                                                                recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.trans.hd.dry[, , i] = mean.perturbed.annual.mat.trans.hd.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.trans.hd.dry = annual.trans.hd.dry[year.index] + (annual.trans.hd.dry[year.index] - mean.trans.hd.dry)
    
    sd.perturbed.annual.mat.trans.hd.dry = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                    survS = annual.surv.sa.dry[year.index], 
                                                    survH = annual.surv.h.dry[year.index],
                                                    survD = annual.surv.d.dry[year.index], 
                                                    transHD = sd.perturbed.trans.hd.dry, 
                                                    emigH = annual.emig.h.dry[year.index], 
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = annual.surv.d.wet[year.index], 
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index], 
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.trans.hd.dry[, , i] = sd.perturbed.annual.mat.trans.hd.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.trans.hd.dry = annual.trans.hd.dry[year.index] + annual.trans.hd.dry[year.index]
    
    vr.value.perturbed.annual.mat.trans.hd.dry = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                          survS = annual.surv.sa.dry[year.index], 
                                                          survH = annual.surv.h.dry[year.index], 
                                                          survD = annual.surv.d.dry[year.index],
                                                          transHD = vr.value.perturbed.trans.hd.dry,
                                                          emigH = annual.emig.h.dry[year.index],
                                                          recruitH = annual.recruit.h.dry[year.index], 
                                                          recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                    survS = annual.surv.sa.wet[year.index],
                                                                                                                    survH = annual.surv.h.wet[year.index],
                                                                                                                    survD = annual.surv.d.wet[year.index], 
                                                                                                                    transHD = annual.trans.hd.wet[year.index], 
                                                                                                                    emigH = annual.emig.h.wet[year.index], 
                                                                                                                    recruitH = annual.recruit.h.wet[year.index],
                                                                                                                    recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.trans.hd.dry[, , i] = vr.value.perturbed.annual.mat.trans.hd.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper-Dominant wet-season transition
    
    # Perturbation with the mean
    mean.perturbed.trans.hd.wet = annual.trans.hd.wet[year.index] + mean.trans.hd.wet
    
    mean.perturbed.annual.mat.trans.hd.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                      survS = annual.surv.sa.dry[year.index],
                                                      survH = annual.surv.h.dry[year.index], 
                                                      survD = annual.surv.d.dry[year.index], 
                                                      transHD = annual.trans.hd.dry[year.index],
                                                      emigH = annual.emig.h.dry[year.index], 
                                                      recruitH = annual.recruit.h.dry[year.index], 
                                                      recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                survS = annual.surv.sa.wet[year.index], 
                                                                                                                survH = annual.surv.h.wet[year.index],
                                                                                                                survD = annual.surv.d.wet[year.index],
                                                                                                                transHD = mean.perturbed.trans.hd.wet, 
                                                                                                                emigH = annual.emig.h.wet[year.index], 
                                                                                                                recruitH = annual.recruit.h.wet[year.index],
                                                                                                                recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.trans.hd.wet[, , i] = mean.perturbed.annual.mat.trans.hd.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.trans.hd.wet = annual.trans.hd.wet[year.index] + (annual.trans.hd.wet[year.index] - mean.trans.hd.wet)
    
    sd.perturbed.annual.mat.trans.hd.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                    survS = annual.surv.sa.dry[year.index],
                                                    survH = annual.surv.h.dry[year.index],
                                                    survD = annual.surv.d.dry[year.index], 
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = annual.emig.h.dry[year.index], 
                                                    recruitH = annual.recruit.h.dry[year.index], 
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                              survS = annual.surv.sa.wet[year.index], 
                                                                                                              survH = annual.surv.h.wet[year.index], 
                                                                                                              survD = annual.surv.d.wet[year.index],
                                                                                                              transHD = sd.perturbed.trans.hd.wet, 
                                                                                                              emigH = annual.emig.h.wet[year.index],
                                                                                                              recruitH = annual.recruit.h.wet[year.index], 
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.trans.hd.wet[, , i] = sd.perturbed.annual.mat.trans.hd.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.trans.hd.wet = annual.trans.hd.wet[year.index] + annual.trans.hd.wet[year.index]
    
    vr.value.perturbed.annual.mat.trans.hd.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                          survS = annual.surv.sa.dry[year.index],
                                                          survH = annual.surv.h.dry[year.index], 
                                                          survD = annual.surv.d.dry[year.index], 
                                                          transHD = annual.trans.hd.dry[year.index], 
                                                          emigH = annual.emig.h.dry[year.index],
                                                          recruitH = annual.recruit.h.dry[year.index],
                                                          recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                    survS = annual.surv.sa.wet[year.index], 
                                                                                                                    survH = annual.surv.h.wet[year.index], 
                                                                                                                    survD = annual.surv.d.wet[year.index], 
                                                                                                                    transHD = vr.value.perturbed.trans.hd.wet, 
                                                                                                                    emigH = annual.emig.h.wet[year.index], 
                                                                                                                    recruitH = annual.recruit.h.wet[year.index],
                                                                                                                    recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.trans.hd.wet[, , i] = vr.value.perturbed.annual.mat.trans.hd.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper dry-season emigration
    
    # Perturbation with the mean
    mean.perturbed.emig.h.dry = annual.emig.h.dry[year.index] + mean.emig.h.dry
    
    mean.perturbed.annual.mat.emig.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                    survS = annual.surv.sa.dry[year.index], 
                                                    survH = annual.surv.h.dry[year.index], 
                                                    survD = annual.surv.d.dry[year.index], 
                                                    transHD = annual.trans.hd.dry[year.index],
                                                    emigH = mean.perturbed.emig.h.dry, 
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                              survS = annual.surv.sa.wet[year.index],
                                                                                                              survH = annual.surv.h.wet[year.index],
                                                                                                              survD = annual.surv.d.wet[year.index],
                                                                                                              transHD = annual.trans.hd.wet[year.index],
                                                                                                              emigH = annual.emig.h.wet[year.index], 
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.emig.h.dry[, , i] = mean.perturbed.annual.mat.emig.h.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.emig.h.dry = annual.emig.h.dry[year.index] + (annual.emig.h.dry[year.index] - mean.emig.h.dry)
    
    sd.perturbed.annual.mat.emig.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                  survS = annual.surv.sa.dry[year.index],
                                                  survH = annual.surv.h.dry[year.index],
                                                  survD = annual.surv.d.dry[year.index],
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = sd.perturbed.emig.h.dry, 
                                                  recruitH = annual.recruit.h.dry[year.index],
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = annual.surv.h.wet[year.index], 
                                                                                                            survD = annual.surv.d.wet[year.index],
                                                                                                            transHD = annual.trans.hd.wet[year.index],
                                                                                                            emigH = annual.emig.h.wet[year.index], 
                                                                                                            recruitH = annual.recruit.h.wet[year.index],
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.emig.h.dry[, , i] = sd.perturbed.annual.mat.emig.h.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.emig.h.dry = annual.emig.h.dry[year.index] + annual.emig.h.dry[year.index]
    
    vr.value.perturbed.annual.mat.emig.h.dry = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                        survS = annual.surv.sa.dry[year.index],
                                                        survH = annual.surv.h.dry[year.index], 
                                                        survD = annual.surv.d.dry[year.index], 
                                                        transHD = annual.trans.hd.dry[year.index], 
                                                        emigH = vr.value.perturbed.emig.h.dry, 
                                                        recruitH = annual.recruit.h.dry[year.index],
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                  survS = annual.surv.sa.wet[year.index], 
                                                                                                                  survH = annual.surv.h.wet[year.index], 
                                                                                                                  survD = annual.surv.d.wet[year.index],
                                                                                                                  transHD = annual.trans.hd.wet[year.index], 
                                                                                                                  emigH = annual.emig.h.wet[year.index], 
                                                                                                                  recruitH = annual.recruit.h.wet[year.index],
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.emig.h.dry[, , i] = vr.value.perturbed.annual.mat.emig.h.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper wet-season emigration
    
    # Perturbation with the mean
    mean.perturbed.emig.h.wet = annual.emig.h.wet[year.index] + mean.emig.h.wet
    
    mean.perturbed.annual.mat.emig.h.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                    survS = annual.surv.sa.dry[year.index], 
                                                    survH = annual.surv.h.dry[year.index],
                                                    survD = annual.surv.d.dry[year.index],
                                                    transHD = annual.trans.hd.dry[year.index], 
                                                    emigH = annual.emig.h.dry[year.index], 
                                                    recruitH = annual.recruit.h.dry[year.index],
                                                    recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                              survS = annual.surv.sa.wet[year.index], 
                                                                                                              survH = annual.surv.h.wet[year.index], 
                                                                                                              survD = annual.surv.d.wet[year.index], 
                                                                                                              transHD = annual.trans.hd.wet[year.index], 
                                                                                                              emigH = mean.perturbed.emig.h.wet, 
                                                                                                              recruitH = annual.recruit.h.wet[year.index],
                                                                                                              recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.emig.h.wet[, , i] = mean.perturbed.annual.mat.emig.h.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.emig.h.wet = annual.emig.h.wet[year.index] + (annual.emig.h.wet[year.index] - mean.emig.h.wet)
    
    sd.perturbed.annual.mat.emig.h.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                  survS = annual.surv.sa.dry[year.index], 
                                                  survH = annual.surv.h.dry[year.index], 
                                                  survD = annual.surv.d.dry[year.index], 
                                                  transHD = annual.trans.hd.dry[year.index],
                                                  emigH = annual.emig.h.dry[year.index], 
                                                  recruitH = annual.recruit.h.dry[year.index], 
                                                  recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                            survS = annual.surv.sa.wet[year.index],
                                                                                                            survH = annual.surv.h.wet[year.index],
                                                                                                            survD = annual.surv.d.wet[year.index], 
                                                                                                            transHD = annual.trans.hd.wet[year.index],
                                                                                                            emigH = sd.perturbed.emig.h.wet, 
                                                                                                            recruitH = annual.recruit.h.wet[year.index], 
                                                                                                            recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.emig.h.wet[, , i] = sd.perturbed.annual.mat.emig.h.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.emig.h.wet = annual.emig.h.wet[year.index] + annual.emig.h.wet[year.index]
    
    vr.value.perturbed.annual.mat.emig.h.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                        survS = annual.surv.sa.dry[year.index], 
                                                        survH = annual.surv.h.dry[year.index],
                                                        survD = annual.surv.d.dry[year.index],
                                                        transHD = annual.trans.hd.dry[year.index],
                                                        emigH = annual.emig.h.dry[year.index], 
                                                        recruitH = annual.recruit.h.dry[year.index], 
                                                        recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                  survS = annual.surv.sa.wet[year.index],
                                                                                                                  survH = annual.surv.h.wet[year.index], 
                                                                                                                  survD = annual.surv.d.wet[year.index], 
                                                                                                                  transHD = annual.trans.hd.wet[year.index], 
                                                                                                                  emigH = vr.value.perturbed.emig.h.wet, 
                                                                                                                  recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                  recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.emig.h.wet[, , i] = vr.value.perturbed.annual.mat.emig.h.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper dry-season recruitment
    
    # Perturbation with the mean
    mean.perturbed.recruit.h.dry = annual.recruit.h.dry[year.index] + mean.recruit.h.dry
    
    mean.perturbed.annual.mat.recruit.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                       survS = annual.surv.sa.dry[year.index],
                                                       survH = annual.surv.h.dry[year.index], 
                                                       survD = annual.surv.d.dry[year.index], 
                                                       transHD = annual.trans.hd.dry[year.index],
                                                       emigH = annual.emig.h.dry[year.index],
                                                       recruitH = mean.perturbed.recruit.h.dry,
                                                       recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                 survS = annual.surv.sa.wet[year.index], 
                                                                                                                 survH = annual.surv.h.wet[year.index], 
                                                                                                                 survD = annual.surv.d.wet[year.index], 
                                                                                                                 transHD = annual.trans.hd.wet[year.index], 
                                                                                                                 emigH = annual.emig.h.wet[year.index],
                                                                                                                 recruitH = annual.recruit.h.wet[year.index],
                                                                                                                 recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.recruit.h.dry[, , i] = mean.perturbed.annual.mat.recruit.h.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.recruit.h.dry = annual.recruit.h.dry[year.index] + (annual.recruit.h.dry[year.index] - mean.recruit.h.dry)
    
    sd.perturbed.annual.mat.recruit.h.dry = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                     survS = annual.surv.sa.dry[year.index], 
                                                     survH = annual.surv.h.dry[year.index], 
                                                     survD = annual.surv.d.dry[year.index],
                                                     transHD = annual.trans.hd.dry[year.index],
                                                     emigH = annual.emig.h.dry[year.index], 
                                                     recruitH = sd.perturbed.recruit.h.dry, 
                                                     recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                               survS = annual.surv.sa.wet[year.index], 
                                                                                                               survH = annual.surv.h.wet[year.index], 
                                                                                                               survD = annual.surv.d.wet[year.index], 
                                                                                                               transHD = annual.trans.hd.wet[year.index],
                                                                                                               emigH = annual.emig.h.wet[year.index],
                                                                                                               recruitH = annual.recruit.h.wet[year.index],
                                                                                                               recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.recruit.h.dry[, , i] = sd.perturbed.annual.mat.recruit.h.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.recruit.h.dry = annual.recruit.h.dry[year.index] + annual.recruit.h.dry[year.index]
    
    vr.value.perturbed.annual.mat.recruit.h.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                           survS = annual.surv.sa.dry[year.index],
                                                           survH = annual.surv.h.dry[year.index], 
                                                           survD = annual.surv.d.dry[year.index], 
                                                           transHD = annual.trans.hd.dry[year.index],
                                                           emigH = annual.emig.h.dry[year.index], 
                                                           recruitH = vr.value.perturbed.recruit.h.dry,
                                                           recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                     survS = annual.surv.sa.wet[year.index],
                                                                                                                     survH = annual.surv.h.wet[year.index], 
                                                                                                                     survD = annual.surv.d.wet[year.index], 
                                                                                                                     transHD = annual.trans.hd.wet[year.index], 
                                                                                                                     emigH = annual.emig.h.wet[year.index], 
                                                                                                                     recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                     recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.recruit.h.dry[, , i] = vr.value.perturbed.annual.mat.recruit.h.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Helper wet-season recruitment
    
    # Perturbation with the mean
    mean.perturbed.recruit.h.wet = annual.recruit.h.wet[year.index] + mean.recruit.h.wet
    
    mean.perturbed.annual.mat.recruit.h.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                       survS = annual.surv.sa.dry[year.index], 
                                                       survH = annual.surv.h.dry[year.index],
                                                       survD = annual.surv.d.dry[year.index],
                                                       transHD = annual.trans.hd.dry[year.index], 
                                                       emigH = annual.emig.h.dry[year.index], 
                                                       recruitH = annual.recruit.h.dry[year.index], 
                                                       recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                 survS = annual.surv.sa.wet[year.index],
                                                                                                                 survH = annual.surv.h.wet[year.index],
                                                                                                                 survD = annual.surv.d.wet[year.index],
                                                                                                                 transHD = annual.trans.hd.wet[year.index],
                                                                                                                 emigH = annual.emig.h.wet[year.index], 
                                                                                                                 recruitH = mean.perturbed.recruit.h.wet, 
                                                                                                                 recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.recruit.h.wet[, , i] = mean.perturbed.annual.mat.recruit.h.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.recruit.h.wet = annual.recruit.h.wet[year.index] + (annual.recruit.h.wet[year.index] - mean.recruit.h.wet)
    
    sd.perturbed.annual.mat.recruit.h.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                     survS = annual.surv.sa.dry[year.index],
                                                     survH = annual.surv.h.dry[year.index],
                                                     survD = annual.surv.d.dry[year.index],
                                                     transHD = annual.trans.hd.dry[year.index],
                                                     emigH = annual.emig.h.dry[year.index], 
                                                     recruitH = annual.recruit.h.dry[year.index], 
                                                     recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                               survS = annual.surv.sa.wet[year.index], 
                                                                                                               survH = annual.surv.h.wet[year.index], 
                                                                                                               survD = annual.surv.d.wet[year.index], 
                                                                                                               transHD = annual.trans.hd.wet[year.index],
                                                                                                               emigH = annual.emig.h.wet[year.index], 
                                                                                                               recruitH = sd.perturbed.recruit.h.wet, 
                                                                                                               recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.recruit.h.wet[, , i] = sd.perturbed.annual.mat.recruit.h.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.recruit.h.wet = annual.recruit.h.wet[year.index] + annual.recruit.h.wet[year.index]
    
    vr.value.perturbed.annual.mat.recruit.h.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                           survS = annual.surv.sa.dry[year.index],
                                                           survH = annual.surv.h.dry[year.index],
                                                           survD = annual.surv.d.dry[year.index], 
                                                           transHD = annual.trans.hd.dry[year.index],
                                                           emigH = annual.emig.h.dry[year.index],
                                                           recruitH = annual.recruit.h.dry[year.index], 
                                                           recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                     survS = annual.surv.sa.wet[year.index],
                                                                                                                     survH = annual.surv.h.wet[year.index], 
                                                                                                                     survD = annual.surv.d.wet[year.index], 
                                                                                                                     transHD = annual.trans.hd.wet[year.index],
                                                                                                                     emigH = annual.emig.h.wet[year.index],
                                                                                                                     recruitH = vr.value.perturbed.recruit.h.wet,
                                                                                                                     recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.recruit.h.wet[, , i] = vr.value.perturbed.annual.mat.recruit.h.wet - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Dominant dry-season recruitment
    
    # Perturbation with the mean
    mean.perturbed.recruit.d.dry = annual.recruit.d.dry[year.index] + mean.recruit.d.dry
    
    mean.perturbed.annual.mat.recruit.d.dry = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                       survS = annual.surv.sa.dry[year.index], 
                                                       survH = annual.surv.h.dry[year.index], 
                                                       survD = annual.surv.d.dry[year.index], 
                                                       transHD = annual.trans.hd.dry[year.index],
                                                       emigH = annual.emig.h.dry[year.index], 
                                                       recruitH = annual.recruit.h.dry[year.index], 
                                                       recruitD = mean.perturbed.recruit.d.dry) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                             survS = annual.surv.sa.wet[year.index], 
                                                                                                             survH = annual.surv.h.wet[year.index], 
                                                                                                             survD = annual.surv.d.wet[year.index],
                                                                                                             transHD = annual.trans.hd.wet[year.index], 
                                                                                                             emigH = annual.emig.h.wet[year.index], 
                                                                                                             recruitH = annual.recruit.h.wet[year.index], 
                                                                                                             recruitD = annual.recruit.d.wet[year.index])
    
    mean.C.mat.recruit.d.dry[, , i] = mean.perturbed.annual.mat.recruit.d.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.recruit.d.dry = annual.recruit.d.dry[year.index] + (annual.recruit.d.dry[year.index] - mean.recruit.d.dry)
    
    sd.perturbed.annual.mat.recruit.d.dry = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                     survS = annual.surv.sa.dry[year.index], 
                                                     survH = annual.surv.h.dry[year.index], 
                                                     survD = annual.surv.d.dry[year.index], 
                                                     transHD = annual.trans.hd.dry[year.index],
                                                     emigH = annual.emig.h.dry[year.index],
                                                     recruitH = annual.recruit.h.dry[year.index],
                                                     recruitD = sd.perturbed.recruit.d.dry) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                         survS = annual.surv.sa.wet[year.index],
                                                                                                         survH = annual.surv.h.wet[year.index], 
                                                                                                         survD = annual.surv.d.wet[year.index], 
                                                                                                         transHD = annual.trans.hd.wet[year.index],
                                                                                                         emigH = annual.emig.h.wet[year.index], 
                                                                                                         recruitH = annual.recruit.h.wet[year.index], 
                                                                                                         recruitD = annual.recruit.d.wet[year.index])
    
    sd.C.mat.recruit.d.dry[, , i] = sd.perturbed.annual.mat.recruit.d.dry - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.recruit.d.dry = annual.recruit.d.dry[year.index] + annual.recruit.d.dry[year.index]
    
    vr.value.perturbed.annual.mat.recruit.d.dry = buildmat(survJ = annual.surv.j.dry[year.index],
                                                           survS = annual.surv.sa.dry[year.index],
                                                           survH = annual.surv.h.dry[year.index], 
                                                           survD = annual.surv.d.dry[year.index], 
                                                           transHD = annual.trans.hd.dry[year.index],
                                                           emigH = annual.emig.h.dry[year.index], 
                                                           recruitH = annual.recruit.h.dry[year.index], 
                                                           recruitD = vr.value.perturbed.recruit.d.dry) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                     survS = annual.surv.sa.wet[year.index],
                                                                                                                     survH = annual.surv.h.wet[year.index], 
                                                                                                                     survD = annual.surv.d.wet[year.index], 
                                                                                                                     transHD = annual.trans.hd.wet[year.index], 
                                                                                                                     emigH = annual.emig.h.wet[year.index], 
                                                                                                                     recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                     recruitD = annual.recruit.d.wet[year.index])
    
    vr.value.C.mat.recruit.d.dry[, , i] = vr.value.perturbed.annual.mat.recruit.d.dry - annual.mat # C matrix with perturbed elements
    
    ################################
    
    # Dominant wet-season recruitment
    
    # Perturbation with the mean
    mean.perturbed.recruit.d.wet = annual.recruit.d.wet[year.index] + mean.recruit.d.wet
    
    mean.perturbed.annual.mat.recruit.d.wet = buildmat(survJ = annual.surv.j.dry[year.index], 
                                                       survS = annual.surv.sa.dry[year.index],
                                                       survH = annual.surv.h.dry[year.index], 
                                                       survD = annual.surv.d.dry[year.index], 
                                                       transHD = annual.trans.hd.dry[year.index], 
                                                       emigH = annual.emig.h.dry[year.index],
                                                       recruitH = annual.recruit.h.dry[year.index],
                                                       recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index],
                                                                                                                 survS = annual.surv.sa.wet[year.index], 
                                                                                                                 survH = annual.surv.h.wet[year.index],
                                                                                                                 survD = annual.surv.d.wet[year.index], 
                                                                                                                 transHD = annual.trans.hd.wet[year.index],
                                                                                                                 emigH = annual.emig.h.wet[year.index], 
                                                                                                                 recruitH = annual.recruit.h.wet[year.index], 
                                                                                                                 recruitD = mean.perturbed.recruit.d.wet)
    
    mean.C.mat.recruit.d.wet[, , i] = mean.perturbed.annual.mat.recruit.d.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the standard deviation
    sd.perturbed.recruit.d.wet = annual.recruit.d.wet[year.index] + (annual.recruit.d.wet[year.index] - mean.recruit.d.wet)
    
    sd.perturbed.annual.mat.recruit.d.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                     survS = annual.surv.sa.dry[year.index],
                                                     survH = annual.surv.h.dry[year.index], 
                                                     survD = annual.surv.d.dry[year.index], 
                                                     transHD = annual.trans.hd.dry[year.index],
                                                     emigH = annual.emig.h.dry[year.index], 
                                                     recruitH = annual.recruit.h.dry[year.index],
                                                     recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                               survS = annual.surv.sa.wet[year.index], 
                                                                                                               survH = annual.surv.h.wet[year.index], 
                                                                                                               survD = annual.surv.d.wet[year.index], 
                                                                                                               transHD = annual.trans.hd.wet[year.index],
                                                                                                               emigH = annual.emig.h.wet[year.index], 
                                                                                                               recruitH = annual.recruit.h.wet[year.index], 
                                                                                                               recruitD = sd.perturbed.recruit.d.wet)
    
    sd.C.mat.recruit.d.wet[, , i] = sd.perturbed.annual.mat.recruit.d.wet - annual.mat # C matrix with perturbed elements
    
    # Perturbation with the vital rate value
    vr.value.perturbed.recruit.d.wet = annual.recruit.d.wet[year.index] + annual.recruit.d.wet[year.index]
    
    vr.value.perturbed.annual.mat.recruit.d.wet = buildmat(survJ = annual.surv.j.dry[year.index],
                                                           survS = annual.surv.sa.dry[year.index],
                                                           survH = annual.surv.h.dry[year.index], 
                                                           survD = annual.surv.d.dry[year.index], 
                                                           transHD = annual.trans.hd.dry[year.index], 
                                                           emigH = annual.emig.h.dry[year.index],
                                                           recruitH = annual.recruit.h.dry[year.index], 
                                                           recruitD = annual.recruit.d.dry[year.index]) %*% buildmat(survJ = annual.surv.j.wet[year.index], 
                                                                                                                     survS = annual.surv.sa.wet[year.index], 
                                                                                                                     survH = annual.surv.h.wet[year.index],
                                                                                                                     survD = annual.surv.d.wet[year.index],
                                                                                                                     transHD = annual.trans.hd.wet[year.index],
                                                                                                                     emigH = annual.emig.h.wet[year.index], 
                                                                                                                     recruitH = annual.recruit.h.wet[year.index],
                                                                                                                     recruitD = vr.value.perturbed.recruit.d.wet)
    
    vr.value.C.mat.recruit.d.wet[, , i] = vr.value.perturbed.annual.mat.recruit.d.wet - annual.mat # C matrix with perturbed elements
    
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
    year.index = which(unique(data.meerkats$year) == year)
    
    # Create annual matrix
    annual.mat = dry.matrices[, , year.index] %*% wet.matrices[, , year.index]
    
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
    
    # Elasticities juvenile dry-season survival
    mean.pert.surv.j.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.j.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.j.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.j.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.j.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.j.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.j.dry[, , x] <- elasticity.mean.surv.j.dry[, , x] + mean.pert.surv.j.dry
    elasticity.sd.surv.j.dry[, , x] <- elasticity.sd.surv.j.dry[, , x] + sd.pert.surv.j.dry
    elasticity.stochastic.surv.j.dry[, , x] <- elasticity.stochastic.surv.j.dry[, , x] + vr.value.pert.surv.j.dry
    
    
    # Elasticities juvenile wet-season survival
    mean.pert.surv.j.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.j.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.j.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.j.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.j.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.j.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.j.wet[, , x] <- elasticity.mean.surv.j.wet[, , x] + mean.pert.surv.j.wet
    elasticity.sd.surv.j.wet[, , x] <- elasticity.sd.surv.j.wet[, , x] + sd.pert.surv.j.wet
    elasticity.stochastic.surv.j.wet[, , x] <- elasticity.stochastic.surv.j.wet[, , x] + vr.value.pert.surv.j.wet
    
    
    # Elasticities subadult dry-season survival
    mean.pert.surv.sa.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.sa.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.sa.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.sa.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.sa.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.sa.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.sa.dry[, , x] <- elasticity.mean.surv.sa.dry[, , x] + mean.pert.surv.sa.dry
    elasticity.sd.surv.sa.dry[, , x] <- elasticity.sd.surv.sa.dry[, , x] + sd.pert.surv.sa.dry
    elasticity.stochastic.surv.sa.dry[, , x] <- elasticity.stochastic.surv.sa.dry[, , x] + vr.value.pert.surv.sa.dry
    
    
    # Elasticities subadult wet-season survival
    mean.pert.surv.sa.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.sa.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.sa.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.sa.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.sa.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.sa.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.sa.wet[, , x] <- elasticity.mean.surv.sa.wet[, , x] + mean.pert.surv.sa.wet
    elasticity.sd.surv.sa.wet[, , x] <- elasticity.sd.surv.sa.wet[, , x] + sd.pert.surv.sa.wet
    elasticity.stochastic.surv.sa.wet[, , x] <- elasticity.stochastic.surv.sa.wet[, , x] + vr.value.pert.surv.sa.wet
    
    
    # Elasticities helper dry-season survival
    mean.pert.surv.h.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.h.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.h.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.h.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.h.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.h.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.h.dry[, , x] <- elasticity.mean.surv.h.dry[, , x] + mean.pert.surv.h.dry
    elasticity.sd.surv.h.dry[, , x] <- elasticity.sd.surv.h.dry[, , x] + sd.pert.surv.h.dry
    elasticity.stochastic.surv.h.dry[, , x] <- elasticity.stochastic.surv.h.dry[, , x] + vr.value.pert.surv.h.dry
    
    
    # Elasticities helper wet-season survival
    mean.pert.surv.h.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.h.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.h.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.h.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.h.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.h.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.h.wet[, , x] <- elasticity.mean.surv.h.wet[, , x] + mean.pert.surv.h.wet
    elasticity.sd.surv.h.wet[, , x] <- elasticity.sd.surv.h.wet[, , x] + sd.pert.surv.h.wet
    elasticity.stochastic.surv.h.wet[, , x] <- elasticity.stochastic.surv.h.wet[, , x] + vr.value.pert.surv.h.wet
    
    
    # Elasticities dominant dry-season survival
    mean.pert.surv.d.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.d.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.d.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.d.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.d.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.d.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.d.dry[, , x] <- elasticity.mean.surv.d.dry[, , x] + mean.pert.surv.d.dry
    elasticity.sd.surv.d.dry[, , x] <- elasticity.sd.surv.d.dry[, , x] + sd.pert.surv.d.dry
    elasticity.stochastic.surv.d.dry[, , x] <- elasticity.stochastic.surv.d.dry[, , x] + vr.value.pert.surv.d.dry
    
    
    # Elasticities dominant wet-season survival
    mean.pert.surv.d.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.surv.d.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.surv.d.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.surv.d.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.surv.d.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.surv.d.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.surv.d.wet[, , x] <- elasticity.mean.surv.d.wet[, , x] + mean.pert.surv.d.wet
    elasticity.sd.surv.d.wet[, , x] <- elasticity.sd.surv.d.wet[, , x] + sd.pert.surv.d.wet
    elasticity.stochastic.surv.d.wet[, , x] <- elasticity.stochastic.surv.d.wet[, , x] + vr.value.pert.surv.d.wet
    
    
    # Elasticities helper-dominant dry-season transition
    mean.pert.trans.hd.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.trans.hd.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.trans.hd.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.trans.hd.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.trans.hd.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.trans.hd.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.trans.hd.dry[, , x] <- elasticity.mean.trans.hd.dry[, , x] + mean.pert.trans.hd.dry
    elasticity.sd.trans.hd.dry[, , x] <- elasticity.sd.trans.hd.dry[, , x] + sd.pert.trans.hd.dry
    elasticity.stochastic.trans.hd.dry[, , x] <- elasticity.stochastic.trans.hd.dry[, , x] + vr.value.pert.trans.hd.dry
    
    
    # Elasticities helper-dominant wet-season transition
    mean.pert.trans.hd.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.trans.hd.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.trans.hd.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.trans.hd.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.trans.hd.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.trans.hd.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.trans.hd.wet[, , x] <- elasticity.mean.trans.hd.wet[, , x] + mean.pert.trans.hd.wet
    elasticity.sd.trans.hd.wet[, , x] <- elasticity.sd.trans.hd.wet[, , x] + sd.pert.trans.hd.wet
    elasticity.stochastic.trans.hd.wet[, , x] <- elasticity.stochastic.trans.hd.wet[, , x] + vr.value.pert.trans.hd.wet
    
    
    # Elasticities helper dry-season emigration
    mean.pert.emig.h.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.emig.h.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.emig.h.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.emig.h.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.emig.h.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.emig.h.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.emig.h.dry[, , x] <- elasticity.mean.emig.h.dry[, , x] + mean.pert.emig.h.dry
    elasticity.sd.emig.h.dry[, , x] <- elasticity.sd.emig.h.dry[, , x] + sd.pert.emig.h.dry
    elasticity.stochastic.emig.h.dry[, , x] <- elasticity.stochastic.emig.h.dry[, , x] + vr.value.pert.emig.h.dry
    
    
    # Elasticities helper wet-season emigration
    mean.pert.emig.h.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.emig.h.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.emig.h.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.emig.h.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.emig.h.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.emig.h.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.emig.h.wet[, , x] <- elasticity.mean.emig.h.wet[, , x] + mean.pert.emig.h.wet
    elasticity.sd.emig.h.wet[, , x] <- elasticity.sd.emig.h.wet[, , x] + sd.pert.emig.h.wet
    elasticity.stochastic.emig.h.wet[, , x] <- elasticity.stochastic.emig.h.wet[, , x] + vr.value.pert.emig.h.wet
    
    
    # Elasticities helper dry-season recruitment
    mean.pert.recruit.h.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.recruit.h.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.recruit.h.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.recruit.h.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.recruit.h.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.recruit.h.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.recruit.h.dry[, , x] <- elasticity.mean.recruit.h.dry[, , x] + mean.pert.recruit.h.dry
    elasticity.sd.recruit.h.dry[, , x] <- elasticity.sd.recruit.h.dry[, , x] + sd.pert.recruit.h.dry
    elasticity.stochastic.recruit.h.dry[, , x] <- elasticity.stochastic.recruit.h.dry[, , x] + vr.value.pert.recruit.h.dry
    
    
    # Elasticities helper wet-season recruitment
    mean.pert.recruit.h.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.recruit.h.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.recruit.h.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.recruit.h.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.recruit.h.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.recruit.h.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.recruit.h.wet[, , x] <- elasticity.mean.recruit.h.wet[, , x] + mean.pert.recruit.h.wet
    elasticity.sd.recruit.h.wet[, , x] <- elasticity.sd.recruit.h.wet[, , x] + sd.pert.recruit.h.wet
    elasticity.stochastic.recruit.h.wet[, , x] <- elasticity.stochastic.recruit.h.wet[, , x] + vr.value.pert.recruit.h.wet
    
    
    # Elasticities dominant dry-season recruitment
    mean.pert.recruit.d.dry <- a * (diag(vvecs[, itime]) %*% mean.C.mat.recruit.d.dry[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.recruit.d.dry <- a * (diag(vvecs[, itime]) %*% sd.C.mat.recruit.d.dry[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.recruit.d.dry <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.recruit.d.dry[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.recruit.d.dry[, , x] <- elasticity.mean.recruit.d.dry[, , x] + mean.pert.recruit.d.dry
    elasticity.sd.recruit.d.dry[, , x] <- elasticity.sd.recruit.d.dry[, , x] + sd.pert.recruit.d.dry
    elasticity.stochastic.recruit.d.dry[, , x] <- elasticity.stochastic.recruit.d.dry[, , x] + vr.value.pert.recruit.d.dry
    
    
    # Elasticities dominant wet-season recruitment
    mean.pert.recruit.d.wet <- a * (diag(vvecs[, itime]) %*% mean.C.mat.recruit.d.wet[, , itime] %*% diag(uvecs[, i1]))
    sd.pert.recruit.d.wet <- a * (diag(vvecs[, itime]) %*% sd.C.mat.recruit.d.wet[, , itime] %*% diag(uvecs[, i1]))
    vr.value.pert.recruit.d.wet <- a * (diag(vvecs[, itime]) %*% vr.value.C.mat.recruit.d.wet[, , itime] %*% diag(uvecs[, i1]))
    
    elasticity.mean.recruit.d.wet[, , x] <- elasticity.mean.recruit.d.wet[, , x] + mean.pert.recruit.d.wet
    elasticity.sd.recruit.d.wet[, , x] <- elasticity.sd.recruit.d.wet[, , x] + sd.pert.recruit.d.wet
    elasticity.stochastic.recruit.d.wet[, , x] <- elasticity.stochastic.recruit.d.wet[, , x] + vr.value.pert.recruit.d.wet
    
  }
  
  elasticity.mean.surv.j.dry[, , x] = elasticity.mean.surv.j.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.j.dry[, , x] = elasticity.sd.surv.j.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.j.dry[, , x] = elasticity.stochastic.surv.j.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.surv.j.wet[, , x] = elasticity.mean.surv.j.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.j.wet[, , x] = elasticity.sd.surv.j.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.j.wet[, , x] = elasticity.stochastic.surv.j.wet[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.surv.sa.dry[, , x] = elasticity.mean.surv.sa.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.sa.dry[, , x] = elasticity.sd.surv.sa.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.sa.dry[, , x] = elasticity.stochastic.surv.sa.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.surv.sa.wet[, , x] = elasticity.mean.surv.sa.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.sa.wet[, , x] = elasticity.sd.surv.sa.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.sa.wet[, , x] = elasticity.stochastic.surv.sa.wet[, , x]/(sim.n.years.asymptotic)

  elasticity.mean.surv.h.dry[, , x] = elasticity.mean.surv.h.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.h.dry[, , x] = elasticity.sd.surv.h.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.h.dry[, , x] = elasticity.stochastic.surv.h.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.surv.h.wet[, , x] = elasticity.mean.surv.h.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.h.wet[, , x] = elasticity.sd.surv.h.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.h.wet[, , x] = elasticity.stochastic.surv.h.wet[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.surv.d.dry[, , x] = elasticity.mean.surv.d.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.d.dry[, , x] = elasticity.sd.surv.d.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.d.dry[, , x] = elasticity.stochastic.surv.d.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.surv.d.wet[, , x] = elasticity.mean.surv.d.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.surv.d.wet[, , x] = elasticity.sd.surv.d.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.surv.d.wet[, , x] = elasticity.stochastic.surv.d.wet[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.trans.hd.dry[, , x] = elasticity.mean.trans.hd.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.trans.hd.dry[, , x] = elasticity.sd.trans.hd.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.trans.hd.dry[, , x] = elasticity.stochastic.trans.hd.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.trans.hd.wet[, , x] = elasticity.mean.trans.hd.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.trans.hd.wet[, , x] = elasticity.sd.trans.hd.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.trans.hd.wet[, , x] = elasticity.stochastic.trans.hd.wet[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.emig.h.dry[, , x] = elasticity.mean.emig.h.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.emig.h.dry[, , x] = elasticity.sd.emig.h.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.emig.h.dry[, , x] = elasticity.stochastic.emig.h.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.emig.h.wet[, , x] = elasticity.mean.emig.h.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.emig.h.wet[, , x] = elasticity.sd.emig.h.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.emig.h.wet[, , x] = elasticity.stochastic.emig.h.wet[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.recruit.h.dry[, , x] = elasticity.mean.recruit.h.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.recruit.h.dry[, , x] = elasticity.sd.recruit.h.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.recruit.h.dry[, , x] = elasticity.stochastic.recruit.h.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.recruit.h.wet[, , x] = elasticity.mean.recruit.h.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.recruit.h.wet[, , x] = elasticity.sd.recruit.h.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.recruit.h.wet[, , x] = elasticity.stochastic.recruit.h.wet[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.recruit.d.dry[, , x] = elasticity.mean.recruit.d.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.recruit.d.dry[, , x] = elasticity.sd.recruit.d.dry[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.recruit.d.dry[, , x] = elasticity.stochastic.recruit.d.dry[, , x]/(sim.n.years.asymptotic)
  
  elasticity.mean.recruit.d.wet[, , x] = elasticity.mean.recruit.d.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.sd.recruit.d.wet[, , x] = elasticity.sd.recruit.d.wet[, , x]/(sim.n.years.asymptotic)
  elasticity.stochastic.recruit.d.wet[, , x] = elasticity.stochastic.recruit.d.wet[, , x]/(sim.n.years.asymptotic)
    
}




###########################################################################
#
# 4. Processing the results ----
#
###########################################################################

## 4.1. Check that the elasticities sum to 0 ----
# ------------------------------------------

elast.sums.vec = array(0, c(16, n.sim))

for(x in 1:n.sim){
  
  elast.sums.vec[1, x] = sum(round(elasticity.mean.surv.j.dry[, , x] + elasticity.sd.surv.j.dry[, , x] - elasticity.stochastic.surv.j.dry[, , x], 10))
  elast.sums.vec[2, x] = sum(round(elasticity.mean.surv.j.wet[, , x] + elasticity.sd.surv.j.wet[, , x] - elasticity.stochastic.surv.j.wet[, , x], 10))
  
  elast.sums.vec[3, x] = sum(round(elasticity.mean.surv.sa.dry[, , x] + elasticity.sd.surv.sa.dry[, , x] - elasticity.stochastic.surv.sa.dry[, , x], 10))
  elast.sums.vec[4, x] = sum(round(elasticity.mean.surv.sa.wet[, , x] + elasticity.sd.surv.sa.wet[, , x] - elasticity.stochastic.surv.sa.wet[, , x], 10))
  
  elast.sums.vec[5, x] = sum(round(elasticity.mean.surv.h.dry[, , x] + elasticity.sd.surv.h.dry[, , x] - elasticity.stochastic.surv.h.dry[, , x], 10))
  elast.sums.vec[6, x] = sum(round(elasticity.mean.surv.h.wet[, , x] + elasticity.sd.surv.h.wet[, , x] - elasticity.stochastic.surv.h.wet[, , x], 10))
  
  elast.sums.vec[7, x] = sum(round(elasticity.mean.surv.d.dry[, , x] + elasticity.sd.surv.d.dry[, , x] - elasticity.stochastic.surv.d.dry[, , x], 10))
  elast.sums.vec[8, x] = sum(round(elasticity.mean.surv.d.wet[, , x] + elasticity.sd.surv.d.wet[, , x] - elasticity.stochastic.surv.d.wet[, , x], 10))
  
  elast.sums.vec[9, x] = sum(round(elasticity.mean.trans.hd.dry[, , x] + elasticity.sd.trans.hd.dry[, , x] - elasticity.stochastic.trans.hd.dry[, , x], 10))
  elast.sums.vec[10, x] = sum(round(elasticity.mean.trans.hd.wet[, , x] + elasticity.sd.trans.hd.wet[, , x] - elasticity.stochastic.trans.hd.wet[, , x], 10))
  
  elast.sums.vec[11, x] = sum(round(elasticity.mean.emig.h.dry[, , x] + elasticity.sd.emig.h.dry[, , x] - elasticity.stochastic.emig.h.dry[, , x], 10))
  elast.sums.vec[12, x] = sum(round(elasticity.mean.emig.h.wet[, , x] + elasticity.sd.emig.h.wet[, , x] - elasticity.stochastic.emig.h.wet[, , x], 10))
  
  elast.sums.vec[13, x] = sum(round(elasticity.mean.recruit.h.dry[, , x] + elasticity.sd.recruit.h.dry[, , x] - elasticity.stochastic.recruit.h.dry[, , x], 10))
  elast.sums.vec[14, x] = sum(round(elasticity.mean.recruit.h.wet[, , x] + elasticity.sd.recruit.h.wet[, , x] - elasticity.stochastic.recruit.h.wet[, , x], 10))
  
  elast.sums.vec[15, x] = sum(round(elasticity.mean.recruit.d.dry[, , x] + elasticity.sd.recruit.d.dry[, , x] - elasticity.stochastic.recruit.d.dry[, , x], 10))
  elast.sums.vec[16, x] = sum(round(elasticity.mean.recruit.d.wet[, , x] + elasticity.sd.recruit.d.wet[, , x] - elasticity.stochastic.recruit.d.wet[, , x], 10))

}


# 4.2. Processing elasticities results ----
# ------------------------------------

## 4.2.1. Creating data frame with elasticities ----
# ---------------------------------------------

elast.df = expand.grid(vr = c("J survival (dry)", "J survival (wet)", "SA survival (dry)", "SA survival (wet)", "H survival (dry)", "H survival (wet)", "D survival (dry)", "D survival (wet)", "H-D transition (dry)", "H-D transition (wet)", "H emigration (dry)", "H emigration (wet)", "H recruitment (dry)", "H recruitment (wet)", "D recruitment (dry)", "D recruitment (wet)"), elast.type = c("Mean", "SD", "Stochastic"), elast.mean = NA, elast.low = NA, elast.up = NA)


## 4.2.2. Juvenile dry-season survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "J survival (dry)")] = c(mean(abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum))),
                                                            mean(abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum))),
                                                            mean(abs(apply(elasticity.stochastic.surv.j.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "J survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum)), probs = 0.025),
                                                           quantile(abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum)), probs = 0.025),
                                                           quantile(abs(apply(elasticity.stochastic.surv.j.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "J survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum)), probs = 0.975),
                                                          quantile(abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum)), probs = 0.975),
                                                          quantile(abs(apply(elasticity.stochastic.surv.j.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.3. Juvenile wet-season survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "J survival (wet)")] = c(mean(abs(apply(elasticity.mean.surv.j.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.j.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.j.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "J survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.j.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.j.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.j.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "J survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.j.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.j.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.j.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.4. Subadult dry-season survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "SA survival (dry)")] = c(mean(abs(apply(elasticity.mean.surv.sa.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.sa.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.sa.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "SA survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.sa.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.sa.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.sa.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "SA survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.sa.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.sa.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.sa.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.5 Subadult wet-season survival ----
# -----------------------------------

elast.df$elast.mean[which(elast.df$vr == "SA survival (wet)")] = c(mean(abs(apply(elasticity.mean.surv.sa.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.sa.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.sa.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "SA survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.sa.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.sa.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.sa.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "SA survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.sa.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.sa.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.sa.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.6. Helper dry-season survival ----
# ----------------------------------

elast.df$elast.mean[which(elast.df$vr == "H survival (dry)")] = c(mean(abs(apply(elasticity.mean.surv.h.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.h.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.h.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.h.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.h.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.h.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.h.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.h.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.h.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.7. Helper wet-season survival ----
# ----------------------------------

elast.df$elast.mean[which(elast.df$vr == "H survival (wet)")] = c(mean(abs(apply(elasticity.mean.surv.h.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.h.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.h.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.h.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.h.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.h.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.h.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.h.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.h.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.8. Dominant dry-season survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "D survival (dry)")] = c(mean(abs(apply(elasticity.mean.surv.d.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.d.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.d.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "D survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.d.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.d.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.d.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "D survival (dry)")] = c(quantile(abs(apply(elasticity.mean.surv.d.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.d.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.d.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.9. Dominant wet-season survival ----
# ------------------------------------

elast.df$elast.mean[which(elast.df$vr == "D survival (wet)")] = c(mean(abs(apply(elasticity.mean.surv.d.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.surv.d.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.surv.d.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "D survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.d.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.surv.d.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.surv.d.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "D survival (wet)")] = c(quantile(abs(apply(elasticity.mean.surv.d.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.surv.d.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.surv.d.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.10. H-D dry-season transition ----
# ----------------------------------

elast.df$elast.mean[which(elast.df$vr == "H-D transition (dry)")] = c(mean(abs(apply(elasticity.mean.trans.hd.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.trans.hd.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.trans.hd.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H-D transition (dry)")] = c(quantile(abs(apply(elasticity.mean.trans.hd.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.trans.hd.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.trans.hd.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H-D transition (dry)")] = c(quantile(abs(apply(elasticity.mean.trans.hd.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.trans.hd.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.trans.hd.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.11. H-D wet-season transition ----
# ----------------------------------

elast.df$elast.mean[which(elast.df$vr == "H-D transition (wet)")] = c(mean(abs(apply(elasticity.mean.trans.hd.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.trans.hd.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.trans.hd.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H-D transition (wet)")] = c(quantile(abs(apply(elasticity.mean.trans.hd.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.trans.hd.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.trans.hd.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H-D transition (wet)")] = c(quantile(abs(apply(elasticity.mean.trans.hd.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.trans.hd.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.trans.hd.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.12. Helper dry-season emigration ----
# -------------------------------------

elast.df$elast.mean[which(elast.df$vr == "H emigration (dry)")] = c(mean(abs(apply(elasticity.mean.emig.h.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.emig.h.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.emig.h.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H emigration (dry)")] = c(quantile(abs(apply(elasticity.mean.emig.h.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.emig.h.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.emig.h.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H emigration (dry)")] = c(quantile(abs(apply(elasticity.mean.emig.h.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.emig.h.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.emig.h.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.13. Helper wet-season emigration ----
# -------------------------------------

elast.df$elast.mean[which(elast.df$vr == "H emigration (wet)")] = c(mean(abs(apply(elasticity.mean.emig.h.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.emig.h.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.emig.h.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H emigration (wet)")] = c(quantile(abs(apply(elasticity.mean.emig.h.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.emig.h.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.emig.h.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H emigration (wet)")] = c(quantile(abs(apply(elasticity.mean.emig.h.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.emig.h.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.emig.h.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.14. Helper dry-season recruitment ----
# --------------------------------------

elast.df$elast.mean[which(elast.df$vr == "H recruitment (dry)")] = c(mean(abs(apply(elasticity.mean.recruit.h.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.recruit.h.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.recruit.h.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H recruitment (dry)")] = c(quantile(abs(apply(elasticity.mean.recruit.h.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.recruit.h.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.recruit.h.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H recruitment (dry)")] = c(quantile(abs(apply(elasticity.mean.recruit.h.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.recruit.h.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.recruit.h.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.15. Helper wet-season recruitment ----
# --------------------------------------

elast.df$elast.mean[which(elast.df$vr == "H recruitment (wet)")] = c(mean(abs(apply(elasticity.mean.recruit.h.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.recruit.h.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.recruit.h.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "H recruitment (wet)")] = c(quantile(abs(apply(elasticity.mean.recruit.h.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.recruit.h.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.recruit.h.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "H recruitment (wet)")] = c(quantile(abs(apply(elasticity.mean.recruit.h.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.recruit.h.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.recruit.h.wet, 3, FUN = sum)), probs = 0.975))


## 4.2.16. Dominant dry-season recruitment ----
# ----------------------------------------

elast.df$elast.mean[which(elast.df$vr == "D recruitment (dry)")] = c(mean(abs(apply(elasticity.mean.recruit.d.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.recruit.d.dry, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.recruit.d.dry, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "D recruitment (dry)")] = c(quantile(abs(apply(elasticity.mean.recruit.d.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.recruit.d.dry, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.recruit.d.dry, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "D recruitment (dry)")] = c(quantile(abs(apply(elasticity.mean.recruit.d.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.recruit.d.dry, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.recruit.d.dry, 3, FUN = sum)), probs = 0.975))


## 4.2.17. Dominant wet-season recruitment ----
# ----------------------------------------

elast.df$elast.mean[which(elast.df$vr == "D recruitment (wet)")] = c(mean(abs(apply(elasticity.mean.recruit.d.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.sd.recruit.d.wet, 3, FUN = sum))),
                                                                  mean(abs(apply(elasticity.stochastic.recruit.d.wet, 3, FUN = sum))))

elast.df$elast.low[which(elast.df$vr == "D recruitment (wet)")] = c(quantile(abs(apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)), probs = 0.025),
                                                                 quantile(abs(apply(elasticity.stochastic.recruit.d.wet, 3, FUN = sum)), probs = 0.025))

elast.df$elast.up[which(elast.df$vr == "D recruitment (wet)")] = c(quantile(abs(apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)), probs = 0.975),
                                                                quantile(abs(apply(elasticity.stochastic.recruit.d.wet, 3, FUN = sum)), probs = 0.975))


## 4.3. Processing relative effect of variability results ----
# -------------------------------------------------------

# Following Morris et al. (2008), we compute the ratio: sum(elasticity to SD of all vital rates in one of the survival/transition/reproductive rates categories) / sum(elasticity to the mean + elasticity to SD of all vital rates)

## 4.3.1. Survival rates ----
# ----------------------

relative.elast.surv = abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                          apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                          apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                          apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                          apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                          apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                          apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                          apply(elasticity.sd.surv.d.wet, 3, FUN = sum)) / (abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                                                                                apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                                                                                apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                                                                                apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.surv.d.wet, 3, FUN = sum) +
                                                                                apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.trans.hd.wet, 3, FUN = sum) +
                                                                                apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                                                                                apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                                                                                apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) + abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.surv.j.wet, 3, FUN = sum) + 
                                                                                                                                        apply(elasticity.mean.surv.sa.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.surv.sa.wet, 3, FUN = sum) + 
                                                                                                                                        apply(elasticity.mean.surv.h.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.surv.h.wet, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.surv.d.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.surv.d.wet, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.trans.hd.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.trans.hd.wet, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.recruit.h.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.recruit.h.wet, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.recruit.d.dry, 3, FUN = sum) +
                                                                                                                                        apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)))


## 4.3.2. Reproductive rates ----
# --------------------------

relative.elast.repro = abs(apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                           apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                           apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                           apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) / (abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                                                                                    apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                                                                                    apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.surv.d.wet, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.trans.hd.wet, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                                                                                    apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) + abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.surv.j.wet, 3, FUN = sum) + 
                                                                                                                                            apply(elasticity.mean.surv.sa.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.surv.sa.wet, 3, FUN = sum) + 
                                                                                                                                            apply(elasticity.mean.surv.h.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.surv.h.wet, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.surv.d.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.surv.d.wet, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.trans.hd.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.trans.hd.wet, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.recruit.h.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.recruit.h.wet, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.recruit.d.dry, 3, FUN = sum) +
                                                                                                                                            apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)))


## 4.3.3. Transition rate ----
# -----------------------

relative.elast.trans = abs(apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                           apply(elasticity.sd.trans.hd.wet, 3, FUN = sum)) / (abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                                                                                   apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                                                                                   apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.surv.d.wet, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.trans.hd.wet, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                                                                                   apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) + abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.surv.j.wet, 3, FUN = sum) + 
                                                                                                                                           apply(elasticity.mean.surv.sa.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.surv.sa.wet, 3, FUN = sum) + 
                                                                                                                                           apply(elasticity.mean.surv.h.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.surv.h.wet, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.surv.d.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.surv.d.wet, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.trans.hd.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.trans.hd.wet, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.recruit.h.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.recruit.h.wet, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.recruit.d.dry, 3, FUN = sum) +
                                                                                                                                           apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)))


## 4.3.4. Emigration ----
# ------------------

relative.elast.emig = abs(apply(elasticity.sd.emig.h.dry, 3, FUN = sum) +
                           apply(elasticity.sd.emig.h.wet, 3, FUN = sum)) / (abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                                                                                       apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                                                                                       apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.surv.d.wet, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.trans.hd.wet, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                                                                                       apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) + abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.surv.j.wet, 3, FUN = sum) + 
                                                                                                                                                 apply(elasticity.mean.surv.sa.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.surv.sa.wet, 3, FUN = sum) + 
                                                                                                                                                 apply(elasticity.mean.surv.h.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.surv.h.wet, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.surv.d.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.surv.d.wet, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.trans.hd.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.trans.hd.wet, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.recruit.h.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.recruit.h.wet, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.recruit.d.dry, 3, FUN = sum) +
                                                                                                                                                 apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)))


## 4.3.5. All vital rates ----
# -----------------------

relative.elast.all.vr = abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                            apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                            apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                            apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                            apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                            apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                            apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                            apply(elasticity.sd.surv.d.wet, 3, FUN = sum) +
                            apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                            apply(elasticity.sd.trans.hd.wet, 3, FUN = sum) +
                            apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                            apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                            apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                            apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) / (abs(apply(elasticity.sd.surv.j.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.surv.j.wet, 3, FUN = sum) + 
                                                                                     apply(elasticity.sd.surv.sa.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.surv.sa.wet, 3, FUN = sum) + 
                                                                                     apply(elasticity.sd.surv.h.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.surv.h.wet, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.surv.d.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.surv.d.wet, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.trans.hd.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.trans.hd.wet, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.recruit.h.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.recruit.h.wet, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.recruit.d.dry, 3, FUN = sum) +
                                                                                     apply(elasticity.sd.recruit.d.wet, 3, FUN = sum)) + abs(apply(elasticity.mean.surv.j.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.surv.j.wet, 3, FUN = sum) + 
                                                                                                                                             apply(elasticity.mean.surv.sa.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.surv.sa.wet, 3, FUN = sum) + 
                                                                                                                                             apply(elasticity.mean.surv.h.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.surv.h.wet, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.surv.d.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.surv.d.wet, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.trans.hd.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.trans.hd.wet, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.recruit.h.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.recruit.h.wet, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.recruit.d.dry, 3, FUN = sum) +
                                                                                                                                             apply(elasticity.mean.recruit.d.wet, 3, FUN = sum)))

relative.elast.df = data.frame(vr = c("Survival", "Reproduction", "Transition", "Emigration", "All vital rates"),
                               relative.elast.mean = c(mean(relative.elast.surv), mean(relative.elast.repro), mean(relative.elast.trans), mean(relative.elast.emig), mean(relative.elast.all.vr)), 
                               relative.elast.low = c(quantile(relative.elast.surv, probs = 0.025), quantile(relative.elast.repro, probs = 0.025), quantile(relative.elast.trans, probs = 0.025), quantile(relative.elast.emig, probs = 0.025), quantile(relative.elast.all.vr, probs = 0.025)), 
                               relative.elast.up = c(quantile(relative.elast.surv, probs = 0.975), quantile(relative.elast.repro, probs = 0.975), quantile(relative.elast.trans, probs = 0.975), quantile(relative.elast.emig, probs = 0.975), quantile(relative.elast.all.vr, probs = 0.975)))




###########################################################################
#
# 5. Plotting the results -----
#
###########################################################################

## 5.1. Preparing plotting ----
# ------------------------

elast.df$elast.type = factor(elast.df$elast.type, levels = c("Stochastic", "Mean", "SD"))
elast.df$elast.type.parsed = c(rep("S\U003BC", 16), rep("S\u03C3", 16), rep("S", 16))


## 5.2. Elasticities ----
# ------------------

png(filename = "Meerkats_ElastMeanSDStoch.png",
    width=7000,
    height=3000,
    units="px",
    bg="white",
    res=300,
    type="cairo")

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

png(filename = "Meerkats_RelativeEffectVariability.png",
    width=4000,
    height=3000,
    units="px",
    bg="white",
    res=300,
    type="cairo")

ggplot(relative.elast.df, aes(x = vr, y = relative.elast.mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = relative.elast.low, ymax = relative.elast.up, width = 0.25)) +
  labs(x = "Vital rate category", y="Relative effect of variability") +
  theme_bw() +
  theme(axis.title.x = element_text(size=25,colour="black", margin = margin(t=18,r=0,b=0,l=0),face="bold"), axis.title.y = element_text(size=25,colour="black", margin = margin(t=0,r=10,b=0,l=0),face="bold"), axis.text.x = element_text(size=22,colour="black", margin = margin(t=10,r=0,b=0,l=0), hjust=0.5), axis.text.y = element_text(size=22,colour="black", margin = margin(t=0,r=10,b=0,l=0)), legend.text = element_text(size=22), legend.title = element_text(face="bold",size=25), legend.position = "right", legend.key.size = unit(6,"lines"),strip.text.x = element_text(size=20,face="bold"))

dev.off()




###########################################################################
#
# 6. Saving the results -----
#
###########################################################################

write.csv(elast.df, "Meerkats_Elasticities.csv", row.names = F)
write.csv(relative.elast.df, "Meerkats_RelativeElasticities.csv", row.names = F)
