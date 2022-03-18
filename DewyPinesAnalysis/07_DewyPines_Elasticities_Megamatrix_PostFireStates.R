##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
#
# This script uses the post-fire state-specific vital rate estimates for the dewy pine population. 
# The aim of this script is to calculate the sensitivity of the population to each post-fire state using 
# a megamatrix, following Pascarella and Horvitz (1998).
#
# Author: Eva Conquet
#
# Pascarella, J. B., and C. C. Horvitz. 1998. Hurricane disturbance and the population
# dynamics of a tropical understory shrub: megamatrix elasticity analysis. 
# Ecology 79: 547–563.
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
  library(Matrix)
  library(popbio)
}

load.librairies()


## 1.3. Loading and preparing data ----
# --------------------------------

data.seedbank = read.csv("DPVR_Seedbank_SitesEF.csv", stringsAsFactors = F, na.strings = c("", "NA"))
data.trans.seedbank.TSF0 = read.csv("DPVR_SeedbankTransitions_TSF0.csv", stringsAsFactors = F, na.strings = c("", "NA"))
data.vr.TSF1.to.TSF3 = read.csv("DPVR_TSF1_to_3.csv", stringsAsFactors = F, na.strings = c("", "NA"))
data.vr.TSF3plus = read.csv("DPVR_TSF3plus.csv", stringsAsFactors = F, na.strings = c("", "NA"))


seeds.per.flower = 9.8

head(data.seedbank)
head(data.trans.seedbank.TSF0)
head(data.vr.TSF1.to.TSF3)
head(data.vr.TSF3plus)




###########################################################################
#
# 2. Preparing the simulations ----
#
###########################################################################

## 2.1. Creating functions to create the MPMs ----
# -------------------------------------------

# TSF0 ----

buildmat.TSF0 <- function(site){
  
  TSF0.mat = matrix(0, nrow = 5, ncol = 5)
  colnames(TSF0.mat) = c("SB", "SD", "J", "SR", "LR")
  rownames(TSF0.mat) = colnames(TSF0.mat)
  
  TSF0.mat["SB", "SB"] = data.seedbank$staySB[data.seedbank$TSF == "zero"]
  
  TSF0.mat["SD", "SB"] = data.seedbank$outSB[data.seedbank$TSF == "zero"] * data.trans.seedbank.TSF0$transSB_SD[which(data.trans.seedbank.TSF0$site == site)]
  
  TSF0.mat["J", "SB"] = data.seedbank$outSB[data.seedbank$TSF == "zero"] * data.trans.seedbank.TSF0$transSB_J[which(data.trans.seedbank.TSF0$site == site)]
  
  return(TSF0.mat)
  
}


# TSF1 ----

buildmat.TSF1 <- function(site){
  
  TSF1.mat = matrix(0, nrow = 5, ncol = 5)
  colnames(TSF1.mat) = c("SB", "SD", "J", "SR", "LR")
  rownames(TSF1.mat) = colnames(TSF1.mat)
  
  # In/out seedbank (SB)
  TSF1.mat["SB", "SB"] = data.seedbank$staySB[data.seedbank$TSF == "one"]
  
  TSF1.mat["SD", "SB"] = data.seedbank$outSB[data.seedbank$TSF == "one"]
  
  # Seedlings (SD) to Juveniles (J)
  TSF1.mat["J", "SD"] = data.vr.TSF1.to.TSF3$survSD[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF1")]
  
  # Juveniles (J) to Reproductive (SR/LR)
  TSF1.mat["SR", "J"] = data.vr.TSF1.to.TSF3$survJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF1")] * (1 - data.vr.TSF1.to.TSF3$transJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF1")])
    
  TSF1.mat["LR", "J"] = data.vr.TSF1.to.TSF3$survJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF1")] * data.vr.TSF1.to.TSF3$transJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF1")]
  
  return(TSF1.mat)
  
}


# TSF2 ----

buildmat.TSF2 <- function(site){
  
  TSF2.mat = matrix(0, nrow = 5, ncol = 5)
  colnames(TSF2.mat) = c("SB", "SD", "J", "SR", "LR")
  rownames(TSF2.mat) = colnames(TSF2.mat)
  
  # In/out seedbank (SB)
  TSF2.mat["SB", "SB"] = data.seedbank$staySB[data.seedbank$TSF == "two"]
  
  TSF2.mat["SD", "SB"] = data.seedbank$outSB[data.seedbank$TSF == "two"]
  
  # Seedlings (SD) to Juveniles (J)
  TSF2.mat["J", "SD"] = data.vr.TSF1.to.TSF3$survSD[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")]

  # Juveniles (J) to Reproductive (SR/LR)
  TSF2.mat["SR", "J"] = data.vr.TSF1.to.TSF3$survJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * (1 - data.vr.TSF1.to.TSF3$transJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")])
    
  TSF2.mat["LR", "J"] = data.vr.TSF1.to.TSF3$survJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$transJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")]
  
  # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
  TSF2.mat["SR", "SR"] = data.vr.TSF1.to.TSF3$survSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * (1 - data.vr.TSF1.to.TSF3$transSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")])
   
  TSF2.mat["LR", "SR"] = data.vr.TSF1.to.TSF3$survSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$transSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")]
  
  # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
  TSF2.mat["LR", "LR"] = data.vr.TSF1.to.TSF3$survLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * (1 - data.vr.TSF1.to.TSF3$transLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")])
    
  TSF2.mat["SR", "LR"] = data.vr.TSF1.to.TSF3$survLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$transLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")]
  
  # Reproduction of Small Reproductive (SR)
  TSF2.mat["SB", "SR"] = data.vr.TSF1.to.TSF3$floweringprobSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfpsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * seeds.per.flower * data.seedbank$goSB[which(data.seedbank$TSF == "two")]
    
  TSF2.mat["SD", "SR"] = data.vr.TSF1.to.TSF3$floweringprobSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfpsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * seeds.per.flower * data.seedbank$goCont[which(data.seedbank$TSF == "two")]
  
  # Reproduction of Large Reproductive (LR)
  TSF2.mat["SB", "LR"] = data.vr.TSF1.to.TSF3$floweringprobLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfpsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * seeds.per.flower * data.seedbank$goSB[which(data.seedbank$TSF == "two")]
  
  TSF2.mat["SD", "LR"] = data.vr.TSF1.to.TSF3$floweringprobLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * data.vr.TSF1.to.TSF3$nbfpsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF2")] * seeds.per.flower * data.seedbank$goCont[which(data.seedbank$TSF == "two")]
  
  return(TSF2.mat)
  
}


# TSF3 ----

buildmat.TSF3 <- function(site){
  
  TSF3.mat = matrix(0, nrow = 5, ncol = 5)
  colnames(TSF3.mat) = c("SB", "SD", "J", "SR", "LR")
  rownames(TSF3.mat) = colnames(TSF3.mat)
  
  # In/out seedbank (SB)
  TSF3.mat["SB", "SB"] = data.seedbank$staySB[data.seedbank$TSF == "three"]
  
  TSF3.mat["SD", "SB"] = data.seedbank$outSB[data.seedbank$TSF == "three"]
  
  # Seedlings (SD) to Juveniles (J)
  TSF3.mat["J", "SD"] = data.vr.TSF1.to.TSF3$survSD[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")]
  
  # Juveniles (J) to Reproductive (SR/LR)
  TSF3.mat["SR", "J"] = data.vr.TSF1.to.TSF3$survJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * (1 - data.vr.TSF1.to.TSF3$transJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")])
  
  TSF3.mat["LR", "J"] = data.vr.TSF1.to.TSF3$survJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$transJ[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")]
  
  # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
  TSF3.mat["SR", "SR"] = data.vr.TSF1.to.TSF3$survSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * (1 - data.vr.TSF1.to.TSF3$transSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")])
  
  TSF3.mat["LR", "SR"] = data.vr.TSF1.to.TSF3$survSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$transSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")]
  
  # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
  TSF3.mat["LR", "LR"] = data.vr.TSF1.to.TSF3$survLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * (1 - data.vr.TSF1.to.TSF3$transLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")])
  
  TSF3.mat["SR", "LR"] = data.vr.TSF1.to.TSF3$survLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$transLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")]
  
  # Reproduction of Small Reproductive (SR)
  TSF3.mat["SB", "SR"] = data.vr.TSF1.to.TSF3$floweringprobSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfpsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * seeds.per.flower * data.seedbank$goSB[which(data.seedbank$TSF == "three")]
  
  TSF3.mat["SD", "SR"] = data.vr.TSF1.to.TSF3$floweringprobSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfpsSR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * seeds.per.flower * data.seedbank$goCont[which(data.seedbank$TSF == "three")]
  
  # Reproduction of Large Reproductive (LR)
  TSF3.mat["SB", "LR"] = data.vr.TSF1.to.TSF3$floweringprobLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfpsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * seeds.per.flower * data.seedbank$goSB[which(data.seedbank$TSF == "three")]
  
  TSF3.mat["SD", "LR"] = data.vr.TSF1.to.TSF3$floweringprobLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * data.vr.TSF1.to.TSF3$nbfpsLR[which(data.vr.TSF1.to.TSF3$site == site & data.vr.TSF1.to.TSF3$TSF == "TSF3")] * seeds.per.flower * data.seedbank$goCont[which(data.seedbank$TSF == "three")]
  
  return(TSF3.mat)
  
}


# TSF>3 ----

buildmat.stochastic <- function(year){
  
  stochastic.mat = matrix(0, nrow = 5, ncol = 5)
  colnames(stochastic.mat) = c("SB", "SD", "J", "SR", "LR")
  rownames(stochastic.mat) = colnames(stochastic.mat)
  
  # In/out seedbank (SB)
  stochastic.mat["SB", "SB"] = data.seedbank$staySB[data.seedbank$TSF == "three"]
  
  stochastic.mat["SD", "SB"] = data.seedbank$outSB[data.seedbank$TSF == "three"]
  
  # Seedlings (SD) to Juveniles (J)
  stochastic.mat["J", "SD"] = data.vr.TSF3plus$survSD[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")]
  
  # Juveniles (J) to Reproductive (SR/LR)
  stochastic.mat["SR", "J"] = data.vr.TSF3plus$survJ[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * (1 - data.vr.TSF3plus$transJ[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")])
  
  stochastic.mat["LR", "J"] = data.vr.TSF3plus$survJ[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$transJ[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")]
  
  # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
  stochastic.mat["SR", "SR"] = data.vr.TSF3plus$survSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * (1 - data.vr.TSF3plus$transSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")])
  
  stochastic.mat["LR", "SR"] = data.vr.TSF3plus$survSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$transSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")]
  
  # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
  stochastic.mat["LR", "LR"] = data.vr.TSF3plus$survLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * (1 - data.vr.TSF3plus$transLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")])
  
  stochastic.mat["SR", "LR"] = data.vr.TSF3plus$survLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$transLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")]
  
  # Reproduction of Small Reproductive (SR)
  stochastic.mat["SB", "SR"] = data.vr.TSF3plus$floweringprobSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfsSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfpsSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * seeds.per.flower * data.seedbank$goSB[which(data.seedbank$TSF == "three")]
  
  stochastic.mat["SD", "SR"] = data.vr.TSF3plus$floweringprobSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfsSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfpsSR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * seeds.per.flower * data.seedbank$goCont[which(data.seedbank$TSF == "three")]
  
  # Reproduction of Large Reproductive (LR)
  stochastic.mat["SB", "LR"] = data.vr.TSF3plus$floweringprobLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfsLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfpsLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * seeds.per.flower * data.seedbank$goSB[which(data.seedbank$TSF == "three")]
  
  stochastic.mat["SD", "LR"] = data.vr.TSF3plus$floweringprobLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfsLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * data.vr.TSF3plus$nbfpsLR[which(data.vr.TSF3plus$year == year & data.vr.TSF3plus$grazing == "low")] * seeds.per.flower * data.seedbank$goCont[which(data.seedbank$TSF == "three")]
  
  return(stochastic.mat)
  
}


## 2.2. Creating the megamatrix and calculating elasticities (following Pascarella and Horvitz, 1998) ----
# ---------------------------------------------------------------------------------------------------

# 5 life-history stages
# 5 post-fire habitat states

# (1) Define a 5x5 identity matrix (I5)
I5 = diag(5)

# (2) Define a 5x5 matrix with only zeros (Z5)
Z5 = matrix(0, nrow = 5, ncol = 5)

# (3) Define matrices for each post-fire habitat state (A1 to A5)
A1_TSF0 = buildmat.TSF0(site = sample(c("siteE", "siteF"), 1))
A2_TSF1 = buildmat.TSF1(site = sample(c("siteE", "siteF"), 1))
A3_TSF2 = buildmat.TSF2(site = sample(c("siteE", "siteF"), 1))
A4_TSF3 = buildmat.TSF3(site = sample(c("siteE", "siteF"), 1))

M.elasticities.TSF0 = rep(0, 1000)
M.elasticities.TSF1 = rep(0, 1000)
M.elasticities.TSF2 = rep(0, 1000)
M.elasticities.TSF3 = rep(0, 1000)
M.elasticities.TSFstochastic = rep(0, 1000)

for(i in 1:1000){
  
  year = sample(data.vr.TSF3plus$year, 1)
  A5_stochastic = buildmat.stochastic(year = year)
  
  # (4) Define a 25x25 matrix composed of the Ai's and Z8's (A15)
  A15 = bdiag(list(A1_TSF0, A2_TSF1, A3_TSF2, A4_TSF3, A5_stochastic))
  
  # (5) Define the habitat state-transition matrix (C)
  
  # Probability of fire when in TSF>2 ~ once every 10 years 
  p = 1/30
  
  # Matrix for the environment transitions
  C = matrix(0, 5, 5)
  
  C[2, 1] = C[3, 2] = C[4, 3] = C[5, 4] = 1
  C[1, 5] = p
  C[5, 5] = 1 - p
  
  colnames(C) = rownames(C) = 1:5
  
  # (6) Define the megamatrix (M)
  M = (C %x% I5) %*% A15
  
  # Getting elasticities from the megamatrix
  
  M.elasticities = elasticity(M)
  
  M.elasticities.TSF0[i] = sum(M.elasticities[, 1:5])
  M.elasticities.TSF1[i] = sum(M.elasticities[, 6:10])
  M.elasticities.TSF2[i] = sum(M.elasticities[, 11:15])
  M.elasticities.TSF3[i] = sum(M.elasticities[, 16:20])
  M.elasticities.TSFstochastic[i] = sum(M.elasticities[, 21:25])

}

M.elasticities.df = data.frame(TSF = c("TSF0", "TSF1", "TSF2", "TSF3", "TSF>3"), 
                               elast.mean = c(mean(M.elasticities.TSF0), mean(M.elasticities.TSF1), mean(M.elasticities.TSF2), mean(M.elasticities.TSF3), mean(M.elasticities.TSFstochastic)),
                               elast.low = c(quantile(M.elasticities.TSF0, probs = 0.025), quantile(M.elasticities.TSF1, probs = 0.025), quantile(M.elasticities.TSF2, probs = 0.025), quantile(M.elasticities.TSF3, probs = 0.025), quantile(M.elasticities.TSFstochastic, probs = 0.025)),
                               elast.up = c(quantile(M.elasticities.TSF0, probs = 0.975), quantile(M.elasticities.TSF1, probs = 0.975), quantile(M.elasticities.TSF2, probs = 0.975), quantile(M.elasticities.TSF3, probs = 0.975), quantile(M.elasticities.TSFstochastic, probs = 0.975)))

M.elasticities.df$TSF = factor(M.elasticities.df$TSF, levels = unique(M.elasticities.df$TSF))




###########################################################################
#
# 3. Plotting the results ----
#
###########################################################################

tiff(filename = "DP_ElasticityToHabitatState.tiff",
    width=3,
    height=3,
    units="in",
    bg="white",
    res=600,
    compression="lzw")

ggplot(M.elasticities.df, aes(x = TSF, y = elast.mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = elast.low, ymax = elast.up, width = 0.25)) +
  labs(x = "Time since fire (TSF)", y = "Elasticity to the post-fire habitat") +
  scale_x_discrete(labels = c(expression("TSF"["0"]), expression("TSF"["1"]), expression("TSF"["2"]), expression("TSF"["3"]), expression("TSF"[">3"]))) +
  ylim(0, 0.6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 9, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        axis.text.x = element_text(size = 7, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0), angle = 0, hjust = 0.5), 
        axis.text.y = element_text(size = 7, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)), 
        legend.text = element_text(size = 7), 
        legend.title = element_text( size = 9), 
        legend.position = "right", 
        legend.key.size = unit(2, "lines"),
        strip.text.x = element_text(size = 10))

dev.off()




###########################################################################
#
# 4. Saving the files ----
#
###########################################################################

write.csv(M.elasticities.df, file = "DewyPine_ElasticitiesToHabitatState.csv")
