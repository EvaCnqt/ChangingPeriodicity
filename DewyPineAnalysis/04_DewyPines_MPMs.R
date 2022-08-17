##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al. 2022, Ecology).
# 
# Eva Conquet, Arpat Ozgul, Daniel T. Blumstein, Kenneth B. Armitage,
# Madan K. Oli, Julien G. A. Martin, Tim H. Clutton-Brock, and Maria Paniw
#
# This script contains the functions needed to build the matrix population models 
# for the dewy pine population.
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


## 1.2. Loading data ----
# ------------------


## 1.2.1. Seed bank and TSF0 constants and number of seeds per flower ----
# -------------------------------------------------------------------
seed.bankSiteC = read.csv("DPVR_Seedbank_SiteC.csv")
seed.bankSitesEF = read.csv("DPVR_Seedbank_SitesEF.csv")

sb.transition.TSF0 = read.csv("DPVR_SeedbankTransitions_TSF0.csv")

load("NbSeedsPerFlower.RData")


## 1.2.2. TSF1 to 3 vital rates tables and models ----
# -----------------------------------------------

vr.per.TSF.per.site = read.csv("DPVR_TSF1_to_3.csv")

load("GLM_survSD_siteC_TSF1.RData")

load("GLM_survSD_siteE_TSF1_3.RData")

load("GLM_survSD_siteF_TSF3.RData")

load("GLM_survJ_siteC_TSF1.RData")

load("GLM_survJ_siteE.RData")

load("GLM_survJ_siteF_TSF1_3.RData")

load("GLM_survSR_siteC.RData")

load("GLM_survSR_siteE_TSF3.RData")

load("GLM_survLR_siteC_TSF3.RData")

load("GLM_survLR_siteE_TSF3.RData")

load("GLM_survLR_siteF.RData")

load("GLM_transJ_siteC_TSF1.RData")

load("GLM_transJ_siteE.RData")

load("GLM_transJ_siteF_TSF1.RData")

load("GLM_transSR_siteE_TSF3.RData")

load("GLM_transLR_siteC_TSF3.RData")

load("GLM_transLR_siteE_TSF3.RData")

load("GLM_transLR_siteF_TSF2.RData")

load("GLM_floweringprobSR_siteC.RData")

load("GLM_floweringprobSR_siteE_TSF3.RData")

load("GLM_floweringprobLR_siteC_TSF3.RData")

load("GLM_floweringprobLR_siteE_TSF3.RData")

load("GLM_floweringprobLR_siteF.RData")

load("GLM_nbfsSR_siteE_TSF3.RData")

load("GLM_nbfsLR_siteC_TSF3.RData")

load("GLM_nbfsLR_siteE_TSF3.RData")

load("GLM_nbfsLR_siteF_TSF3.RData")

load("GLM_nbfpsSR_siteE_TSF3.RData")

load("GLM_nbfpsLR_siteE_TSF3.RData")

load("GLM_nbfpsLR_siteF_TSF3.RData")


## 1.2.3. TSF>3 vital rates tables and models ----
# -------------------------------------------

vr.per.year.TSF3plus = read.csv("DPVR_TSF3plus.csv")

load("GLMM_survSD_LG.RData")

load("GLMM_survSD_HG.RData")

load("GLMM_survJ_LG.RData")

load("GLMM_survJ_HG.RData")

load("GLMM_survSR_LG.RData")

load("GLMM_survSR_HG.RData")

load("GLMM_survLR_LG.RData")

load("GLMM_survLR_HG.RData")

load("GLMM_transJ_LG.RData")

load("GLMM_transJ_HG.RData")

load("GLMM_transSR_LG.RData")

load("GLMM_transSR_LG.RData")

load("GLMM_transLR_LG.RData")

load("GLMM_transLR_HG.RData")

load("GLMM_floweringprobSR_LG.RData")

load("GLMM_floweringprobSR_HG.RData")

load("GLMM_floweringprobLR_LG.RData")

load("GLMM_floweringprobLR_HG.RData")

load("GLMM_nbfsLR_LG.RData")

load("GLMM_nbfpsSR_HG.RData")

load("GLMM_nbfpsLR_LG.RData")

load("GLMM_nbfpsLR_HG.RData")




###########################################################################
#
# 2. Functions to build matrices for TSF 0 to 3 in each site ----
#
###########################################################################

## 2.1. Site C - TSF0 ----
# -------------------

dewy.buildmat.TSF0.siteC <- function(){
 
 # Creating the matrix
 TSF0.mat.siteC = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF0.mat.siteC) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF0.mat.siteC) = colnames(TSF0.mat.siteC)
 
 # Staying in the seedbank
 TSF0.mat.siteC["SB", "SB"] = seed.bankSiteC$staySB[seed.bankSiteC$TSF == "zero"]
 
 # Germinating to seedling
 TSF0.mat.siteC["SD", "SB"] = seed.bankSiteC$outSB[seed.bankSiteC$TSF == "zero"] * sb.transition.TSF0$transSB_SD[sb.transition.TSF0$site == "siteC"]
 
 # Germinating and growing to juvenile
 TSF0.mat.siteC["J", "SB"] = seed.bankSiteC$outSB[seed.bankSiteC$TSF == "zero"] * sb.transition.TSF0$transSB_J[sb.transition.TSF0$site == "siteC"]
 
 return(TSF0.mat.siteC)
}


## 2.2. Site C - TSF1 ----
# -------------------

dewy.buildmat.TSF1.siteC <- function(density){
  
 # Creating the matrix
 TSF1.mat.siteC = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF1.mat.siteC) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF1.mat.siteC) = colnames(TSF1.mat.siteC)
 
 # In/out seedbank (SB)
 TSF1.mat.siteC["SB", "SB"] = seed.bankSiteC$staySB[seed.bankSiteC$TSF == "one"]
 TSF1.mat.siteC["SD", "SB"] = seed.bankSiteC$outSB[seed.bankSiteC$TSF == "one"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF1.mat.siteC["J", "SD"] = predict(survSD.siteC.TSF1, newdata = data.frame(density = density), type = "response")
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF1.mat.siteC["SR", "J"] = unique(predict(survJ.siteC.TSF1, type = "response")) * (1 - unique(predict(transJ.siteC.TSF1, type = "response")))
 TSF1.mat.siteC["LR", "J"] = unique(predict(survJ.siteC.TSF1, type = "response")) * unique(predict(transJ.siteC.TSF1, type = "response"))
 
 return(TSF1.mat.siteC)
}


## 2.3. Site C - TSF2 ----
# -------------------

dewy.buildmat.TSF2.siteC <- function(density){
  
 # Creating the matrix
 TSF2.mat.siteC = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF2.mat.siteC) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF2.mat.siteC) = colnames(TSF2.mat.siteC)
 
 # In/out seedbank (SB)
 TSF2.mat.siteC["SB", "SB"] = seed.bankSiteC$staySB[seed.bankSiteC$TSF == "two"]
 TSF2.mat.siteC["SD", "SB"] = seed.bankSiteC$outSB[seed.bankSiteC$TSF == "two"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF2.mat.siteC["J", "SD"] = vr.per.TSF.per.site$survSD[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF2.mat.siteC["SR", "J"] = vr.per.TSF.per.site$survJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * (1 - vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteC["LR", "J"] = vr.per.TSF.per.site$survJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
 TSF2.mat.siteC["SR", "SR"] = predict(survSR.siteC, newdata = data.frame(TSF = "two", density = density), type = "response") * (1 - vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteC["LR", "SR"] = predict(survSR.siteC, newdata = data.frame(TSF = "two", density = density), type = "response") * vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
 TSF2.mat.siteC["LR", "LR"] = vr.per.TSF.per.site$survLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * (1 - vr.per.TSF.per.site$transLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteC["SR", "LR"] = vr.per.TSF.per.site$survLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$transLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Reproduction of Small Reproductive (SR)
 TSF2.mat.siteC["SB", "SR"] = predict(floweringprobSR.siteC, newdata = data.frame(TSF = "two"), type = "response") * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSiteC$goSB[seed.bankSiteC$TSF == "two"]
 TSF2.mat.siteC["SD", "SR"] = predict(floweringprobSR.siteC, newdata = data.frame(TSF = "two"), type = "response") * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSiteC$goCont[seed.bankSiteC$TSF == "two"]
 
 # Reproduction of Large Reproductive (LR)
 TSF2.mat.siteC["SB", "LR"] = vr.per.TSF.per.site$floweringprobLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSiteC$goSB[seed.bankSiteC$TSF == "two"]
 TSF2.mat.siteC["SD", "LR"] = vr.per.TSF.per.site$floweringprobLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSiteC$goCont[seed.bankSiteC$TSF == "two"]
 
 return(TSF2.mat.siteC)
}


## 2.4. Site C - TSF3 ----
# -------------------

dewy.buildmat.TSF3.siteC <- function(density){

 # Create the matrix
 TSF3.mat.siteC = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF3.mat.siteC) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF3.mat.siteC) = colnames(TSF3.mat.siteC)
 
 # In/out seedbank (SB)
 TSF3.mat.siteC["SB", "SB"] = seed.bankSiteC$staySB[seed.bankSiteC$TSF == "three"]
 TSF3.mat.siteC["SD", "SB"] = seed.bankSiteC$outSB[seed.bankSiteC$TSF == "three"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF3.mat.siteC["J", "SD"] = vr.per.TSF.per.site$survSD[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"]
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF3.mat.siteC["SR", "J"] = vr.per.TSF.per.site$survJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * (1 - vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"])
 TSF3.mat.siteC["LR", "J"] = vr.per.TSF.per.site$survJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
 TSF3.mat.siteC["SR", "SR"] = predict(survSR.siteC, newdata = data.frame(TSF = "three", density = density), type = "response") * (1 - vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"])
 TSF3.mat.siteC["LR", "SR"] = predict(survSR.siteC, newdata = data.frame(TSF = "three", density = density), type = "response") * vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"]
 
 # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
 TSF3.mat.siteC["LR", "LR"] = predict(transLR.siteC.TSF3, newdata = data.frame(density = density), type = "response") * (1 - predict(transLR.siteC.TSF3, newdata = data.frame(density = density), type = "response"))
 TSF3.mat.siteC["SR", "LR"] = predict(transLR.siteC.TSF3, newdata = data.frame(density = density), type = "response") * predict(transLR.siteC.TSF3, newdata = data.frame(density = density), type = "response")
 
 # Reproduction of Small Reproductive (SR)
 TSF3.mat.siteC["SB", "SR"] = predict(floweringprobSR.siteC, newdata = data.frame(TSF = "three"), type = "response") * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * seeds.per.flower * seed.bankSiteC$goSB[seed.bankSiteC$TSF == "three"]
 TSF3.mat.siteC["SD", "SR"] = predict(floweringprobSR.siteC, newdata = data.frame(TSF = "three"), type = "response") * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * seeds.per.flower * seed.bankSiteC$goCont[seed.bankSiteC$TSF == "three"]
 
 # Reproduction of Large Reproductive (LR)
 TSF3.mat.siteC["SB", "LR"] = unique(predict(floweringprobLR.siteC.TSF3, type = "response")) * unique(predict(nbfsLR.siteC.TSF3, type = "response")) * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * seeds.per.flower * seed.bankSiteC$goSB[seed.bankSiteC$TSF == "three"]
 TSF3.mat.siteC["SD", "LR"] = unique(predict(floweringprobLR.siteC.TSF3, type = "response")) * unique(predict(nbfsLR.siteC.TSF3, type = "response")) * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteC" & vr.per.TSF.per.site$TSF == "TSF3"] * seeds.per.flower * seed.bankSiteC$goCont[seed.bankSiteC$TSF == "three"]
 
 return(TSF3.mat.siteC)
}


## 2.5. Site C - TSF>3 ----
# --------------------

dewy.buildmat.stochastic.siteC <- function(density, year){
  
  # Creating the matrix
  TSF3stoch.mat.siteC = matrix(0, nrow = 5, ncol = 5)
  colnames(TSF3stoch.mat.siteC) = c("SB", "SD", "J", "SR", "LR")
  rownames(TSF3stoch.mat.siteC) = colnames(TSF3stoch.mat.siteC)
  
  # In/out seedbank (SB)
  TSF3stoch.mat.siteC["SB", "SB"] = seed.bankSiteC$staySB[seed.bankSiteC$TSF == "three"]
  TSF3stoch.mat.siteC["SD", "SB"] = seed.bankSiteC$outSB[seed.bankSiteC$TSF == "three"]
  
  # Seedlings (SD) to Juveniles (J)
  TSF3stoch.mat.siteC["J", "SD"] = predict(survSD.HG, newdata = data.frame(density = density, time = year), type = "response")
  
  # Juveniles (J) to Reproductive (SR/LR)
  TSF3stoch.mat.siteC["SR", "J"] = predict(survJ.HG, newdata = data.frame(time = year), type = "response") * (1 - predict(transJ.HG, newdata = data.frame(density = density, time = year), type = "response"))
  TSF3stoch.mat.siteC["LR", "J"] = predict(survJ.HG, newdata = data.frame(time = year), type = "response") * predict(transJ.HG, newdata = data.frame(density = density, time = year), type = "response")
  
  # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
  TSF3stoch.mat.siteC["SR", "SR"] = predict(survSR.HG, newdata = data.frame(density = density, time = year), type = "response") * (1 - predict(transSR.HG, newdata = data.frame(time = year), type = "response"))
  TSF3stoch.mat.siteC["LR", "SR"] = predict(survSR.HG, newdata = data.frame(density = density, time = year), type = "response") * predict(transSR.HG, newdata = data.frame(time = year), type = "response")
  
  # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
  TSF3stoch.mat.siteC["LR", "LR"] = predict(survLR.HG, newdata = data.frame(time = year), type = "response") * (1 - predict(transLR.HG, newdata = data.frame(density = density), type = "response"))
  TSF3stoch.mat.siteC["SR", "LR"] = predict(survLR.HG, newdata = data.frame(time = year), type = "response") * predict(transLR.HG, newdata = data.frame(density = density), type = "response")
  
  # Reproduction of Small Reproductive (SR)
  TSF3stoch.mat.siteC["SB", "SR"] = predict(floweringprobSR.HG, newdata = data.frame(time = year), type = "response") * vr.per.year.TSF3plus$nbfsSR[vr.per.year.TSF3plus$grazing == "high" & vr.per.year.TSF3plus$year == year] * predict(nbfpsSR.HG, newdata = data.frame(density = density, time = year), type = "response") * seeds.per.flower * seed.bankSiteC$goSB[seed.bankSiteC$TSF == "three"]
  TSF3stoch.mat.siteC["SD", "SR"] = predict(floweringprobSR.HG, newdata = data.frame(time = year), type = "response") * vr.per.year.TSF3plus$nbfsSR[vr.per.year.TSF3plus$grazing == "high" & vr.per.year.TSF3plus$year == year] * predict(nbfpsSR.HG, newdata = data.frame(density = density, time = year), type = "response") * seeds.per.flower * seed.bankSiteC$goCont[seed.bankSiteC$TSF == "three"]
  
  # Reproduction of Large Reproductive (LR)
  TSF3stoch.mat.siteC["SB", "LR"] = predict(floweringprobLR.HG, newdata = data.frame(density = density, time = year), type = "response") * vr.per.year.TSF3plus$nbfsLR[vr.per.year.TSF3plus$grazing == "high" & vr.per.year.TSF3plus$year == year] * predict(nbfpsLR.HG, newdata = data.frame(time = year), type = "response") * seeds.per.flower * seed.bankSiteC$goSB[seed.bankSiteC$TSF == "three"]
  TSF3stoch.mat.siteC["SD", "LR"] = predict(floweringprobLR.HG, newdata = data.frame(density = density, time = year), type = "response") * vr.per.year.TSF3plus$nbfsLR[vr.per.year.TSF3plus$grazing == "high" & vr.per.year.TSF3plus$year == year] * predict(nbfpsLR.HG, newdata = data.frame(time = year), type = "response") * seeds.per.flower * seed.bankSiteC$goCont[seed.bankSiteC$TSF == "three"]
  
  return(TSF3stoch.mat.siteC)
}



## 2.6. Site E - TSF0 ----
# -------------------

dewy.buildmat.TSF0.siteE <- function(){
  
 # Creating the matrix
 TSF0.mat.siteE = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF0.mat.siteE) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF0.mat.siteE) = colnames(TSF0.mat.siteE)
 
 # Stay in the seed bank
 TSF0.mat.siteE["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "zero"]
 
 # Germinating to seedling
 TSF0.mat.siteE["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "zero"] * sb.transition.TSF0$transSB_SD[sb.transition.TSF0$site == "siteE"]
 
 # Germinating and growing to juvenile
 TSF0.mat.siteE["J", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "zero"] * sb.transition.TSF0$transSB_J[sb.transition.TSF0$site == "siteE"]
 
 return(TSF0.mat.siteE)
}


## 2.7. Site E - TSF1 ----
# -------------------

dewy.buildmat.TSF1.siteE <- function(density){
  
 # Creating the matrix
 TSF1.mat.siteE = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF1.mat.siteE) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF1.mat.siteE) = colnames(TSF1.mat.siteE)
 
 # In/out seedbank (SB)
 TSF1.mat.siteE["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "one"]
 TSF1.mat.siteE["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "one"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF1.mat.siteE["J", "SD"] = predict(survSD.siteE.TSF1_3, newdata = data.frame(TSF = "one", density = density), type = "response")
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF1.mat.siteE["SR", "J"] = predict(survJ.siteE, newdata = data.frame(TSF = "one", density = density), type = "response") * (1 - predict(transJ.siteE, newdata = data.frame(TSF = "one", density = density), type = "response"))
 TSF1.mat.siteE["LR", "J"] = predict(survJ.siteE, newdata = data.frame(TSF = "one", density = density), type = "response") * predict(transJ.siteE, newdata = data.frame(TSF = "one", density = density), type = "response")
 
 return(TSF1.mat.siteE)
}


## 2.8. Site E - TSF2 ----
# -------------------

dewy.buildmat.TSF2.siteE <- function(density){

 # Creating the matrix
 TSF2.mat.siteE = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF2.mat.siteE) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF2.mat.siteE) = colnames(TSF2.mat.siteE)
 
 # In/out seedbank (SB)
 TSF2.mat.siteE["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "two"]
 TSF2.mat.siteE["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "two"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF2.mat.siteE["J", "SD"] = vr.per.TSF.per.site$survSD[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF2.mat.siteE["SR", "J"] = predict(survJ.siteE, newdata = data.frame(TSF = "two", density = density), type = "response") * (1 - predict(transJ.siteE, newdata = data.frame(TSF = "two", density = density), type = "response"))
 TSF2.mat.siteE["LR", "J"] = predict(survJ.siteE, newdata = data.frame(TSF = "two", density = density), type = "response") * predict(transJ.siteE, newdata = data.frame(TSF = "two", density = density), type = "response")
 
 # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
 TSF2.mat.siteE["SR", "SR"] = vr.per.TSF.per.site$survSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * (1 - vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteE["LR", "SR"] = vr.per.TSF.per.site$survSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
 TSF2.mat.siteE["LR", "LR"] = vr.per.TSF.per.site$survLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * (1 - vr.per.TSF.per.site$transLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteE["SR", "LR"] = vr.per.TSF.per.site$survLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$transLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Reproduction of Small Reproductive (SR)
 TSF2.mat.siteE["SB", "SR"] = vr.per.TSF.per.site$floweringprobSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "two"]
 TSF2.mat.siteE["SD", "SR"] = vr.per.TSF.per.site$floweringprobSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "two"]
 
 # Reproduction of Large Reproductive (LR)
 TSF2.mat.siteE["SB", "LR"] = vr.per.TSF.per.site$floweringprobLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "two"]
 TSF2.mat.siteE["SD", "LR"] = vr.per.TSF.per.site$floweringprobLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteE" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "two"]
 
 return(TSF2.mat.siteE)
}


## 2.9. Site E - TSF3 ----
# -------------------

dewy.buildmat.TSF3.siteE <- function(density){
  
 # Creating the matrix
 TSF3.mat.siteE = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF3.mat.siteE) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF3.mat.siteE) = colnames(TSF3.mat.siteE)
 
 # In/out seedbank (SB)
 TSF3.mat.siteE["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "three"]
 TSF3.mat.siteE["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "three"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF3.mat.siteE["J", "SD"] = predict(survSD.siteE.TSF1_3, newdata = data.frame(TSF = "three", density = density), type = "response")
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF3.mat.siteE["SR", "J"] = predict(survJ.siteE, newdata = data.frame(TSF = "three", density = density), type = "response") * (1 - predict(transJ.siteE, newdata = data.frame(TSF = "three", density = density), type = "response"))
 TSF3.mat.siteE["LR", "J"] = predict(survJ.siteE, newdata = data.frame(TSF = "three", density = density), type = "response") * predict(transJ.siteE, newdata = data.frame(TSF = "three", density = density), type = "response")
 
 # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
 TSF3.mat.siteE["SR", "SR"] = predict(survSR.siteE.TSF3, newdata = data.frame(density = density), type = "response") * (1 - unique(predict(transSR.siteE.TSF3, type = "response")))
 TSF3.mat.siteE["LR", "SR"] = predict(survSR.siteE.TSF3, newdata = data.frame(density = density), type = "response") * unique(predict(transSR.siteE.TSF3, type = "response"))
 
 # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
 TSF3.mat.siteE["LR", "LR"] = predict(survLR.siteE.TSF3, newdata = data.frame(density = density), type = "response") * (1 - predict(transLR.siteE.TSF3, newdata = data.frame(density), type = "response"))
 TSF3.mat.siteE["SR", "LR"] = predict(survLR.siteE.TSF3, newdata = data.frame(density = density), type = "response") * predict(transLR.siteE.TSF3, newdata = data.frame(density = density), type = "response")
 
 # Reproduction of Small Reproductive (SR)
 TSF3.mat.siteE["SB", "SR"] = unique(predict(floweringprobSR.siteE.TSF3, type = "response")) * unique(predict(nbfsSR.siteE.TSF3, type = "response")) * unique(predict(nbfpsSR.siteE.TSF3, type = "response")) * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "three"]
 TSF3.mat.siteE["SD", "SR"] = unique(predict(floweringprobSR.siteE.TSF3, type = "response")) * unique(predict(nbfsSR.siteE.TSF3, type = "response")) * unique(predict(nbfpsSR.siteE.TSF3, type = "response")) * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "three"]
 
 # Reproduction of Large Reproductive (LR)
 TSF3.mat.siteE["SB", "LR"] = unique(predict(floweringprobLR.siteE.TSF3, type = "response")) * unique(predict(nbfsLR.siteE.TSF3, type = "response")) * unique(predict(nbfpsLR.siteE.TSF3, type = "response")) * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "three"]
 TSF3.mat.siteE["SD", "LR"] = unique(predict(floweringprobLR.siteE.TSF3, type = "response")) * unique(predict(nbfsLR.siteE.TSF3, type = "response")) * unique(predict(nbfpsLR.siteE.TSF3, type = "response")) * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "three"]
 
 return(TSF3.mat.siteE)
}


## 2.10. Site F - TSF 0 ----
# --------------------

dewy.buildmat.TSF0.siteF <- function(){
  
 # Creating the matrix
 TSF0.mat.siteF = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF0.mat.siteF) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF0.mat.siteF) = colnames(TSF0.mat.siteF)
 
 # Stay in the seed bank
 TSF0.mat.siteF["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "zero"]
 
 # Germinating to seedling
 TSF0.mat.siteF["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "zero"] * sb.transition.TSF0$transSB_SD[sb.transition.TSF0$site == "siteF"]
 
 # Germinating and growing to juvenile
 TSF0.mat.siteF["J", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "zero"] * sb.transition.TSF0$transSB_J[sb.transition.TSF0$site == "siteF"]
 
 return(TSF0.mat.siteF)
}


## 2.11. Site F - TSF1 ----
# --------------------

dewy.buildmat.TSF1.siteF <- function(density){
  
 # Creating the matrix
 TSF1.mat.siteF = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF1.mat.siteF) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF1.mat.siteF) = colnames(TSF1.mat.siteF)
 
 # In/out seedbank (SB)
 TSF1.mat.siteF["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "one"]
 TSF1.mat.siteF["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "one"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF1.mat.siteF["J", "SD"] = vr.per.TSF.per.site$survSD[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF1"]
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF1.mat.siteF["SR", "J"] = unique(predict(survJ.siteF.TSF1_3, type = "response")) * (1 - unique(predict(transJ.siteF.TSF1, type = "response")))
 TSF1.mat.siteF["LR", "J"] = unique(predict(survJ.siteF.TSF1_3, type = "response")) * unique(predict(transJ.siteF.TSF1, type = "response"))
 
 return(TSF1.mat.siteF)
}


## 2.12. Site F - TSF2 ----
# --------------------

dewy.buildmat.TSF2.siteF <- function(){
  
 # Creating the matrix
 TSF2.mat.siteF = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF2.mat.siteF) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF2.mat.siteF) = colnames(TSF2.mat.siteF)
 
 # In/out seedbank (SB)
 TSF2.mat.siteF["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "two"]
 TSF2.mat.siteF["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "two"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF2.mat.siteF["J", "SD"] = vr.per.TSF.per.site$survSD[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF2.mat.siteF["SR", "J"] = vr.per.TSF.per.site$survJ[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * (1 - vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteF["LR", "J"] = vr.per.TSF.per.site$survJ[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
 TSF2.mat.siteF["SR", "SR"] = vr.per.TSF.per.site$survSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * (1 - vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"])
 TSF2.mat.siteF["LR", "SR"] = vr.per.TSF.per.site$survSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"]
 
 # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
 TSF2.mat.siteF["LR", "LR"] = predict(survLR.siteF, newdata = data.frame(TSF = "two"), type = "response") * (1 - unique(predict(transLR.siteF.TSF2, type = "response")))
 TSF2.mat.siteF["SR", "LR"] = predict(survLR.siteF, newdata = data.frame(TSF = "two"), type = "response") * unique(predict(transLR.siteF.TSF2, type = "response"))

 # Reproduction of Small Reproductive (SR)
 TSF2.mat.siteF["SB", "SR"] = vr.per.TSF.per.site$floweringprobSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "two"]
 TSF2.mat.siteF["SD", "SR"] = vr.per.TSF.per.site$floweringprobSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "two"]
 
 # Reproduction of Large Reproductive (LR)
 TSF2.mat.siteF["SB", "LR"] = unique(predict(floweringprobLR.siteF, type = "response")) * vr.per.TSF.per.site$nbfsLR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "two"]
 TSF2.mat.siteF["SD", "LR"] = unique(predict(floweringprobLR.siteF, type = "response")) * vr.per.TSF.per.site$nbfsLR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * vr.per.TSF.per.site$nbfpsLR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF2"] * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "two"]
 
 return(TSF2.mat.siteF)
}


## 2.13. Site F - TSF3 ----
# --------------------

dewy.buildmat.TSF3.siteF <- function(density){
  
 # Creating the matrix
 TSF3.mat.siteF = matrix(0, nrow = 5, ncol = 5)
 colnames(TSF3.mat.siteF) = c("SB", "SD", "J", "SR", "LR")
 rownames(TSF3.mat.siteF) = colnames(TSF3.mat.siteF)
 
 # In/out seedbank (SB)
 TSF3.mat.siteF["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "three"]
 TSF3.mat.siteF["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "three"]
 
 # Seedlings (SD) to Juveniles (J)
 TSF3.mat.siteF["J", "SD"] = predict(survSD.siteF.TSF3, newdata = data.frame(density = density), type = "response")
 
 # Juveniles (J) to Reproductive (SR/LR)
 TSF3.mat.siteF["SR", "J"] = unique(predict(survJ.siteF.TSF1_3, type = "response")) * (1 - vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"])
 TSF3.mat.siteF["LR", "J"] = unique(predict(survJ.siteF.TSF1_3, type = "response")) * vr.per.TSF.per.site$transJ[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"]
 
 # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
 TSF3.mat.siteF["SR", "SR"] = vr.per.TSF.per.site$survSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * (1 - vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"])
 TSF3.mat.siteF["LR", "SR"] = vr.per.TSF.per.site$survSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$transSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"]
 
 # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
 TSF3.mat.siteF["LR", "LR"] = predict(survLR.siteF, newdata = data.frame(TSF = "three"), type = "response") * (1 - vr.per.TSF.per.site$transLR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"])
 TSF3.mat.siteF["SR", "LR"] = predict(survLR.siteF, newdata = data.frame(TSF = "three"), type = "response") * vr.per.TSF.per.site$transLR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"]
 
 # Reproduction of Small Reproductive (SR)
 TSF3.mat.siteF["SB", "SR"] = vr.per.TSF.per.site$floweringprobSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "three"]
 TSF3.mat.siteF["SD", "SR"] = vr.per.TSF.per.site$floweringprobSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$nbfsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * vr.per.TSF.per.site$nbfpsSR[vr.per.TSF.per.site$site == "siteF" & vr.per.TSF.per.site$TSF == "TSF3"] * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "three"]
 
 # Reproduction of Large Reproductive (LR)
 TSF3.mat.siteF["SB", "LR"] = unique(predict(floweringprobLR.siteF, type = "response")) * predict(nbfsLR.siteF.TSF3, newdata = data.frame(density = density), type = "response") * unique(predict(nbfpsLR.siteF.TSF3, type = "response")) * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "three"]
 TSF3.mat.siteF["SD", "LR"] = unique(predict(floweringprobLR.siteF, type = "response")) * predict(nbfsLR.siteF.TSF3, newdata = data.frame(density = density), type = "response") * unique(predict(nbfpsLR.siteF.TSF3, type = "response")) * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "three"]
 
 return(TSF3.mat.siteF)
}


## 2.14. Sites E and F - TSF>3 ----
# ----------------------------

dewy.buildmat.stochastic.sitesEF <- function(density, year){
  
  # Creating the matrix
  TSF3stoch.mat.siteEF = matrix(0, nrow = 5, ncol = 5)
  colnames(TSF3stoch.mat.siteEF) = c("SB", "SD", "J", "SR", "LR")
  rownames(TSF3stoch.mat.siteEF) = colnames(TSF3stoch.mat.siteEF)
  
  # In/out seedbank (SB)
  TSF3stoch.mat.siteEF["SB", "SB"] = seed.bankSitesEF$staySB[seed.bankSitesEF$TSF == "three"]
  TSF3stoch.mat.siteEF["SD", "SB"] = seed.bankSitesEF$outSB[seed.bankSitesEF$TSF == "three"]
  
  # Seedlings (SD) to Juveniles (J)
  TSF3stoch.mat.siteEF["J", "SD"] = predict(survSD.LG, newdata = data.frame(density = density, time = year), type = "response")
  
  # Juveniles (J) to Reproductive (SR/LR)
  TSF3stoch.mat.siteEF["SR", "J"] = predict(survJ.LG, newdata = data.frame(density = density, time = year), type = "response") * (1 - predict(transJ.LG, newdata = data.frame(time = year), type = "response"))
  TSF3stoch.mat.siteEF["LR", "J"] = predict(survJ.LG, newdata = data.frame(density = density, time = year), type = "response") * predict(transJ.LG, newdata = data.frame(time = year), type = "response")
  
  # Small reproductive (SR) to Large Reproductive (LR) or stasis (SR)
  TSF3stoch.mat.siteEF["SR", "SR"] = predict(survSR.LG, newdata = data.frame(time = year), type = "response") * (1 - vr.per.year.TSF3plus$transSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == year])
  TSF3stoch.mat.siteEF["LR", "SR"] = predict(survSR.LG, newdata = data.frame(time = year), type = "response") * vr.per.year.TSF3plus$transSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == year]
  
  # Large reproductive (LR) to Small Reproductive (SR) or stasis (LR)
  TSF3stoch.mat.siteEF["LR", "LR"] = predict(survLR.LG, newdata = data.frame(time = year), type = "response") * (1 - predict(transLR.LG, newdata = data.frame(time = year), type = "response"))
  TSF3stoch.mat.siteEF["SR", "LR"] = predict(survLR.LG, newdata = data.frame(time = year), type = "response") * predict(transLR.LG, newdata = data.frame(time = year), type = "response")
  
  # Reproduction of Small Reproductive (SR)
  TSF3stoch.mat.siteEF["SB", "SR"] = ifelse(year != 2016, predict(floweringprobSR.LG, newdata = data.frame(density = density, time = year), type = "response"), mean(predict(floweringprobSR.LG, newdata = data.frame(density = density, time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * vr.per.year.TSF3plus$nbfsSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == year] * vr.per.year.TSF3plus$nbfpsSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == year] * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "three"]
  TSF3stoch.mat.siteEF["SD", "SR"] = ifelse(year != 2016, predict(floweringprobSR.LG, newdata = data.frame(density = density, time = year), type = "response"), mean(predict(floweringprobSR.LG, newdata = data.frame(density = density, time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * vr.per.year.TSF3plus$nbfsSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == year] * vr.per.year.TSF3plus$nbfpsSR[vr.per.year.TSF3plus$grazing == "low" & vr.per.year.TSF3plus$year == year] * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "three"]
  
  # Reproduction of Large Reproductive (LR)
  TSF3stoch.mat.siteEF["SB", "LR"] = ifelse(year != 2016, predict(floweringprobLR.LG, newdata = data.frame(density = density, time = year), type = "response"), mean(predict(floweringprobLR.LG, newdata = data.frame(density = density, time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * ifelse(year != 2016, predict(nbfsLR.LG, newdata = data.frame(time = year), type = "response"), mean(predict(nbfsLR.LG, newdata = data.frame(time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * ifelse(year != 2016, predict(nbfpsLR.LG, newdata = data.frame(time = year), type = "response"), mean(predict(nbfpsLR.LG, newdata = data.frame(time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * seeds.per.flower * seed.bankSitesEF$goSB[seed.bankSitesEF$TSF == "three"]
  TSF3stoch.mat.siteEF["SD", "LR"] = ifelse(year != 2016, predict(floweringprobLR.LG, newdata = data.frame(density = density, time = year), type = "response"), mean(predict(floweringprobLR.LG, newdata = data.frame(density = density, time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * ifelse(year != 2016, predict(nbfsLR.LG, newdata = data.frame(time = year), type = "response"), mean(predict(nbfsLR.LG, newdata = data.frame(time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * ifelse(year != 2016, predict(nbfpsLR.LG, newdata = data.frame(time = year), type = "response"), mean(predict(nbfpsLR.LG, newdata = data.frame(time = vr.per.year.TSF3plus$year[- which(vr.per.year.TSF3plus$year == 2016)]), type = "response"))) * seeds.per.flower * seed.bankSitesEF$goCont[seed.bankSitesEF$TSF == "three"]
  
  return(TSF3stoch.mat.siteEF)
}
