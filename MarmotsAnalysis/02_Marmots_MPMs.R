##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
#
# This script contains the functions needed to build the matrix population models 
# for the marmot population.
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
}

load.librairies()


## 1.3. Loading data ----
# ------------------

marmots.vr = read.csv("Marmots_VitalRatesEstimates.csv")




###########################################################################
#
# 2. Functions to build summer and winter transition matrices ----
#
###########################################################################


## 3.1. Summer ----
# ------------

marmots.buildmat.summer <- function(marmots.vr, year){
  
  # Creating the matrix
  summer.mat = matrix(0, nrow = 4, ncol = 4)
  colnames(summer.mat) = c("J", "Y", "NR", "R")
  rownames(summer.mat) = colnames(summer.mat)
  
  # Yearlings
  summer.mat["Y", "Y"] = marmots.vr$survY_summer[which(marmots.vr$year == year)]
  
  # Non reproductive adults
  summer.mat["NR", "NR"] = marmots.vr$survNR_summer[which(marmots.vr$year == year)]
  
  # Reproductive adults
  summer.mat["R", "R"] = marmots.vr$survR_summer[which(marmots.vr$year == year)]
  
  # Reproduction of reproductive adults
  summer.mat["J", "R"] = (marmots.vr$survR_summer[which(marmots.vr$year == year)] * marmots.vr$recruit[which(marmots.vr$year == year)]) / 2
  
  return(summer.mat)
}


## 3.2. Winter ----
# ------------

marmots.buildmat.winter <- function(marmots.vr, year){
  
  # Creating the matrix
  winter.mat = matrix(0, nrow = 4, ncol = 4)
  colnames(winter.mat) = c("J", "Y", "NRA", "RA")
  rownames(winter.mat) = colnames(winter.mat)
  
  # Juveniles
  winter.mat["Y", "J"] = marmots.vr$survJ[which(marmots.vr$year == year)]
  
  # Yearlings
  winter.mat["RA", "Y"] = marmots.vr$survY_winter[which(marmots.vr$year == year)] * marmots.vr$transYR[which(marmots.vr$year == year)] 
  winter.mat["NRA", "Y"] = marmots.vr$survY_winter[which(marmots.vr$year == year)] * (1 - marmots.vr$transYR[which(marmots.vr$year == year)])
  
  # Non reproductive adults
  winter.mat["NRA", "NRA"] = marmots.vr$survNR_winter[which(marmots.vr$year == year)] * (1 - marmots.vr$transNRR[which(marmots.vr$year == year)])
  winter.mat["RA", "NRA"] = marmots.vr$survNR_winter[which(marmots.vr$year == year)] * marmots.vr$transNRR[which(marmots.vr$year == year)]
  
  # Reproductive adults
  winter.mat["RA", "RA"] = marmots.vr$survR_winter[which(marmots.vr$year == year)] * marmots.vr$transRR[which(marmots.vr$year == year)]
  winter.mat["NRA" ,"RA"] =  marmots.vr$survR_winter[which(marmots.vr$year == year)] * (1 - marmots.vr$transRR[which(marmots.vr$year == year)])
  
  return(winter.mat)
}