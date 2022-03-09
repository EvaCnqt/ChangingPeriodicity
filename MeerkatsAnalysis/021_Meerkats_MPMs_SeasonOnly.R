##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script contains the functions needed to build the matrix population models 
# for the meerkat population, without the density effect.
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

load("GLMM_survJ_SeasonOnly.RData")
load("GLMM_survS_SeasonOnly.RData")
load("GLMM_survH_SeasonOnly.RData")
load("GLMM_survD_SeasonOnly.RData")
load("GLMM_transHD_SeasonOnly.RData")
load("GLMM_emig_SeasonOnly.RData")
load("GLMM_recruitH_SeasonOnly.RData")
load("GLMM_recruitD_SeasonOnly.RData")

predict.recruitH = function(year, season){
  
  coefs.recruitment.H = coef(recruitment.H)
  recruitH.all.years = data.frame(year = rownames(coefs.recruitment.H), pred = NA)
  
  if(season == "dry"){
    
    recruitH.all.years$pred = exp(coefs.recruitment.H$`(Intercept)`)
    
  }

  else{
    
    recruitH.all.years$pred = exp(coefs.recruitment.H$`(Intercept)` + coefs.recruitment.H$seasonrain)
    
  }
  
  recruitH = recruitH.all.years$pred[which(recruitH.all.years$year == year)]
    
  return(recruitH)
  
}


predict.recruitD = function(year, season){
  
  coefs.recruitment.D = coef(recruitment.D)
  recruitD.all.years = data.frame(year = rownames(coefs.recruitment.D), pred = NA)
  
  if(season == "dry"){
    
    recruitD.all.years$pred = exp(coefs.recruitment.D$`(Intercept)`)

  }
  
  else{
    
    recruitD.all.years$pred = exp(coefs.recruitment.D$`(Intercept)` + coefs.recruitment.D$seasonrain)
    
  }
  
  recruitD = recruitD.all.years$pred[which(recruitD.all.years$year == year)]
  
  return(recruitD)
  
}




###########################################################################
#
# 2. Functions to create seasonal matrices ----
#
###########################################################################

meerkats.buildmat <- function(year, season){
  
  # Creating the matrix
  seas.mat = matrix(0, nrow =  4, ncol = 4)
  colnames(seas.mat) = c("J", "S", "H", "D")
  rownames(seas.mat) = colnames(seas.mat)
  
  # Juveniles
  seas.mat["J", "H"] = predict.recruitH(year = year,
                                        season = season)
  
  seas.mat["J", "D"] = predict.recruitD(year = year,
                                        season = season)
  
  # Subadults
  seas.mat["S", "J"] = predict(survJ, newdata = data.frame(year = year, 
                                                           season = season),
                               type = "response")
  
  # Helpers
  seas.mat["H", "S"] = predict(survS, newdata = data.frame(year = year,
                                                           season = season),
                               type = "response")
  
  seas.mat["H", "H"] = predict(survH, newdata = data.frame(year = year,
                                                           season = season),
                               type = "response") * (1 - predict(emig, newdata = data.frame(year = year,
                                                                                            season = season),
                                                                 type = "response")) * (1 - predict(transition, newdata = data.frame(year = year, 
                                                                                                                                     season = season),
                                                                                                    type = "response"))
  
  # Dominants 
  seas.mat["D", "H"] = predict(survH, newdata = data.frame(year = year,
                                                           season = season),
                               type = "response") * (1 - predict(emig, newdata = data.frame(year = year,
                                                                                            season = season),
                                                                 type = "response")) * predict(transition, newdata = data.frame(year = year,
                                                                                                                                season = season), 
                                                                                               type = "response")
  seas.mat["D", "D"] = predict(survD, newdata = data.frame(year = year,
                                                           season = season), 
                               type = "response")
  
  return(seas.mat)
}
