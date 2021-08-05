############################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script uses the data of a dewy pine population in South-Eastern Spain collected between 2012 and 2019. 
# The aim of this script is to define vital rates for the seed bank, number of seeds per flower, and compute seed bank transition rates to seedlings or juveniles for TSF0.
#
# Author: Eva Conquet
#
###########################################################################


###########################################################################
#
# 1. House keeping and loading libraries and data ----
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading and preparing data ----
# --------------------------------

data.droso = read.csv("DewyPinesTSF1_to_3_Data.csv")
head(data.droso)




###########################################################################
#
# 2. Seed bank and number of seeds per flower constants ----
#
###########################################################################

## 2.1. Seed bank constants ----
# -------------------------


# Site C

seed.bankSiteC = data.frame(TSF = c("zero", "one", "two", "three"),
                          goSB = c(0, 0, 0.93, 0.84),
                          goCont = c(0, 0, 0.04, 0.13),
                          staySB = c(0.1, 0.05, 0.60, 0.60),
                          outSB = c(0.81 * 0.45, 0.054 * 0.45, 0.03 * 0.35, 0.05 * 0.35))


# Site E and F

seed.bankSitesEF = data.frame(TSF = c("zero", "one", "two", "three"),
                            goSB = c(0, 0, 0.93, 0.93),
                            goCont = c(0, 0, 0.04, 0.04),
                            staySB = c(0.1, 0.05, 0.85, 0.85),
                            outSB = c(0.81 * 0.84, 0.054 * 0.84, 0.03 * 0.82, 0.03 * 0.8))


## 2.2. Number of seeds per flower ----
# --------------------------------

seeds.per.flower = 9.8




###########################################################################
#
# 3. TSF0 transition rates ----
#
###########################################################################

## 3.1. Site C ----
# ------------

propJ.siteC = nrow(data.droso[data.droso$site == "siteC" & data.droso$TSF == "one" & data.droso$stage == "J", ]) / nrow(data.droso[data.droso$site == "siteC" & data.droso$TSF == "one", ])

propSD.siteC = nrow(data.droso[data.droso$site == "siteC" & data.droso$TSF == "one" & data.droso$stage == "SD", ]) / nrow(data.droso[data.droso$site == "siteC" & data.droso$TSF == "one", ])


## 3.2. Site E ----
# ------------

propJ.siteE = nrow(data.droso[data.droso$site == "siteE" & data.droso$TSF == "one" & data.droso$stage == "J", ]) / nrow(data.droso[data.droso$site == "siteE" & data.droso$TSF == "one", ])

propSD.siteE = nrow(data.droso[data.droso$site == "siteE" & data.droso$TSF == "one" & data.droso$stage == "SD", ]) / nrow(data.droso[data.droso$site == "siteE" & data.droso$TSF == "one", ])


## 3.3. Site F ----
# ------------

propJ.siteF = nrow(data.droso[data.droso$site == "siteF" & data.droso$TSF == "one" & data.droso$stage == "J", ]) / nrow(data.droso[data.droso$site == "siteF" & data.droso$TSF == "one", ])

propSD.siteF = nrow(data.droso[data.droso$site == "siteF" & data.droso$TSF == "one" & data.droso$stage == "SD", ]) / nrow(data.droso[data.droso$site == "siteF" & data.droso$TSF == "one", ])


# Creating a table for all sites

sb.transition.TSF0 = data.frame(site = c("siteC", "siteE", "siteF"),
                                transSB_SD = c(propSD.siteC, propSD.siteE, propSD.siteF), 
                                transSB_J = c(propJ.siteC, propJ.siteE, propJ.siteF))




###########################################################################
#
# 4. Saving tables ----
#
###########################################################################

write.csv(seed.bankSiteC, "DPVR_Seedbank_SiteC.csv", row.names = F)
write.csv(seed.bankSitesEF, "DPVR_Seedbank_SitesEF.csv", row.names = F)
write.csv(sb.transition.TSF0, "DPVR_SeedbankTransitions_TSF0.csv", row.names = F)
save(seeds.per.flower, file = "NbSeedsPerFlower.RData")
