##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity 
# (Conquet et al., under review at Ecology).
#
# This script contains the function needed to project the dewy pine population dynamics with or
# without perturbed vital rates in each post-fire habitat state (time since fire, TSF) and under
# each fire regime: stochastic or periodic fires occurring every 15 or 30 years. 
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

source("04_DewyPines_MPMs.R")
mean.densities = read.csv("MeanDensities.csv")
data.grazing = read.csv("DewyPinesTSF3plus_Data.csv")
nbsquares = read.csv("NbSquares_PerSites_TSF1_3.csv")

# Getting the number of squares
nbsquares.siteE.TSF1 = nbsquares$nbsquares[nbsquares$site == "siteE" & nbsquares$TSF == "one"]
nbsquares.siteE.TSF2 = nbsquares$nbsquares[nbsquares$site == "siteE" & nbsquares$TSF == "two"]
nbsquares.siteE.TSF3 = nbsquares$nbsquares[nbsquares$site == "siteE" & nbsquares$TSF == "three"]

nbsquares.siteF.TSF1 = nbsquares$nbsquares[nbsquares$site == "siteF" & nbsquares$TSF == "one"]
nbsquares.siteF.TSF2 = nbsquares$nbsquares[nbsquares$site == "siteF" & nbsquares$TSF == "two"]
nbsquares.siteF.TSF3 = nbsquares$nbsquares[nbsquares$site == "siteF" & nbsquares$TSF == "three"]




###########################################################################
#
# 2. Projection function ----
#
###########################################################################

stoch.sim.droso <- function(nsimul = 500, 
                            nyears = 100, 
                            n0, 
                            disturbance.type = "stochastic",
                            fire.prob = 1/30,
                            fire.freq = 15,
                            perturbed.state = 0, 
                            density = T){
 
 # (1) Control case = Simulating what happens in the normal conditions: We use the populations in the fire-disturbed sites E and F as natural populations. Each year, we randomly pick an MPM between site E and F.
 # (2) Perturbed periodicity case = Perturbing the periodicity in post-fire habitats in each post-fire state independently through inclusion of human-disturbance conditions, while keeping fire-return probability to 0.033 (one fire every 30 years). To do so, we use the given state matrix of the human-disturbed site C.
 
 # Preparing arrays, vectors, and lists to record lambdas, extinction, density, population vectors, and post-fire habitat states
 
 log.lambda = array(NA, c(nsimul, nyears))
 ext = rep(0, nsimul)
 ext.year = rep(0, nsimul)
 dens = array(NA, c(nsimul, nyears))
 pop.vec = list()
 
 post.fire.habitat = list()
 
 # Building the Markov chain to create the succession of post-fire environmental states
 
 # If we include fire disturbance, we use fire.prob as fire probability to go back to state 1. Otherwise, the population stays in state 5 until the next perturbation.
 p = ifelse(disturbance.type == "stochastic", 1 - fire.prob, 1)

 envS = matrix(0, 5, 5)
 
 envS[2, 1] = envS[3, 2] = envS[4, 3] = envS[5, 4] = 1
 envS[1, 5] = 1 - p
 envS[5, 5] = p
 
 colnames(envS) = rownames(envS) = 1:5
 
 # Simulations
 for(i in 1:nsimul){
  
  if(nsimul >= 20){
   if((i %% round(nsimul / 20)) == 0)cat(i, "simulations\n")
  }
  
  # Starting with the same n1/dens1 for each simulation
  n1 = n0
  dens1 = 0
  
  # Building the post-fire environmental states succession vector
  states = numeric(nyears)
  env.at.t = 1
  
  # If we assume periodic fires, going back to state 1 is defined by the frequency of that disturbance. Otherwise only the first state is defined as 1, the rest depends on the fire-return probability
  if(disturbance.type == "stochastic"){
    
    states[1] = 1
    
    for(x in 2:nyears){
      states[x] = env.at.t = sample(ncol(envS), 1, pr = envS[, env.at.t])
    }
  }
  
  else if(disturbance.type == "periodic"){
    
    states[seq(1, nyears, fire.freq)] = 1
    states[which(states == 1) + 1] = 2
    states[which(states == 1) + 2] = 3
    states[which(states == 1) + 3] = 4
    states[which(states == 0)] = 5
    
  }
  
  post.fire.habitat[[i]] = states
  
  # Preparing lists to record yearly population vector for each simulation
  pop.vec.sim = list()
  
  # Counter for the number of fires
  nb.fires = 0
  
  # Yearly time step
  for(j in 1:nyears){
   
   # Recording population vector at the beginning of the year to be able to compute lambda
   pop.vec1 = n1
   
   # Year state from the state vector
   state = states[j]
   
   ###########################
   
   ## TSF0 ----
   
   if(state == 1){
    
    nb.fires = nb.fires + 1 # Adding 1 to the fires counter
    
    random.mat = sample(c(1, 2), 1) # Picking a matrix between the two undisturbed sites E and F.
    
    # Build MPM for state 1 (TSF0)
    
    # Case (2) - Perturbed periodicity case 
    
    if(any(perturbed.state == 1)){
     
     # Use the MPM of the human-perturbed site C in state 1
     pop.mat = dewy.buildmat.TSF0.siteC()
    }
    
    # Case (1) - Unperturbed periodicity case
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      pop.mat = dewy.buildmat.TSF0.siteE()
     }
     
     else if(random.mat == 2){
      pop.mat = dewy.buildmat.TSF0.siteF()
     }
    }
    
    # Compute new population vector from the MPMs
    n1 = pop.mat %*% n1
    
    # Capping the seedbank to 4000 to avoid population density explosion (twice the maximum observed amount of seeds in the soil)
    if(n1[1] > 4000) n1[1] = 4000
    
    # Capping the number of seedlings to 1700 to avoid population density explosion (twice the maximum observed amount of seedlings in the population)
    if(n1[2] > 1700) n1[2] = 1700
    
    # Compute above-ground density (i.e. without seedbank)
    if(density == T){
     
     if(random.mat == 1){
      dens1 = sum(n1[-1]) / nbsquares.siteE.TSF1
     }
     
     else if(random.mat == 2){
      dens1 = sum(n1[-1]) / nbsquares.siteF.TSF1
     }
    } 
   }
   
   ###########################
   
   ## TSF1 case ----
   
   if(state == 2){
    
    random.mat = sample(c(1, 2), 1) # Picking a matrix between the two undisturbed sites E and F.
    
    # Build MPM for state 2 (TSF1)
    
    # Case (2) - Perturbed periodicity case 
    
    if(any(perturbed.state == 2)){
     
     # Use the MPM of the human-perturbed site C in state 2
     if(density == T){ # For density-dependent simulation, use the updated density
      pop.mat = dewy.buildmat.TSF1.siteC(dens1)
     }
     
     else if(density == F){ # For average density simulation, use the average density in site C, TSF1
      
      dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "one"]
      pop.mat = dewy.buildmat.TSF1.siteC(dens1)
     }
    }
    
    
    # Case (1) - Unperturbed periodicity case 
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      
      if(density == T){ # For density-dependent simulation, use the updated density
       pop.mat = dewy.buildmat.TSF1.siteE(dens1)
      }
      
      else if(density == F){# For average density simulation, use the average density in site E, TSF1
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteE" & mean.densities$TSF == "one"]
       pop.mat = dewy.buildmat.TSF1.siteE(dens1)
      }
     }
     
     else if(random.mat == 2){
      
      if(density == T){ # For density-dependent simulation, use the updated density
       pop.mat = dewy.buildmat.TSF1.siteF(dens1)
      }
      
      else if(density == F){ # For average density simulation, use the average density in site F, TSF1
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteF" & mean.densities$TSF == "one"]
       pop.mat = dewy.buildmat.TSF1.siteF(dens1)
      }
     }
    }
    
    # Compute new population vector from the matrices
    n1 = pop.mat %*% n1
    
    # Capping the seedbank to 4000 to avoid population density explosion (twice the maximum observed amount of seeds in the soil)
    if(n1[1] > 4000) n1[1] = 4000
    
    # Capping the number of seedlings to 1700 to avoid population density explosion (twice the maximum observed amount of seedlings in the population)
    if(n1[2] > 1700) n1[2] = 1700
    
    # Compute above-ground density (without seedbank)
    if(density == T){
     
     if(random.mat == 1){
      dens1 = sum(n1[-1]) / nbsquares.siteE.TSF2
     }
     
     else if(random.mat == 2){
      dens1 = sum(n1[-1]) / nbsquares.siteF.TSF2
     }
    }
   }
   
   ###########################
   
   ## TSF2 case ----
   
   if(state == 3){
    
    random.mat = sample(c(1, 2), 1) # Picking a matrix between the two undisturbed sites E and F.
    
    # Build MPM for state 3 (TSF2)
    
    # Case (2) - Perturbed periodicity case 
    
    if(any(perturbed.state == 3)){
     
     # Use the MPM of the human-perturbed site C in state 3
     if(density == T){ # For density-dependent simulation, use the updated density
 
      pop.mat = dewy.buildmat.TSF2.siteC(dens1)
     }
     
     else if(density == F){ # For average density simulation, use the average density in site C, TSF2
      
      dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "two"]
      pop.mat = dewy.buildmat.TSF2.siteC(dens1)
     }
    }
    
    
    # Case (1) - Unperturbed periodicity case 
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      
      if(density == T){ # For density-dependent simulation, use the updated density
       pop.mat = dewy.buildmat.TSF2.siteE(dens1)
      }
      
      else if(density == F){ # For average density simulation, use the average density in site E, TSF2
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteE" & mean.densities$TSF == "two"]
       pop.mat = dewy.buildmat.TSF2.siteE(dens1)
      }
     }
     
     else if(random.mat == 2){
 
      pop.mat = dewy.buildmat.TSF2.siteF()
     }
    }
    
    # Compute new population vector from the matrices
    n1 = pop.mat %*% n1
    
    # Capping the seedbank to 4000 to avoid population explosion(twice the maximum observed amount of seeds in the soil)
    if(n1[1] > 4000) n1[1] = 4000
    
    # Capping the number of seedlings to 1700 to avoid population explosion (twice the maximum observed amount of seedlings in the population)
    if(n1[2] > 1700) n1[2] = 1700
    
    # Compute above-ground density (without seedbank)
    if(density == T){
     
     if(random.mat == 1){
      
      dens1 = sum(n1[-1]) / nbsquares.siteE.TSF3
     }
     
     else if(random.mat == 2){
      
      dens1 = sum(n1[-1]) / nbsquares.siteF.TSF3
     }
    }
   }
   
   ###########################
   
   ## TSF3 case ----
   
   if(state == 4){
    
    random.mat = sample(c(1, 2), 1) # Picking a matrix between unperturbed sites E and F
    
    # Build MPM for state 4 (TSF>2)
    
    # Case (2) - Perturbed periodicity case 
    
    if(any(perturbed.state == 4)){
     
     # Use the MPM of the human-perturbed site C in state 4
     if(density == T){ # For density-dependent simulation, use the updated density
      
      pop.mat = dewy.buildmat.TSF3.siteC(dens1)
     }
     
     else if(density == F){ # For average density simulation, use the average density in site C, TSF3
      
      dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "three"]
      pop.mat = dewy.buildmat.TSF3.siteC(dens1)
     }
    } 
    
    # Case (1) - Unperturbed periodicity case
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      
      if(density == T){ # For density-dependent simulation, use the updated density
       pop.mat = dewy.buildmat.TSF3.siteE(dens1)
      }
      
      else if(density == F){ # For average density simulation, use the average density in site E, TSF>2
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteE" & mean.densities$TSF == "three"]
       pop.mat = dewy.buildmat.TSF3.siteE(dens1)
      }
     }
     
     else if(random.mat == 2){
      
      if(density == T){ # For density-dependent simulation, use the updated density
       
       pop.mat = dewy.buildmat.TSF3.siteF(dens1)
      }
      
      else if(density == F){ # For average density simulation, use the average density in site E, TSF>2
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteF" & mean.densities$TSF == "three"]
       pop.mat = dewy.buildmat.TSF3.siteF(dens1)
      }
     }
    }
    
    # Compute new population vector from the matrices
    n1 = pop.mat %*% n1
    
    # Capping the seedbank to 4000 to avoid population density explosion (twice the maximum observed amount of seeds in the soil)
    if(n1[1] > 4000) n1[1] = 4000
    
    # Capping the number of seedlings to 1700 to avoid population density explosion (twice the maximum observed amount of seedlings in the population)
    if(n1[2] > 1700) n1[2] = 1700
    
    # Compute above-ground density (without seedbank)
    if(density == T){
     
     if(random.mat == 1){
      
      dens1 = sum(n1[-1]) / nbsquares.siteE.TSF3
     }
     
     else if(random.mat == 2){
      
      dens1 = sum(n1[-1]) / nbsquares.siteF.TSF3
     }
    }
   }
   
   ###########################   
   
   ## TSF>3 - Stochastic case ----
   
   # For TSF>3, there are two possible perturbation scenarios: Either we perturb all the years in state 5 until a new fire occurs, or we only perturb a given proportion of those years.
   
   if(state == 5){
    
    random.mat = sample(c(1, 2), 1) # Picking a matrix between unperturbed sites E and F
    
    # Sample random year
    year = sample(unique(data.grazing$time)[1:8], 1)
    
    # Build MPM for state 5 (TSF>3)
    
    # Case (2) - Perturbed periodicity case
    
    if(any(perturbed.state == 5)){
     
      if(density == T){ # For density-dependent simulation, use the updated density
       
        pop.mat = dewy.buildmat.stochastic.siteC(dens1, year)
      }
      
      else if(density == F){ # For average density simulation, use the average density in site C, TSF>3
       
        dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "three"]
        pop.mat = dewy.buildmat.stochastic.siteC(dens1, year)
      }
    } 
    
    
    # Case (1) - Unperturbed periodicity case
    
    else{
     
     if(density == T){ # For density-dependent simulation, use the updated density
      
      pop.mat = dewy.buildmat.stochastic.sitesEF(dens1, year)
     }
     
     else if(density == F){ # For average density simulation, use the average density in site E or F, TSF>2
      
      if(random.mat == 1){
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteE" & mean.densities$TSF == "three"]
       pop.mat = dewy.buildmat.stochastic.sitesEF(dens1, year)
      }
      
      else if(random.mat == 2){
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteF" & mean.densities$TSF == "three"]
       pop.mat = dewy.buildmat.stochastic.sitesEF(dens1, year)
      }
     }
    }  
    
    # Compute new population vector from the matrices
    n1 = pop.mat %*% n1

    # Capping the seedbank to 4000 to avoid population density explosion (twice the maximum observed amount of seeds in the soil)
    if(n1[1] > 4000) n1[1] = 4000
    
    # Capping the number of seedlings to 1700 to avoid population density explosion (twice the maximum observed amount of seedlings in the population)
    if(n1[2] > 1700) n1[2] = 1700
    
    # Compute above-ground density (without seedbank)
    if(density == T){
     
     if(random.mat == 1){
      
      dens1 = sum(n1[-1]) / nbsquares.siteE.TSF3
     }
     
     else if(random.mat == 2){
      
      dens1 = sum(n1[-1]) / nbsquares.siteF.TSF3
     }
    }
   }
   
   # Record population vectors
   pop.vec.sim[[j]] = n1
   
   # Record year above-ground density
   dens[i, j] = dens1
   
   # Compute and record year lambda
   log.lambda[i, j] = log(sum(n1) / sum(pop.vec1))
   
   # Extinction case
   if(sum(n1[-1]) < 5 & n1[1] < 50){
    
    ext[i] = 1
    ext.year[i] = j
    break
   }
  }
  
  pop.vec[[i]] = pop.vec.sim 
  
 }
 return(list(log.lambda, dens, ext, ext.year, pop.vec, post.fire.habitat))
}




###########################################################################
#
# 3. Performing the simulations ----
#
###########################################################################

# We look at the effect of human-perturbed vital rates in different post-fire states (TSF0, TSF1, TSF2, TSF3, and TSF>3). 
# We look at the effect of human-perturbed vital rates in different proportions (0.1, 0.3, 0.5, 0.7, and 0.9) of the years spent in TSF>3.

# We fix the fire return to p = 0.033, i.e. once every 30 years.

# We perform 500 simulations of 100-year projections.


# Initial population vectors for each site
n0 = c(1000, 0, 0, 0, 0)


## 3.1. Tests ----
# -----------

## 3.1.1. Control test - Stochastic fires ----
# --------------------------------

test.sim0 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30)

# Densities
matplot(t(test.sim0[[2]]), type = "l", col = rainbow(nrow(test.sim0[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim0[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim0.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, density = F)


## 3.1.2. State 5 perturbed test - Stochastic fires ----
# -------------------------------------

test.sim1 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = 5)

# Densities
matplot(t(test.sim1[[2]]), type = "l", col = rainbow(nrow(test.sim1[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim1[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim1.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = 5, density = F)


## 3.1.3. States 4+5 perturbed test - Stochastic fires ----
# ----------------------------------------

test.sim2 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(4, 5))

# Densities
matplot(t(test.sim2[[2]]), type = "l", col = rainbow(nrow(test.sim2[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim2[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim2.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(4, 5), density = F)


## 3.1.4. States 3+4+5 perturbed test - Stochastic fires ----
# ------------------------------------------

test.sim3 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(3, 4, 5))

# Densities
matplot(t(test.sim3[[2]]), type = "l", col = rainbow(nrow(test.sim3[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim3[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim3.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(3, 4, 5), density = F)


## 3.1.5. States 2+3+4+5 perturbed test - Stochastic fires ----
# --------------------------------------------

test.sim4 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(2, 3, 4, 5))

# Densities
matplot(t(test.sim4[[2]]), type = "l", col = rainbow(nrow(test.sim4[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim4[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim4.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(2, 3, 4, 5), density = F)


## 3.1.6. States 1+2+3+4+5 perturbed test - Stochastic fires ----
# ----------------------------------------------

test.sim5 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(1, 2, 3, 4, 5))

# Densities
matplot(t(test.sim5[[2]]), type = "l", col = rainbow(nrow(test.sim5[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim5[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim5.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(1, 2, 3, 4, 5), density = F)



## 3.1.1. Control test - Periodic fires ----
# -----------------------------------------

test.sim0.periodicFire = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15)

# Densities
matplot(t(test.sim0.periodicFire[[2]]), type = "l", col = rainbow(nrow(test.sim0.periodicFire[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim0.periodicFire[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim0.periodicFire.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, density = F)


## 3.1.2. State 5 perturbed test - Periodic fires ----
# ---------------------------------------------------

test.sim1.periodicFire = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = 5)

# Densities
matplot(t(test.sim1.periodicFire[[2]]), type = "l", col = rainbow(nrow(test.sim1.periodicFire[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim1.periodicFire[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim1.periodicFire.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = 5, density = F)


## 3.1.3. States 4+5 perturbed test - Periodic fires ----
# ------------------------------------------------------

test.sim2.periodicFire = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(4, 5))

# Densities
matplot(t(test.sim2.periodicFire[[2]]), type = "l", col = rainbow(nrow(test.sim2.periodicFire[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim2.periodicFire[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim2.periodicFire.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(4, 5), density = F)


## 3.1.4. States 3+4+5 perturbed test - Periodic fires ----
# --------------------------------------------------------

test.sim3.periodicFire = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(3, 4, 5))

# Densities
matplot(t(test.sim3.periodicFire[[2]]), type = "l", col = rainbow(nrow(test.sim3.periodicFire[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim3.periodicFire[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim3.periodicFire.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(3, 4, 5), density = F)


## 3.1.5. States 2+3+4+5 perturbed test - Periodic fires ----
# ----------------------------------------------------------

test.sim4.periodicFire = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(2, 3, 4, 5))

# Densities
matplot(t(test.sim4.periodicFire[[2]]), type = "l", col = rainbow(nrow(test.sim4.periodicFire[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim4.periodicFire[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim4.periodicFire.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(2, 3, 4, 5), density = F)


## 3.1.6. States 1+2+3+4+5 perturbed test - Periodic fires ----
# ------------------------------------------------------------

test.sim5.periodicFire = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(1, 2, 3, 4, 5))

# Densities
matplot(t(test.sim5.periodicFire[[2]]), type = "l", col = rainbow(nrow(test.sim5.periodicFire[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim5.periodicFire[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim5.periodicFire.avgdens = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(1, 2, 3, 4, 5), density = F)



## 3.2. With density dependence and stochastic fire disturbance (p = 1/30) ----
# ------------------------------------------------------------------------

## 3.2.1. Unperturbed periodicity - Control case
# ----------------------------------------------

print("CONTROL (UNPERTURBED PERIODICITY) - p = 1/30")
stochasticFire30.control.sim = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30)

# Simulation-wise stochastic lambdas
stochasticFire30.control.sim.lambda = apply(stochasticFire30.control.sim[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.control.stoch.lambda = mean(stochasticFire30.control.sim.lambda)

# Variance in log lambda across years
stochasticFire30.control.var.lambda = apply(stochasticFire30.control.sim[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.control.ext.prob = mean(stochasticFire30.control.sim[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.control.sim[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.control.sim[[1]])), main = "Simulations yearly log lambda - Control scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.control.sim[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.control.sim[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.control.sim[[2]])), main = "Density over seasonal time steps - Control scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.control.sim[[2]], 2, mean, na.rm = T), lwd = 2) # Average density accross simulations for each time step


## 3.2.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("PERTURBED STATE 5 - p = 1/30")
stochasticFire30.perturb.state5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = 5)

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state5.sim.lambda = apply(stochasticFire30.perturb.state5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state5.stoch.lambda = mean(stochasticFire30.perturb.state5.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state5.var.lambda = apply(stochasticFire30.perturb.state5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state5.ext.prob = mean(stochasticFire30.perturb.state5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state5[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state5[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state5[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state5[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.3. Perturbed state = 4 + 5 (TSF3 and TSF>3) ----
# ------------------------------------------------

print("PERTURBED STATES 4 + 5 - p = 1/30")
stochasticFire30.perturb.state4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(4, 5))

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state4_5.sim.lambda = apply(stochasticFire30.perturb.state4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state4_5.stoch.lambda = mean(stochasticFire30.perturb.state4_5.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state4_5.var.lambda = apply(stochasticFire30.perturb.state4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state4_5.ext.prob = mean(stochasticFire30.perturb.state4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state4_5[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state4_5[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state4_5[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state4_5[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("PERTURBED STATES 3 + 4 + 5 - p = 1/30")
stochasticFire30.perturb.state3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(3, 4, 5))

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state3_4_5.sim.lambda = apply(stochasticFire30.perturb.state3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state3_4_5.stoch.lambda = mean(stochasticFire30.perturb.state3_4_5.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state3_4_5.var.lambda = apply(stochasticFire30.perturb.state3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state3_4_5.ext.prob = mean(stochasticFire30.perturb.state3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state3_4_5[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state3_4_5[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state3_4_5[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state3_4_5[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("PERTURBED STATES 2 + 3 + 4 + 5 - p = 1/30")
stochasticFire30.perturb.state2_3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(2, 3, 4, 5))

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state2_3_4_5.sim.lambda = apply(stochasticFire30.perturb.state2_3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state2_3_4_5.stoch.lambda = mean(stochasticFire30.perturb.state2_3_4_5.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state2_3_4_5.var.lambda = apply(stochasticFire30.perturb.state2_3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state2_3_4_5.ext.prob = mean(stochasticFire30.perturb.state2_3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state2_3_4_5[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state2_3_4_5[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state2_3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state2_3_4_5[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state2_3_4_5[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state2_3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.6. Perturbed states =  all ----
# -------------------------------

print("PERTURBED STATES ALL - p = 1/30")
stochasticFire30.perturb.stateAll = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(1, 2, 3, 4, 5))

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.stateAll.sim.lambda = apply(stochasticFire30.perturb.stateAll[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.stateAll.stoch.lambda = mean(stochasticFire30.perturb.stateAll.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.stateAll.var.lambda = apply(stochasticFire30.perturb.stateAll[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.stateAll.ext.prob = mean(stochasticFire30.perturb.stateAll[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.stateAll[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.stateAll[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.stateAll[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.stateAll[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.stateAll[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (stochastic fire, p = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.stateAll[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3. With density dependence and periodic fire disturbance (f = 1/15) ----
# ----------------------------------------------------------------------

## 3.3.1. Unperturbed periodicity - Control case
# ----------------------------------------------

print("CONTROL (UNPERTURBED PERIODICITY) - freq = 15 years")
periodicFire15.control.sim = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15)

# Simulation-wise stochastic lambdas
periodicFire15.control.sim.lambda = apply(periodicFire15.control.sim[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.control.stoch.lambda = mean(periodicFire15.control.sim.lambda)

# Variance in log lambda across years
periodicFire15.control.var.lambda = apply(periodicFire15.control.sim[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.control.ext.prob = mean(periodicFire15.control.sim[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.control.sim[[1]]), type = "l", col = rainbow(nrow(periodicFire15.control.sim[[1]])), main = "Simulations yearly log lambda - Control scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.control.sim[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.control.sim[[2]]), type = "l", col = rainbow(nrow(periodicFire15.control.sim[[2]])), main = "Density over seasonal time steps - Control scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.control.sim[[2]], 2, mean, na.rm = T), lwd = 2) # Average density accross simulations for each time step


## 3.3.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("PERTURBED STATE 5 - freq = 15 years")
periodicFire15.perturb.state5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = 5)

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state5.sim.lambda = apply(periodicFire15.perturb.state5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state5.stoch.lambda = mean(periodicFire15.perturb.state5.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state5.var.lambda = apply(periodicFire15.perturb.state5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state5.ext.prob = mean(periodicFire15.perturb.state5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state5[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state5[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state5[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state5[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.3. Perturbed state = 4 + 5 (TSF3 and TSF>3) ----
# ------------------------------------------------

print("PERTURBED STATES 4 + 5 - freq = 15 years")
periodicFire15.perturb.state4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(4, 5))

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state4_5.sim.lambda = apply(periodicFire15.perturb.state4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state4_5.stoch.lambda = mean(periodicFire15.perturb.state4_5.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state4_5.var.lambda = apply(periodicFire15.perturb.state4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state4_5.ext.prob = mean(periodicFire15.perturb.state4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state4_5[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state4_5[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state4_5[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state4_5[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("PERTURBED STATES 3 + 4 + 5 - freq = 15 years")
periodicFire15.perturb.state3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(3, 4, 5))

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state3_4_5.sim.lambda = apply(periodicFire15.perturb.state3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state3_4_5.stoch.lambda = mean(periodicFire15.perturb.state3_4_5.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state3_4_5.var.lambda = apply(periodicFire15.perturb.state3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state3_4_5.ext.prob = mean(periodicFire15.perturb.state3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state3_4_5[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state3_4_5[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state3_4_5[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state3_4_5[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("PERTURBED STATES 2 + 3 + 4 + 5 - freq = 15 years")
periodicFire15.perturb.state2_3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(2, 3, 4, 5))

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state2_3_4_5.sim.lambda = apply(periodicFire15.perturb.state2_3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state2_3_4_5.stoch.lambda = mean(periodicFire15.perturb.state2_3_4_5.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state2_3_4_5.var.lambda = apply(periodicFire15.perturb.state2_3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state2_3_4_5.ext.prob = mean(periodicFire15.perturb.state2_3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state2_3_4_5[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state2_3_4_5[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state2_3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state2_3_4_5[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state2_3_4_5[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state2_3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.6. Perturbed states =  all ----
# -------------------------------

print("PERTURBED STATES ALL - freq = 15 years")
periodicFire15.perturb.stateAll = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(1, 2, 3, 4, 5))

# Simulation-wise stochastic lambdas
periodicFire15.perturb.stateAll.sim.lambda = apply(periodicFire15.perturb.stateAll[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.stateAll.stoch.lambda = mean(periodicFire15.perturb.stateAll.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.stateAll.var.lambda = apply(periodicFire15.perturb.stateAll[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.stateAll.ext.prob = mean(periodicFire15.perturb.stateAll[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.stateAll[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.stateAll[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.stateAll[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.stateAll[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.stateAll[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.stateAll[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.4. At average density and stochastic fire disturbance (p = 1/30) ----
# -------------------------------------------------------------------

## 3.4.1. Unperturbed periodicity - Control case ----
# ----------------------------------------------

print("AVG DENSITY - CONTROL (UNPERTURBED PERIODICITY) - p = 1/30")
stochasticFire30.control.sim.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, density = F)

# Simulation-wise stochastic lambdas
stochasticFire30.control.avgdens.sim.lambda = apply(stochasticFire30.control.sim.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.control.avgdens.stoch.lambda = mean(stochasticFire30.control.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire30.control.avgdens.var.lambda = apply(stochasticFire30.control.sim.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.control.avgdens.ext.prob = mean(stochasticFire30.control.sim.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.control.sim.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.control.sim.avgdens[[1]])), main = "Simulations yearly log lambda - Control scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.control.sim.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.control.sim.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.control.sim.avgdens[[2]])), main = "Density over seasonal time steps - Control scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.control.sim.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.4.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("AVG DENSITY - PERTURBED STATE 5 - p = 1/30")
stochasticFire30.perturb.state5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = 5, density = F)

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state5.avgdens.sim.lambda = apply(stochasticFire30.perturb.state5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state5.avgdens.stoch.lambda = mean(stochasticFire30.perturb.state5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state5.avgdens.var.lambda = apply(stochasticFire30.perturb.state5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state5.avgdens.ext.prob = mean(stochasticFire30.perturb.state5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state5.avgdens[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state5.avgdens[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.4.3. Perturbed states = 4 + 5 (TSF3 and TSF>3) ----
# -------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 4 + 5 - p = 1/30")
stochasticFire30.perturb.state4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state4_5.avgdens.sim.lambda = apply(stochasticFire30.perturb.state4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state4_5.avgdens.stoch.lambda = mean(stochasticFire30.perturb.state4_5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state4_5.avgdens.var.lambda = apply(stochasticFire30.perturb.state4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state4_5.avgdens.ext.prob = mean(stochasticFire30.perturb.state4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.4.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 3 + 4 + 5 - p = 1/30")
stochasticFire30.perturb.state3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state3_4_5.avgdens.sim.lambda = apply(stochasticFire30.perturb.state3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state3_4_5.avgdens.stoch.lambda = mean(stochasticFire30.perturb.state3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state3_4_5.avgdens.var.lambda = apply(stochasticFire30.perturb.state3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state3_4_5.avgdens.ext.prob = mean(stochasticFire30.perturb.state3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.4.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATE 2 + 3 + 4 + 5 - p = 1/30")
stochasticFire30.perturb.state2_3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.state2_3_4_5.avgdens.sim.lambda = apply(stochasticFire30.perturb.state2_3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.state2_3_4_5.avgdens.stoch.lambda = mean(stochasticFire30.perturb.state2_3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.state2_3_4_5.avgdens.var.lambda = apply(stochasticFire30.perturb.state2_3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.state2_3_4_5.avgdens.ext.prob = mean(stochasticFire30.perturb.state2_3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.state2_3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state2_3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.state2_3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.state2_3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.state2_3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.state2_3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.4.6. Perturbed states = all ----
# ------------------------------

print("AVG DENSITY - PERTURBED STATES ALL - p = 1/30")
stochasticFire30.perturb.stateAll.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/30, perturbed.state = c(1, 2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire30.perturb.stateAll.avgdens.sim.lambda = apply(stochasticFire30.perturb.stateAll.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire30.perturb.stateAll.avgdens.stoch.lambda = mean(stochasticFire30.perturb.stateAll.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire30.perturb.stateAll.avgdens.var.lambda = apply(stochasticFire30.perturb.stateAll.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire30.perturb.stateAll.avgdens.ext.prob = mean(stochasticFire30.perturb.stateAll.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire30.perturb.stateAll.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.stateAll.avgdens[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire30.perturb.stateAll.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire30.perturb.stateAll.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire30.perturb.stateAll.avgdens[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (stochastic fire, p = 30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire30.perturb.stateAll.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step



## 3.5. At average density and periodic fire disturbance (f = 15) ----
# ---------------------------------------------------------------

## 3.5.1. Unperturbed periodicity - Control case ----
# ----------------------------------------------

print("AVG DENSITY - CONTROL (UNPERTURBED PERIODICITY) - freq = 15 years")
periodicFire15.control.sim.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, density = F)

# Simulation-wise stochastic lambdas
periodicFire15.control.avgdens.sim.lambda = apply(periodicFire15.control.sim.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.control.avgdens.stoch.lambda = mean(periodicFire15.control.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire15.control.avgdens.var.lambda = apply(periodicFire15.control.sim.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.control.avgdens.ext.prob = mean(periodicFire15.control.sim.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.control.sim.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire15.control.sim.avgdens[[1]])), main = "Simulations yearly log lambda - Control scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.control.sim.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.control.sim.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire15.control.sim.avgdens[[2]])), main = "Density over seasonal time steps - Control scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.control.sim.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.5.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("AVG DENSITY - PERTURBED STATE 5 - freq = 15 years")
periodicFire15.perturb.state5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = 5, density = F)

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state5.avgdens.sim.lambda = apply(periodicFire15.perturb.state5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state5.avgdens.stoch.lambda = mean(periodicFire15.perturb.state5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state5.avgdens.var.lambda = apply(periodicFire15.perturb.state5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state5.avgdens.ext.prob = mean(periodicFire15.perturb.state5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state5.avgdens[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state5.avgdens[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.5.3. Perturbed states = 4 + 5 (TSF3 and TSF>3) ----
# -------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 4 + 5 - freq = 15 years")
periodicFire15.perturb.state4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state4_5.avgdens.sim.lambda = apply(periodicFire15.perturb.state4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state4_5.avgdens.stoch.lambda = mean(periodicFire15.perturb.state4_5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state4_5.avgdens.var.lambda = apply(periodicFire15.perturb.state4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state4_5.avgdens.ext.prob = mean(periodicFire15.perturb.state4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.5.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 3 + 4 + 5 - freq = 15 years")
periodicFire15.perturb.state3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state3_4_5.avgdens.sim.lambda = apply(periodicFire15.perturb.state3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state3_4_5.avgdens.stoch.lambda = mean(periodicFire15.perturb.state3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state3_4_5.avgdens.var.lambda = apply(periodicFire15.perturb.state3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state3_4_5.avgdens.ext.prob = mean(periodicFire15.perturb.state3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.5.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATE 2 + 3 + 4 + 5 - freq = 15 years")
periodicFire15.perturb.state2_3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire15.perturb.state2_3_4_5.avgdens.sim.lambda = apply(periodicFire15.perturb.state2_3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.state2_3_4_5.avgdens.stoch.lambda = mean(periodicFire15.perturb.state2_3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.state2_3_4_5.avgdens.var.lambda = apply(periodicFire15.perturb.state2_3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.state2_3_4_5.avgdens.ext.prob = mean(periodicFire15.perturb.state2_3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.state2_3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state2_3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.state2_3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.state2_3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.state2_3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.state2_3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.5.6. Perturbed states = all ----
# ------------------------------

print("AVG DENSITY - PERTURBED STATES ALL - freq = 15 years")
periodicFire15.perturb.stateAll.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 15, perturbed.state = c(1, 2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire15.perturb.stateAll.avgdens.sim.lambda = apply(periodicFire15.perturb.stateAll.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire15.perturb.stateAll.avgdens.stoch.lambda = mean(periodicFire15.perturb.stateAll.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire15.perturb.stateAll.avgdens.var.lambda = apply(periodicFire15.perturb.stateAll.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire15.perturb.stateAll.avgdens.ext.prob = mean(periodicFire15.perturb.stateAll.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire15.perturb.stateAll.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.stateAll.avgdens[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire15.perturb.stateAll.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire15.perturb.stateAll.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire15.perturb.stateAll.avgdens[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (periodic fire, f = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire15.perturb.stateAll.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.6. With density dependence and stochastic fire disturbance (p = 15) ----
# ----------------------------------------------------------------------

## 3.6.1. Unperturbed periodicity - Control case
# ----------------------------------------------

print("CONTROL (UNPERTURBED PERIODICITY) - p = 1/15")
stochasticFire15.control.sim = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15)

# Simulation-wise stochastic lambdas
stochasticFire15.control.sim.lambda = apply(stochasticFire15.control.sim[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.control.stoch.lambda = mean(stochasticFire15.control.sim.lambda)

# Variance in log lambda across years
stochasticFire15.control.var.lambda = apply(stochasticFire15.control.sim[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.control.ext.prob = mean(stochasticFire15.control.sim[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.control.sim[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.control.sim[[1]])), main = "Simulations yearly log lambda - Control scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.control.sim[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.control.sim[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.control.sim[[2]])), main = "Density over seasonal time steps - Control scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.control.sim[[2]], 2, mean, na.rm = T), lwd = 2) # Average density accross simulations for each time step


## 3.6.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("PERTURBED STATE 5 - p = 1/15")
stochasticFire15.perturb.state5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = 5)

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state5.sim.lambda = apply(stochasticFire15.perturb.state5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state5.stoch.lambda = mean(stochasticFire15.perturb.state5.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state5.var.lambda = apply(stochasticFire15.perturb.state5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state5.ext.prob = mean(stochasticFire15.perturb.state5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state5[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state5[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state5[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state5[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.6.3. Perturbed state = 4 + 5 (TSF3 and TSF>3) ----
# ------------------------------------------------

print("PERTURBED STATES 4 + 5 - p = 1/15")
stochasticFire15.perturb.state4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(4, 5))

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state4_5.sim.lambda = apply(stochasticFire15.perturb.state4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state4_5.stoch.lambda = mean(stochasticFire15.perturb.state4_5.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state4_5.var.lambda = apply(stochasticFire15.perturb.state4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state4_5.ext.prob = mean(stochasticFire15.perturb.state4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state4_5[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state4_5[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state4_5[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state4_5[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.6.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("PERTURBED STATES 3 + 4 + 5 - p = 1/15")
stochasticFire15.perturb.state3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(3, 4, 5))

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state3_4_5.sim.lambda = apply(stochasticFire15.perturb.state3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state3_4_5.stoch.lambda = mean(stochasticFire15.perturb.state3_4_5.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state3_4_5.var.lambda = apply(stochasticFire15.perturb.state3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state3_4_5.ext.prob = mean(stochasticFire15.perturb.state3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state3_4_5[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state3_4_5[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state3_4_5[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state3_4_5[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.6.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("PERTURBED STATES 2 + 3 + 4 + 5 - p = 1/15")
stochasticFire15.perturb.state2_3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(2, 3, 4, 5))

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state2_3_4_5.sim.lambda = apply(stochasticFire15.perturb.state2_3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state2_3_4_5.stoch.lambda = mean(stochasticFire15.perturb.state2_3_4_5.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state2_3_4_5.var.lambda = apply(stochasticFire15.perturb.state2_3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state2_3_4_5.ext.prob = mean(stochasticFire15.perturb.state2_3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state2_3_4_5[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state2_3_4_5[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state2_3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state2_3_4_5[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state2_3_4_5[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state2_3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.6.6. Perturbed states =  all ----
# -------------------------------

print("PERTURBED STATES ALL - p = 1/15")
stochasticFire15.perturb.stateAll = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(1, 2, 3, 4, 5))

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.stateAll.sim.lambda = apply(stochasticFire15.perturb.stateAll[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.stateAll.stoch.lambda = mean(stochasticFire15.perturb.stateAll.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.stateAll.var.lambda = apply(stochasticFire15.perturb.stateAll[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.stateAll.ext.prob = mean(stochasticFire15.perturb.stateAll[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.stateAll[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.stateAll[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.stateAll[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.stateAll[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.stateAll[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.stateAll[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.7. At average density and stochastic fire disturbance (p = 1/15) ----
# -------------------------------------------------------------------

## 3.7.1. Unperturbed periodicity - Control case ----
# ----------------------------------------------

print("AVG DENSITY - CONTROL (UNPERTURBED PERIODICITY) - p = 1/15")
stochasticFire15.control.sim.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, density = F)

# Simulation-wise stochastic lambdas
stochasticFire15.control.avgdens.sim.lambda = apply(stochasticFire15.control.sim.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.control.avgdens.stoch.lambda = mean(stochasticFire15.control.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire15.control.avgdens.var.lambda = apply(stochasticFire15.control.sim.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.control.avgdens.ext.prob = mean(stochasticFire15.control.sim.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.control.sim.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.control.sim.avgdens[[1]])), main = "Simulations yearly log lambda - Control scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.control.sim.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.control.sim.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.control.sim.avgdens[[2]])), main = "Density over seasonal time steps - Control scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.control.sim.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.7.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("AVG DENSITY - PERTURBED STATE 5 - p = 1/15")
stochasticFire15.perturb.state5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = 5, density = F)

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state5.avgdens.sim.lambda = apply(stochasticFire15.perturb.state5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state5.avgdens.stoch.lambda = mean(stochasticFire15.perturb.state5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state5.avgdens.var.lambda = apply(stochasticFire15.perturb.state5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state5.avgdens.ext.prob = mean(stochasticFire15.perturb.state5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state5.avgdens[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state5.avgdens[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.7.3. Perturbed states = 4 + 5 (TSF3 and TSF>3) ----
# -------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 4 + 5 - p = 1/15")
stochasticFire15.perturb.state4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state4_5.avgdens.sim.lambda = apply(stochasticFire15.perturb.state4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state4_5.avgdens.stoch.lambda = mean(stochasticFire15.perturb.state4_5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state4_5.avgdens.var.lambda = apply(stochasticFire15.perturb.state4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state4_5.avgdens.ext.prob = mean(stochasticFire15.perturb.state4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.7.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 3 + 4 + 5 - p = 1/15")
stochasticFire15.perturb.state3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state3_4_5.avgdens.sim.lambda = apply(stochasticFire15.perturb.state3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state3_4_5.avgdens.stoch.lambda = mean(stochasticFire15.perturb.state3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state3_4_5.avgdens.var.lambda = apply(stochasticFire15.perturb.state3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state3_4_5.avgdens.ext.prob = mean(stochasticFire15.perturb.state3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.7.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATE 2 + 3 + 4 + 5 - p = 1/15")
stochasticFire15.perturb.state2_3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.state2_3_4_5.avgdens.sim.lambda = apply(stochasticFire15.perturb.state2_3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.state2_3_4_5.avgdens.stoch.lambda = mean(stochasticFire15.perturb.state2_3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.state2_3_4_5.avgdens.var.lambda = apply(stochasticFire15.perturb.state2_3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.state2_3_4_5.avgdens.ext.prob = mean(stochasticFire15.perturb.state2_3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.state2_3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state2_3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.state2_3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.state2_3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.state2_3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.state2_3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.7.6. Perturbed states = all ----
# ------------------------------

print("AVG DENSITY - PERTURBED STATES ALL - p = 1/15")
stochasticFire15.perturb.stateAll.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "stochastic", fire.prob = 1/15, perturbed.state = c(1, 2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
stochasticFire15.perturb.stateAll.avgdens.sim.lambda = apply(stochasticFire15.perturb.stateAll.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
stochasticFire15.perturb.stateAll.avgdens.stoch.lambda = mean(stochasticFire15.perturb.stateAll.avgdens.sim.lambda)

# Variance in log lambda across years
stochasticFire15.perturb.stateAll.avgdens.var.lambda = apply(stochasticFire15.perturb.stateAll.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
stochasticFire15.perturb.stateAll.avgdens.ext.prob = mean(stochasticFire15.perturb.stateAll.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(stochasticFire15.perturb.stateAll.avgdens[[1]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.stateAll.avgdens[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(stochasticFire15.perturb.stateAll.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(stochasticFire15.perturb.stateAll.avgdens[[2]]), type = "l", col = rainbow(nrow(stochasticFire15.perturb.stateAll.avgdens[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (stochastic fire, p = 1/15)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(stochasticFire15.perturb.stateAll.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.8. With density dependence and periodic fire disturbance (f = 1/30) ----
# ----------------------------------------------------------------------

## 3.8.1. Unperturbed periodicity - Control case
# ----------------------------------------------

print("CONTROL (UNPERTURBED PERIODICITY) - freq = 30 years")
periodicFire30.control.sim = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30)

# Simulation-wise stochastic lambdas
periodicFire30.control.sim.lambda = apply(periodicFire30.control.sim[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.control.stoch.lambda = mean(periodicFire30.control.sim.lambda)

# Variance in log lambda across years
periodicFire30.control.var.lambda = apply(periodicFire30.control.sim[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.control.ext.prob = mean(periodicFire30.control.sim[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.control.sim[[1]]), type = "l", col = rainbow(nrow(periodicFire30.control.sim[[1]])), main = "Simulations yearly log lambda - Control scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.control.sim[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.control.sim[[2]]), type = "l", col = rainbow(nrow(periodicFire30.control.sim[[2]])), main = "Density over seasonal time steps - Control scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.control.sim[[2]], 2, mean, na.rm = T), lwd = 2) # Average density accross simulations for each time step


## 3.8.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("PERTURBED STATE 5 - freq = 30 years")
periodicFire30.perturb.state5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = 5)

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state5.sim.lambda = apply(periodicFire30.perturb.state5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state5.stoch.lambda = mean(periodicFire30.perturb.state5.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state5.var.lambda = apply(periodicFire30.perturb.state5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state5.ext.prob = mean(periodicFire30.perturb.state5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state5[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state5[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state5[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state5[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.8.3. Perturbed state = 4 + 5 (TSF3 and TSF>3) ----
# ------------------------------------------------

print("PERTURBED STATES 4 + 5 - freq = 30 years")
periodicFire30.perturb.state4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(4, 5))

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state4_5.sim.lambda = apply(periodicFire30.perturb.state4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state4_5.stoch.lambda = mean(periodicFire30.perturb.state4_5.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state4_5.var.lambda = apply(periodicFire30.perturb.state4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state4_5.ext.prob = mean(periodicFire30.perturb.state4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state4_5[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state4_5[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state4_5[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state4_5[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.8.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("PERTURBED STATES 3 + 4 + 5 - freq = 30 years")
periodicFire30.perturb.state3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(3, 4, 5))

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state3_4_5.sim.lambda = apply(periodicFire30.perturb.state3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state3_4_5.stoch.lambda = mean(periodicFire30.perturb.state3_4_5.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state3_4_5.var.lambda = apply(periodicFire30.perturb.state3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state3_4_5.ext.prob = mean(periodicFire30.perturb.state3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state3_4_5[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state3_4_5[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state3_4_5[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state3_4_5[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.8.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("PERTURBED STATES 2 + 3 + 4 + 5 - freq = 30 years")
periodicFire30.perturb.state2_3_4_5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(2, 3, 4, 5))

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state2_3_4_5.sim.lambda = apply(periodicFire30.perturb.state2_3_4_5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state2_3_4_5.stoch.lambda = mean(periodicFire30.perturb.state2_3_4_5.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state2_3_4_5.var.lambda = apply(periodicFire30.perturb.state2_3_4_5[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state2_3_4_5.ext.prob = mean(periodicFire30.perturb.state2_3_4_5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state2_3_4_5[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state2_3_4_5[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state2_3_4_5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state2_3_4_5[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state2_3_4_5[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state2_3_4_5[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.8.6. Perturbed states =  all ----
# -------------------------------

print("PERTURBED ALL STATES - freq = 30 years")
periodicFire30.perturb.stateAll = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(1, 2, 3, 4, 5))

# Simulation-wise stochastic lambdas
periodicFire30.perturb.stateAll.sim.lambda = apply(periodicFire30.perturb.stateAll[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.stateAll.stoch.lambda = mean(periodicFire30.perturb.stateAll.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.stateAll.var.lambda = apply(periodicFire30.perturb.stateAll[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.stateAll.ext.prob = mean(periodicFire30.perturb.stateAll[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.stateAll[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.stateAll[[1]])), main = "Simulations yearly log lambda - States all perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.stateAll[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.stateAll[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.stateAll[[2]])), main = "Density over seasonal time steps - States all perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.stateAll[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.9. At average density and periodic fire disturbance (f = 1/30) ----
# -----------------------------------------------------------------

## 3.9.1. Unperturbed periodicity - Control case ----
# ----------------------------------------------

print("AVG DENSITY - CONTROL (UNPERTURBED PERIODICITY) - freq = 30 years")
periodicFire30.control.sim.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, density = F)

# Simulation-wise stochastic lambdas
periodicFire30.control.avgdens.sim.lambda = apply(periodicFire30.control.sim.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.control.avgdens.stoch.lambda = mean(periodicFire30.control.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire30.control.avgdens.var.lambda = apply(periodicFire30.control.sim.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.control.avgdens.ext.prob = mean(periodicFire30.control.sim.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.control.sim.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire30.control.sim.avgdens[[1]])), main = "Simulations yearly log lambda - Control scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.control.sim.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.control.sim.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire30.control.sim.avgdens[[2]])), main = "Density over seasonal time steps - Control scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.control.sim.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.9.2. Perturbed state = 5 (TSF>3) ----
# ----------------------------------

print("AVG DENSITY - PERTURBED STATE 5 - freq = 30 years")
periodicFire30.perturb.state5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = 5, density = F)

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state5.avgdens.sim.lambda = apply(periodicFire30.perturb.state5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state5.avgdens.stoch.lambda = mean(periodicFire30.perturb.state5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state5.avgdens.var.lambda = apply(periodicFire30.perturb.state5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state5.avgdens.ext.prob = mean(periodicFire30.perturb.state5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state5.avgdens[[1]])), main = "Simulations yearly log lambda - State 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state5.avgdens[[2]])), main = "Density over seasonal time steps - State 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.9.3. Perturbed states = 4 + 5 (TSF3 and TSF>3) ----
# -------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 4 + 5 - freq = 30 years")
periodicFire30.perturb.state4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state4_5.avgdens.sim.lambda = apply(periodicFire30.perturb.state4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state4_5.avgdens.stoch.lambda = mean(periodicFire30.perturb.state4_5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state4_5.avgdens.var.lambda = apply(periodicFire30.perturb.state4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state4_5.avgdens.ext.prob = mean(periodicFire30.perturb.state4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.9.4. Perturbed states = 3 + 4 + 5 (TSF2, TSF3, and TSF>3) ----
# ------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATES 3 + 4 + 5 - freq = 30 years")
periodicFire30.perturb.state3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state3_4_5.avgdens.sim.lambda = apply(periodicFire30.perturb.state3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state3_4_5.avgdens.stoch.lambda = mean(periodicFire30.perturb.state3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state3_4_5.avgdens.var.lambda = apply(periodicFire30.perturb.state3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state3_4_5.avgdens.ext.prob = mean(periodicFire30.perturb.state3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.9.5. Perturbed states = 2 + 3 + 4 + 5 (TSF1, TSF2, TSF3, and TSF>3) ----
# ----------------------------------------------------------------------

print("AVG DENSITY - PERTURBED STATE 2 + 3 + 4 + 5 - freq = 30 years")
periodicFire30.perturb.state2_3_4_5.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire30.perturb.state2_3_4_5.avgdens.sim.lambda = apply(periodicFire30.perturb.state2_3_4_5.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.state2_3_4_5.avgdens.stoch.lambda = mean(periodicFire30.perturb.state2_3_4_5.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.state2_3_4_5.avgdens.var.lambda = apply(periodicFire30.perturb.state2_3_4_5.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.state2_3_4_5.avgdens.ext.prob = mean(periodicFire30.perturb.state2_3_4_5.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.state2_3_4_5.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state2_3_4_5.avgdens[[1]])), main = "Simulations yearly log lambda - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.state2_3_4_5.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.state2_3_4_5.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.state2_3_4_5.avgdens[[2]])), main = "Density over seasonal time steps - States 2 + 3 + 4 + 5 perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.state2_3_4_5.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.9.6. Perturbed states = all ----
# ------------------------------

print("AVG DENSITY - PERTURBED STATES ALL - freq = 30 years")
periodicFire30.perturb.stateAll.avgdens = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, disturbance.type = "periodic", fire.freq = 30, perturbed.state = c(1, 2, 3, 4, 5), density = F)

# Simulation-wise stochastic lambdas
periodicFire30.perturb.stateAll.avgdens.sim.lambda = apply(periodicFire30.perturb.stateAll.avgdens[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
periodicFire30.perturb.stateAll.avgdens.stoch.lambda = mean(periodicFire30.perturb.stateAll.avgdens.sim.lambda)

# Variance in log lambda across years
periodicFire30.perturb.stateAll.avgdens.var.lambda = apply(periodicFire30.perturb.stateAll.avgdens[[1]], 1, var, na.rm = T)

# Mean extinction probability
periodicFire30.perturb.stateAll.avgdens.ext.prob = mean(periodicFire30.perturb.stateAll.avgdens[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(periodicFire30.perturb.stateAll.avgdens[[1]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.stateAll.avgdens[[1]])), main = "Simulations yearly log lambda - All states perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(periodicFire30.perturb.stateAll.avgdens[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(periodicFire30.perturb.stateAll.avgdens[[2]]), type = "l", col = rainbow(nrow(periodicFire30.perturb.stateAll.avgdens[[2]])), main = "Density over seasonal time steps - All states perturbed scenario (periodic fire, f = 1/30)", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(periodicFire30.perturb.stateAll.avgdens[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step





###########################################################################
#
# 4. Saving files ----
#
###########################################################################

save.image("DewyPine_NewProjections_Results.RData")
