##########################################################################################################
#
# RScript complementing the article Demographic consequences of changes in environmental periodicity (Conquet et al., under review at Ecology).
#
# This script contains the function needed to project the dewy pine population dynamics with or
# without perturbed vital rates in each post-fire habitat state (time since fire, TSF). 
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
              fire.prob = 0.033, 
              perturbed.state = 0, 
              density = T, 
              state5.perturbed.years.prop = NA){
 
 # (1) Control case = Simulating what happens in the normal conditions: We use the populations in the fire-disturbed sites E and F as natural populations. Each year, we randomly pick an MPM between site E and F.
 # (2) Perturbed periodicity case = Perturbing the periodicity in post-fire habitats in each post-fire state independently through inclusion of human-disturbance conditions, while keeping fire-return probability to 0.033 (one fire every 30 years). To do so, we use the given state matrix of the human-disturbed site C.
 
 # Preparing arrays, vectors, and lists to record lambdas, extinction, density, population vectors, and post-fire habitat states
 
 log.lambda = array(NA, c(nsimul, nyears))
 ext = seq(0, 0, length.out = nsimul)
 dens = array(NA, c(nsimul, nyears))
 pop.vec = list()
 
 post.fire.habitat = list()
 
 # Building the Markov chain to create the succession of post-fire environmental states
 p = fire.prob
 envS = matrix(0, 5, 5)
 
 envS[2, 1] = envS[3, 2] = envS[4, 3] = envS[5, 4] = 1
 envS[1, 5] = p
 envS[5, 5] = 1 - p
 
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
  states[1] = 1
  
  for(x in 2:nyears){
   states[x] = env.at.t = sample(ncol(envS), 1, pr = envS[, env.at.t])
  }
  
  post.fire.habitat[[i]] = states
  
  # Computing the number of years to be perturbed in state 5
  
  nb.perturbed.years = round(unlist(lapply(split(states[states != 1], cumsum(states == 1)[states != 1]), FUN = function(x) length(which(x == 5))), use.names = F) * state5.perturbed.years.prop, 0)
  
  
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
     if(density == T){ # For density-dependent simulation, use the updated abundance
      pop.mat = dewy.buildmat.TSF1.siteC(dens1)
     }
     
     else if(density == F){ # For average abundance simulation, use the average abundance in site C, TSF1
      
      dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "one"]
      pop.mat = dewy.buildmat.TSF1.siteC(dens1)
     }
    }
    
    
    # Case (1) - Unperturbed periodicity case 
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      
      if(density == T){ # For density-dependent simulation, use the updated abundance
       pop.mat = dewy.buildmat.TSF1.siteE(dens1)
      }
      
      else if(density == F){# For average abundance simulation, use the average abundance in site E, TSF1
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteE" & mean.densities$TSF == "one"]
       pop.mat = dewy.buildmat.TSF1.siteE(dens1)
      }
     }
     
     else if(random.mat == 2){
      
      if(density == T){ # For density-dependent simulation, use the updated abundance
       pop.mat = dewy.buildmat.TSF1.siteF(dens1)
      }
      
      else if(density == F){ # For average abundance simulation, use the average abundance in site F, TSF1
       
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
     if(density == T){ # For density-dependent simulation, use the updated abundance
 
      pop.mat = dewy.buildmat.TSF2.siteC(dens1)
     }
     
     else if(density == F){ # For average abundance simulation, use the average abundance in site C, TSF2
      
      dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "two"]
      pop.mat = dewy.buildmat.TSF2.siteC(dens1)
     }
    }
    
    
    # Case (1) - Unperturbed periodicity case 
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      
      if(density == T){ # For density-dependent simulation, use the updated abundance
       pop.mat = dewy.buildmat.TSF2.siteE(dens1)
      }
      
      else if(density == F){ # For average abundance simulation, use the average abundance in site E, TSF2
       
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
     if(density == T){ # For density-dependent simulation, use the updated abundance
      
      pop.mat = dewy.buildmat.TSF3.siteC(dens1)
     }
     
     else if(density == F){ # For average abundance simulation, use the average abundance in site C, TSF3
      
      dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "three"]
      pop.mat = dewy.buildmat.TSF3.siteC(dens1)
     }
    } 
    
    # Case (1) - Unperturbed periodicity case
    
    else{
     
     # Use an undisturbed MPM between the undisturbed sites E and F (fire-disturbed, natural populations)
     if(random.mat == 1){
      
      if(density == T){ # For density-dependent simulation, use the updated abundance
       pop.mat = dewy.buildmat.TSF3.siteE(dens1)
      }
      
      else if(density == F){ # For average abundance simulation, use the average abundance in site E, TSF>2
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteE" & mean.densities$TSF == "three"]
       pop.mat = dewy.buildmat.TSF3.siteE(dens1)
      }
     }
     
     else if(random.mat == 2){
      
      if(density == T){ # For density-dependent simulation, use the updated abundance
       
       pop.mat = dewy.buildmat.TSF3.siteF(dens1)
      }
      
      else if(density == F){ # For average abundance simulation, use the average abundance in site E, TSF>2
       
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
     
     if(is.na(state5.perturbed.years.prop)){
      
      if(density == T){ # For density-dependent simulation, use the updated abundance
       
       pop.mat = dewy.buildmat.stochastic.siteC(dens1, year)
      }
      
      else if(density == F){ # For average abundance simulation, use the average abundance in site C, TSF>3
       
       dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "three"]
       pop.mat = dewy.buildmat.stochastic.siteC(dens1, year)
      }
     }
     
     else{
      
      if(states[j - 1] != 5){ # The first time the population enters state 5 since the last fire.
       
       nb.years.in.state5 = 1 # Keeping track of the number of years spent in state 5 since the last fire, to be able to stop the perturbation when reaching the right year.
       
       nb.ptbd.yrs = nb.perturbed.years[nb.fires] # Getting the number of years to perturb. 
      }
      
      else{nb.years.in.state5 = nb.years.in.state5 + 1}
      
      if(nb.years.in.state5 <= nb.ptbd.yrs){ # If the year we are in is below the number of years to perturb, we use the MPM of the human-perturbed site C.
       
       # Use the MPM of the human-perturbed site C in state 5
       if(density == T){ # For density-dependent simulation, use the updated abundance
        
        pop.mat = dewy.buildmat.stochastic.siteC(dens1, year)
       }
       
       else if(density == F){ # For average abundance simulation, use the average abundance in site C, TSF>3
        
        dens1 = mean.densities$mean.density[mean.densities$site == "siteC" & mean.densities$TSF == "three"]
        pop.mat = dewy.buildmat.stochastic.siteC(dens1, year)
       }
      }
      
      else{ # If the year we are in is above the number of years to perturb, we use the undisturbed MPM for sites E and F.
       
       if(density == T){ # For density-dependent simulation, use the updated abundance
        
        pop.mat = dewy.buildmat.stochastic.sitesEF(dens1, year)
       }
       
       else if(density == F){ # For average abundance simulation, use the average abundance in site E or F, TSF>3
        
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
     }
    } 
    
    
    # Case (1) - Unperturbed periodicity case
    
    else{
     
     if(density == T){ # For density-dependent simulation, use the updated abundance
      
      pop.mat = dewy.buildmat.stochastic.sitesEF(dens1, year)
     }
     
     else if(density == F){ # For average abundance simulation, use the average abundance in site E or F, TSF>2
      
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
    break
   }
  }
  
  pop.vec[[i]] = pop.vec.sim 
  
 }
 return(list(log.lambda, dens, ext, pop.vec, post.fire.habitat))
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

## (0) Control test

test.sim0 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033)

# Densities
matplot(t(test.sim0[[2]]), type = "l", col = rainbow(nrow(test.sim0[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim0[[2]], 2, mean, na.rm = T), lwd = 2)
test.sim0.avgab = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, density = F)


## (1) State 1 perturbed test

test.sim1 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 1)

# Densities

matplot(t(test.sim1[[2]]), type = "l", col = rainbow(nrow(test.sim1[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim1[[2]], 2, mean, na.rm = T), lwd = 2)
test.sim1.avgab = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 1, density = F)


## (2) State 2 perturbed test

test.sim2 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 2)

# Densities

matplot(t(test.sim2[[2]]), type = "l", col = rainbow(nrow(test.sim2[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim2[[2]], 2, mean, na.rm = T), lwd = 2)
test.sim2.avgab = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 2, density = F)


## (3) State 3 perturbed test

test.sim3 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 3)

# Densities

matplot(t(test.sim3[[2]]), type = "l", col = rainbow(nrow(test.sim3[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim3[[2]], 2, mean, na.rm = T), lwd = 2)
test.sim3.avgab = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 3, density = F)


## (4) State 4 perturbed test

test.sim4 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 4)

# Densities

matplot(t(test.sim4[[2]]), type = "l", col = rainbow(nrow(test.sim4[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim4[[2]], 2, mean, na.rm = T), lwd = 2)
test.sim4.avgab = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 4, density = F)


# (5) State 5 perturbed test

test.sim5 = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 5)

# Densities

matplot(t(test.sim5[[2]]), type = "l", col = rainbow(nrow(test.sim5[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim5[[2]], 2, mean, na.rm = T), lwd = 2)
test.sim5.avgab = stoch.sim.droso(nsimul = 10, nyears = 10, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F)


# (6-10) State 5 perturbed proportion test

test.sim6 = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.1)

# Densities

matplot(t(test.sim6[[2]]), type = "l", col = rainbow(nrow(test.sim6[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim6[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim6.avgab = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F, state5.perturbed.years.prop = 0.1)

test.sim7 = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.3)

# Densities

matplot(t(test.sim7[[2]]), type = "l", col = rainbow(nrow(test.sim7[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim7[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim7.avgab = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F, state5.perturbed.years.prop = 0.3)

test.sim8 = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.5)

# Densities

matplot(t(test.sim8[[2]]), type = "l", col = rainbow(nrow(test.sim8[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim8[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim8.avgab = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F, state5.perturbed.years.prop = 0.5)

test.sim9 = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.7)

# Densities

matplot(t(test.sim9[[2]]), type = "l", col = rainbow(nrow(test.sim9[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim9[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim9.avgab = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F, state5.perturbed.years.prop = 0.7)

test.sim10 = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.9)

# Densities

matplot(t(test.sim10[[2]]), type = "l", col = rainbow(nrow(test.sim10[[2]])), main = "Density over yearly time steps", xlab = "Time steps (years)", ylab = "Density")
lines(apply(test.sim10[[2]], 2, mean, na.rm = T), lwd = 2)

test.sim10.avgab = stoch.sim.droso(nsimul = 10, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F, state5.perturbed.years.prop = 0.9)



## 3.2. With density dependence ----
# -----------------------------

## 3.2.1. Unperturbed periodicity - Control case
# ----------------------------------------------

print("CONTROL (UNPERTURBED PERIODICITY) - p = 0.033")
control.sim = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033)

# Simulation-wise stochastic lambdas
control.sim.lambda = apply(control.sim[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
control.stoch.lambda = mean(control.sim.lambda)

# Variance in log lambda across years
control.var.lambda = apply(control.sim[[1]], 1, var, na.rm = T)

# Mean extinction probability
control.ext.prob = mean(control.sim[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(control.sim[[1]]), type = "l", col = rainbow(nrow(control.sim[[1]])), main = "Simulations yearly log lambda - Control scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(control.sim[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(control.sim[[2]]), type = "l", col = rainbow(nrow(control.sim[[2]])), main = "Density over seasonal time steps - Control scenario", xlab = "Time steps (years)", ylab = "Density") # Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(control.sim[[2]], 2, mean), lwd = 2) # Average density accross simulations for each time step


## 3.2.2. Perturbed state = 1 (TSF0) ----
# ----------------------------------

print("PERTURBED STAGE 1 - p = 0.033")
perturb.state1 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 1)

# Simulation-wise stochastic lambdas
perturb.state1.sim.lambda = apply(perturb.state1[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state1.stoch.lambda = mean(perturb.state1.sim.lambda)

# Variance in log lambda across years
perturb.state1.var.lambda = apply(perturb.state1[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state1.ext.prob = mean(perturb.state1[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state1[[1]]), type = "l", col = rainbow(nrow(perturb.state1[[1]])), main = "Simulations yearly log lambda - perturb.state1 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state1[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state1[[2]]), type = "l", col = rainbow(nrow(perturb.state1[[2]])), main = "Density over seasonal time steps - perturb.state1 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state1[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.2.3. Perturbed state = 2 (TSF1) ----
# ----------------------------------

print("PERTURBED STAGE 2 - p = 0.033")
perturb.state2 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 2)

# Simulation-wise stochastic lambdas
perturb.state2.sim.lambda = apply(perturb.state2[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state2.stoch.lambda = mean(perturb.state2.sim.lambda)

# Variance in log lambda across years
perturb.state2.var.lambda = apply(perturb.state2[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state2.ext.prob = mean(perturb.state2[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state2[[1]]), type = "l", col = rainbow(nrow(perturb.state2[[1]])), main = "Simulations yearly log lambda - perturb.state2 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state2[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state2[[2]]), type = "l", col = rainbow(nrow(perturb.state2[[2]])), main = "Density over seasonal time steps - perturb.state2 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state2[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.2.4. Perturbed state = 3 (TSF2) ----
# ----------------------------------

print("PERTURBED STAGE 3 - p = 0.033")
perturb.state3 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 3)

# Simulation-wise stochastic lambdas
perturb.state3.sim.lambda = apply(perturb.state3[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state3.stoch.lambda = mean(perturb.state3.sim.lambda)

# Variance in log lambda across years
perturb.state3.var.lambda = apply(perturb.state3[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state3.ext.prob = mean(perturb.state3[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state3[[1]]), type = "l", col = rainbow(nrow(perturb.state3[[1]])), main = "Simulations yearly log lambda - perturb.state3 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state3[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state3[[2]]), type = "l", col = rainbow(nrow(perturb.state3[[2]])), main = "Density over seasonal time steps - perturb.state3 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state3[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.2.5. Perturbed state = 4 (TSF3) ----
# ----------------------------------

print("PERTURBED STAGE 4 - p = 0.033")
perturb.state4 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 4)

# Simulation-wise stochastic lambdas
perturb.state4.sim.lambda = apply(perturb.state4[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state4.stoch.lambda = mean(perturb.state4.sim.lambda)

# Variance in log lambda across years
perturb.state4.var.lambda = apply(perturb.state4[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state4.ext.prob = mean(perturb.state4[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state4[[1]]), type = "l", col = rainbow(nrow(perturb.state4[[1]])), main = "Simulations yearly log lambda - perturb.state4 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state4[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state4[[2]]), type = "l", col = rainbow(nrow(perturb.state4[[2]])), main = "Density over seasonal time steps - perturb.state4 scenario ", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state4[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.2.6. Perturbed state = 5 (TSF>3) ----
# -----------------------------------

print("PERTURBED STAGE 5 - p = 0.033")
perturb.state5 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5)

# Simulation-wise stochastic lambdas
perturb.state5.sim.lambda = apply(perturb.state5[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.stoch.lambda = mean(perturb.state5.sim.lambda)

# Variance in log lambda across years
perturb.state5.var.lambda = apply(perturb.state5[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.ext.prob = mean(perturb.state5[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5[[1]]), type = "l", col = rainbow(nrow(perturb.state5[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5[[2]]), type = "l", col = rainbow(nrow(perturb.state5[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.2.7. Proportion of the years perturbed in state 5 = 0.1 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.1")
perturb.state5.proportion.10 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.1)

# Simulation-wise stochastic lambdas
perturb.state5.proportion.10.sim.lambda = apply(perturb.state5.proportion.10[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.proportion.10.stoch.lambda = mean(perturb.state5.proportion.10.sim.lambda)

# Variance in log lambda across years
perturb.state5.proportion.10.var.lambda = apply(perturb.state5.proportion.10[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.proportion.10.ext.prob = mean(perturb.state5.proportion.10[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.proportion.10[[1]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.10[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.proportion.10[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.proportion.10[[2]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.10[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.proportion.10[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.2.8. Proportion of the years perturbed in state 5 = 0.3 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.3")
perturb.state5.proportion.30 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.3)

# Simulation-wise stochastic lambdas
perturb.state5.proportion.30.sim.lambda = apply(perturb.state5.proportion.30[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.proportion.30.stoch.lambda = mean(perturb.state5.proportion.30.sim.lambda)

# Variance in log lambda across years
perturb.state5.proportion.30.var.lambda = apply(perturb.state5.proportion.30[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.proportion.30.ext.prob = mean(perturb.state5.proportion.30[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.proportion.30[[1]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.30[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.proportion.30[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.proportion.30[[2]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.30[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.proportion.30[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.9. Proportion of the years perturbed in state 5 = 0.5 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.5")
perturb.state5.proportion.50 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.5)

# Simulation-wise stochastic lambdas
perturb.state5.proportion.50.sim.lambda = apply(perturb.state5.proportion.50[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.proportion.50.stoch.lambda = mean(perturb.state5.proportion.50.sim.lambda)

# Variance in log lambda across years
perturb.state5.proportion.50.var.lambda = apply(perturb.state5.proportion.50[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.proportion.50.ext.prob = mean(perturb.state5.proportion.50[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.proportion.50[[1]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.50[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.proportion.50[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.proportion.50[[2]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.50[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.proportion.50[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.10. Proportion of the years perturbed in state 5 = 0.7 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.7")
perturb.state5.proportion.70 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.7)

# Simulation-wise stochastic lambdas
perturb.state5.proportion.70.sim.lambda = apply(perturb.state5.proportion.70[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.proportion.70.stoch.lambda = mean(perturb.state5.proportion.70.sim.lambda)

# Variance in log lambda across years
perturb.state5.proportion.70.var.lambda = apply(perturb.state5.proportion.70[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.proportion.70.ext.prob = mean(perturb.state5.proportion.70[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.proportion.70[[1]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.70[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.proportion.70[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.proportion.70[[2]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.70[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.proportion.70[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.2.11. Proportion of the years perturbed in state 5 = 0.9 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.9")
perturb.state5.proportion.90 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.9)

# Simulation-wise stochastic lambdas
perturb.state5.proportion.90.sim.lambda = apply(perturb.state5.proportion.90[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.proportion.90.stoch.lambda = mean(perturb.state5.proportion.90.sim.lambda)

# Variance in log lambda across years
perturb.state5.proportion.90.var.lambda = apply(perturb.state5.proportion.90[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.proportion.90.ext.prob = mean(perturb.state5.proportion.90[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.proportion.90[[1]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.90[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.proportion.90[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.proportion.90[[2]]), type = "l", col = rainbow(nrow(perturb.state5.proportion.90[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.proportion.90[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step




## 3.3. At average density ----
# ------------------------

## 3.3.1. Unperturbed periodicity - Control case ----
# ----------------------------------------------

print("AVG ABUNDANCE - CONTROL (UNPERTURBED PERIODICITY) - p = 0.033")
control.sim.avgab = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, density = F)

# Simulation-wise stochastic lambdas
control.avgab.sim.lambda = apply(control.sim.avgab[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
control.avgab.stoch.lambda = mean(control.avgab.sim.lambda)

# Variance in log lambda across years
control.avgab.var.lambda = apply(control.sim.avgab[[1]], 1, var, na.rm = T)

# Mean extinction probability
control.avgab.ext.prob = mean(control.sim.avgab[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(control.sim.avgab[[1]]), type = "l", col = rainbow(nrow(control.sim.avgab[[1]])), main = "Simulations yearly log lambda - Control scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(control.sim.avgab[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(control.sim.avgab[[2]]), type = "l", col = rainbow(nrow(control.sim.avgab[[2]])), main = "Density over seasonal time steps - Control scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(control.sim.avgab[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.2. Perturbed state = 1 (TSF0) ----
# ----------------------------------

print("AVG ABUNDANCE - PERTURBED STAGE 1 - p = 0.033")
perturb.state1.avgab = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 1, density = F)

# Simulation-wise stochastic lambdas
perturb.state1.avgab.sim.lambda = apply(perturb.state1.avgab[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state1.avgab.stoch.lambda = mean(perturb.state1.avgab.sim.lambda)

# Variance in log lambda across years
perturb.state1.avgab.var.lambda = apply(perturb.state1.avgab[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state1.avgab.ext.prob = mean(perturb.state1.avgab[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state1.avgab[[1]]), type = "l", col = rainbow(nrow(perturb.state1.avgab[[1]])), main = "Simulations yearly log lambda - perturb.state1 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state1.avgab[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state1.avgab[[2]]), type = "l", col = rainbow(nrow(perturb.state1.avgab[[2]])), main = "Density over seasonal time steps - perturb.state1 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state1.avgab[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.3. Perturbed state = 2 (TSF1) ----
# ----------------------------------

print("AVG ABUNDANCE - PERTURBED STAGE 2 - p = 0.033")
perturb.state2.avgab = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 2, density = F)

# Simulation-wise stochastic lambdas
perturb.state2.avgab.sim.lambda = apply(perturb.state2.avgab[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state2.avgab.stoch.lambda = mean(perturb.state2.avgab.sim.lambda)

# Variance in log lambda across years
perturb.state2.avgab.var.lambda = apply(perturb.state2.avgab[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state2.avgab.ext.prob = mean(perturb.state2.avgab[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state2.avgab[[1]]), type = "l", col = rainbow(nrow(perturb.state2.avgab[[1]])), main = "Simulations yearly log lambda - perturb.state2 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state2.avgab[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state2.avgab[[2]]), type = "l", col = rainbow(nrow(perturb.state2.avgab[[2]])), main = "Density over seasonal time steps - perturb.state2 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state2.avgab[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.4. Perturbed state = 3 (TSF2) ----
# ----------------------------------

print("AVG ABUNDANCE - PERTURBED STAGE 3 - p = 0.033")
perturb.state3.avgab = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 3, density = F)

# Simulation-wise stochastic lambdas
perturb.state3.avgab.sim.lambda = apply(perturb.state3.avgab[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state3.avgab.stoch.lambda = mean(perturb.state3.avgab.sim.lambda)

# Variance in log lambda across years
perturb.state3.avgab.var.lambda = apply(perturb.state3.avgab[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state3.avgab.ext.prob = mean(perturb.state3.avgab[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state3.avgab[[1]]), type = "l", col = rainbow(nrow(perturb.state3.avgab[[1]])), main = "Simulations yearly log lambda - perturb.state3 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state3.avgab[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state3.avgab[[2]]), type = "l", col = rainbow(nrow(perturb.state3.avgab[[2]])), main = "Density over seasonal time steps - perturb.state3 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state3.avgab[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.5. Perturbed state = 4 (TSF3) ----
# ----------------------------------

print("AVG ABUNDANCE - PERTURBED STAGE 4 - p = 0.033")
perturb.state4.avgab = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 4, density = F)

# Simulation-wise stochastic lambdas
perturb.state4.avgab.sim.lambda = apply(perturb.state4.avgab[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state4.avgab.stoch.lambda = mean(perturb.state4.avgab.sim.lambda)

# Variance in log lambda across years
perturb.state4.avgab.var.lambda = apply(perturb.state4.avgab[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state4.avgab.ext.prob = mean(perturb.state4.avgab[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state4.avgab[[1]]), type = "l", col = rainbow(nrow(perturb.state4.avgab[[1]])), main = "Simulations yearly log lambda - perturb.state4 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state4.avgab[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state4.avgab[[2]]), type = "l", col = rainbow(nrow(perturb.state4.avgab[[2]])), main = "Density over seasonal time steps - perturb.state4 scenario ", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state4.avgab[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.6. Perturbed state = 5 (TSF>3) ----
# -----------------------------------

print("AVG ABUNDANCE - PERTURBED STAGE 5 - p = 0.033")
perturb.state5.avgab = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, density = F)

# Simulation-wise stochastic lambdas
perturb.state5.avgab.sim.lambda = apply(perturb.state5.avgab[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.avgab.stoch.lambda = mean(perturb.state5.avgab.sim.lambda)

# Variance in log lambda across years
perturb.state5.avgab.var.lambda = apply(perturb.state5.avgab[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.avgab.ext.prob = mean(perturb.state5.avgab[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.avgab[[1]]), type = "l", col = rainbow(nrow(perturb.state5.avgab[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.avgab[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.avgab[[2]]), type = "l", col = rainbow(nrow(perturb.state5.avgab[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.avgab[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.7. Proportion of the years perturbed in state 5 = 0.1 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.1")
perturb.state5.avgab.proportion.10 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.1, density = F)

# Simulation-wise stochastic lambdas
perturb.state5.avgab.proportion.10.sim.lambda = apply(perturb.state5.avgab.proportion.10[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.avgab.proportion.10.stoch.lambda = mean(perturb.state5.avgab.proportion.10.sim.lambda)

# Variance in log lambda across years
perturb.state5.avgab.proportion.10.var.lambda = apply(perturb.state5.avgab.proportion.10[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.avgab.proportion.10.ext.prob = mean(perturb.state5.avgab.proportion.10[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.avgab.proportion.10[[1]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.10[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.avgab.proportion.10[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.avgab.proportion.10[[2]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.10[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.avgab.proportion.10[[2]], 2, mean), lwd = 2) #Average density accross simulations for each time step


## 3.3.8. Proportion of the years perturbed in state 5 = 0.3 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.3")
perturb.state5.avgab.proportion.30 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.3, density = F)

# Simulation-wise stochastic lambdas
perturb.state5.avgab.proportion.30.sim.lambda = apply(perturb.state5.avgab.proportion.30[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.avgab.proportion.30.stoch.lambda = mean(perturb.state5.avgab.proportion.30.sim.lambda)

# Variance in log lambda across years
perturb.state5.avgab.proportion.30.var.lambda = apply(perturb.state5.avgab.proportion.30[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.avgab.proportion.30.ext.prob = mean(perturb.state5.avgab.proportion.30[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.avgab.proportion.30[[1]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.30[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.avgab.proportion.30[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.avgab.proportion.30[[2]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.30[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.avgab.proportion.30[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.9. Proportion of the years perturbed in state 5 = 0.5 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.5")
perturb.state5.avgab.proportion.50 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.5, density = F)

# Simulation-wise stochastic lambdas
perturb.state5.avgab.proportion.50.sim.lambda = apply(perturb.state5.avgab.proportion.50[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.avgab.proportion.50.stoch.lambda = mean(perturb.state5.avgab.proportion.50.sim.lambda)

# Variance in log lambda across years
perturb.state5.avgab.proportion.50.var.lambda = apply(perturb.state5.avgab.proportion.50[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.avgab.proportion.50.ext.prob = mean(perturb.state5.avgab.proportion.50[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.avgab.proportion.50[[1]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.50[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.avgab.proportion.50[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.avgab.proportion.50[[2]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.50[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.avgab.proportion.50[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.10. Proportion of the years perturbed in state 5 = 0.7 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.7")
perturb.state5.avgab.proportion.70 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.7, density = F)

# Simulation-wise stochastic lambdas
perturb.state5.avgab.proportion.70.sim.lambda = apply(perturb.state5.avgab.proportion.70[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.avgab.proportion.70.stoch.lambda = mean(perturb.state5.avgab.proportion.70.sim.lambda)

# Variance in log lambda across years
perturb.state5.avgab.proportion.70.var.lambda = apply(perturb.state5.avgab.proportion.70[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.avgab.proportion.70.ext.prob = mean(perturb.state5.avgab.proportion.70[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.avgab.proportion.70[[1]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.70[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.avgab.proportion.70[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.avgab.proportion.70[[2]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.70[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.avgab.proportion.70[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step


## 3.3.11. Proportion of the years perturbed in state 5 = 0.9 ----
# ----------------------------------------------------------

print("PERTURBED STAGE 5 - prop = 0.9")
perturb.state5.avgab.proportion.90 = stoch.sim.droso(nsimul = 500, nyears = 100, n0 = n0, fire.prob = 0.033, perturbed.state = 5, state5.perturbed.years.prop = 0.9, density = F)

# Simulation-wise stochastic lambdas
perturb.state5.avgab.proportion.90.sim.lambda = apply(perturb.state5.avgab.proportion.90[[1]], 1, mean, na.rm = T)

# Mean stochastic lambda
perturb.state5.avgab.proportion.90.stoch.lambda = mean(perturb.state5.avgab.proportion.90.sim.lambda)

# Variance in log lambda across years
perturb.state5.avgab.proportion.90.var.lambda = apply(perturb.state5.avgab.proportion.90[[1]], 1, var, na.rm = T)

# Mean extinction probability
perturb.state5.avgab.proportion.90.ext.prob = mean(perturb.state5.avgab.proportion.90[[3]])

# Yearly log lambdas for each simulation, and mean
matplot(t(perturb.state5.avgab.proportion.90[[1]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.90[[1]])), main = "Simulations yearly log lambda - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Log lambda")
lines(apply(perturb.state5.avgab.proportion.90[[1]], 2, mean, na.rm = T), lwd = 2)

# Density plots
matplot(t(perturb.state5.avgab.proportion.90[[2]]), type = "l", col = rainbow(nrow(perturb.state5.avgab.proportion.90[[2]])), main = "Density over seasonal time steps - perturb.state5 scenario", xlab = "Time steps (years)", ylab = "Density") #Matplot plots columns of a matrix so we need to use t() for it to plot the rows
lines(apply(perturb.state5.avgab.proportion.90[[2]], 2, mean, na.rm = T), lwd = 2) #Average density accross simulations for each time step




###########################################################################
#
# 4. Saving files ----
#
###########################################################################

save.image("DewyPineProjectionsResults.RData")

