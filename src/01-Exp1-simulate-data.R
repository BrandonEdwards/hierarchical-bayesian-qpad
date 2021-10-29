####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 01-Exp1-simulate-data.R
# Created October 2021
# Last Updated October 2021

####### Import Libraries and External Files #######

library(dirmult)
library(bSims)
library(detect)
library(doParallel)
library(foreach)

# This is the (fairly) generic function to simulate point counts
source("src/functions/sim-pc.R")

####### Set Constants #############################

set.seed(seed = 6846,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")

n_cores <- 30
n_obs <- c(75,500, 1000)
n_sim <- 400

phi <- 0.4
tau <- 100
den <- 10
n_protocols <- 4
n_time_bins <- c(3, 4, 10, 7)
n_dist_bins <- c(3, 4, 7, 11)
max_times <- matrix(data = c(3, 4, 5, NA, NA, NA, NA, NA, NA, NA, 
                             2, 4, 5, 6, NA, NA, NA, NA, NA, NA,
                             1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                             2, 3, 4, 5, 6, 7, 8, NA, NA, NA),
                    nrow = 4, byrow = TRUE)
max_dist <- matrix(data = c(50, 100, 400, NA, NA, NA, NA, NA, NA, NA, NA,
                            25, 50, 100, 150, NA, NA, NA, NA, NA, NA, NA, 
                            25, 50, 75, 100, 125, 150, 400, NA, NA, NA, NA, 
                            10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150),
                   nrow = 4, byrow = TRUE)

####### Simulate Datasets #########################

for (n in n_obs)
{
  for (s in 1:n_sim)
  {
    sim_data <- sim_pc(n_obs = n,
                       phi = phi,
                       tau = tau,
                       den = den,
                       n_protocols = n_protocols,
                       n_time_bins = n_time_bins,
                       n_dist_bins = n_dist_bins,
                       max_times = max_times,
                       max_dist = max_dist,
                       n_cores = n_cores)
    
    save(sim_data, file = paste0("output/exp1/sim_data/n",
                                 n,
                                 "-",
                                 s,
                                 ".rda"))
  }  
}
