####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 1-Exp1.R
# Created October 2021
# Last Updated October 2021

####### Import Libraries and External Files #######

library(dirmult)

# This is the (fairly) generic function to simulate point counts
source("src/functions/sim-pc.R")

####### Set Constants #############################

set.seed(seed = 6846,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")

phi <- 0.3
tau <- 100
n_obs <- 1000
n_sim <- 50

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
                            25, 50, 75, 100, 125, 150, Inf, NA, NA, NA, NA, 
                            10, 20, 30, 40, 50, 60, 70, 80, 90, 100, Inf),
                   nrow = 4, byrow = TRUE)

####### Simulate Datasets #########################

sim_data <- vector(mode = "list", length = n_sim)

start_time <- Sys.time()

for (s in 1:n_sim)
{
  sim_data[[s]] <- sim_pc(n_obs = n_obs,
                          phi = phi,
                          tau = tau,
                          n_protocols = n_protocols,
                          n_time_bins = n_time_bins,
                          n_dist_bins = n_dist_bins,
                          max_times = max_times,
                          max_dist = max_dist,
                          poisson_lambda = 10)
}

end_time <- Sys.time()
end_time - start_time
