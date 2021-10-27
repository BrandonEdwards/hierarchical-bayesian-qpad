####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 1-Exp1.R
# Created October 2021
# Last Updated October 2021

####### Import Libraries and External Files #######

library(dirmult)
library(detect)
library(doParallel)
library(foreach)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# This is the (fairly) generic function to simulate point counts
source("src/functions/sim-pc.R")

####### Set Constants #############################

set.seed(seed = 6846,
         kind = "Mersenne-Twister",
         normal.kind = "Inversion")

n_cores <- 4
n_obs <- 300
n_sim <- 4

phi <- 0.3
tau <- 200
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

sim_data <- vector(mode = "list", length = n_sim)

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
                          poisson_lambda = 20)
}

####### Maximum Likelihood ########################

cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)

foreach (i = 1:n_sim, .packages = 'detect') %dopar%
{
  x <- sim_data[[i]]
  
  removal <- cmulti(as.matrix(x$rem[, 3:ncol(x$rem)]) | as.matrix(x$time_design) ~ 1, type = "rem")
  distance <- cmulti(as.matrix(x$dis[, 3:ncol(x$dis)]) | as.matrix(x$dist_design) ~ 1, type = "dis")
  
  model_list <- list(removal, distance)
  save(model_list, file = paste0("/output/exp1/mle_", i, ".rda"))
}

stopCluster(cluster)

####### Bayesian Model (Stan) #####################

i <- 1
x <- sim_data[[i]]

# Removal Modelling
rem <- x$rem[, 3:(max(n_time_bins) + 2)]

#' Corresponds with "bands_per_sample" in removal.stan
time_bands_per_sample <- unname(apply(rem,
                                      1,
                                      function (x) sum(!is.na(x))))

#' Total abundance per sampling event.
#' I.e., this is the sum of Y_ij over j
#' Corresponds with "abund_per_sample" in removal.stan
total_abund_per_sample <- unname(apply(rem, 1, function(x) sum(x, na.rm = TRUE)))

# Design matrix, just NULL model for now
X <- matrix(data = 1, nrow = n_obs, ncol = 1)

#' Corresponds with "abund_per_band" in removal.stan
abundance_per_band <- rem
abundance_per_band[is.na(abundance_per_band)] <- 0

#' Corresponds with "max_time" in removal.stan
time_design <- x$time_design
time_design[is.na(time_design)] <- 0

stan_data_rem <- list(n_samples = n_obs,
                      n_covariates = ncol(X),
                      max_intervals = ncol(rem),
                      abund_per_band = abundance_per_band,
                      abund_per_sample = total_abund_per_sample,
                      bands_per_sample = time_bands_per_sample,
                      max_time = time_design,
                      X = X)

model <- stan_model(file = "stan/removal.stan")

stime = system.time(stan_fit <- 
                      sampling(model,
                               data = stan_data_rem,
                               verbose = TRUE,
                               chains = 3,
                               iter = 200,
                               warmup = 100,
                               cores = 3,
                               pars = c("gamma"),
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15)))
