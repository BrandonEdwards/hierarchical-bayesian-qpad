####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 05-Exp2-run-models.R
# Created March 2022
# Last Updated March 2022

####### Import Libraries and External Files #######

library(doParallel)
library(foreach)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


####### Set Constants #############################

n_obs <- c(75,500, 1000)
n_sim <- 1

rem_stan <- stan_model(file = "stan/removal_strong-prior.stan")
dis_stan <- stan_model(file = "stan/distance_strong-prior.stan")

for (n in n_obs)
{
  for (i in 1:n_sim)
  {
    # Using the same simulated data from Experiment 1
    file_name_sim <- paste0("output/exp1/sim_data/n",
                            n,
                            "-",
                            i,
                            ".rda")
    file_name_dis <- paste0("output/exp2/dis_mods/n",
                            n,
                            "-",
                            i,
                            "_dis_bayes.rda")
    file_name_rem <- paste0("output/exp2/rem_mods/n",
                            n,
                            "-",
                            i,
                            "_rem_bayes.rda")
    load(file_name_sim)
    
    ############ Removal Modelling loop ################
    rem <- sim_data$rem[, 3:ncol(sim_data$rem)]
    
    #' Corresponds with "bands_per_sample" in removal.stan
    time_bands_per_sample <- unname(apply(rem,
                                          1,
                                          function (x) sum(!is.na(x))))
    
    #' Total abundance per sampling event.
    #' I.e., this is the sum of Y_ij over j
    #' Corresponds with "abund_per_sample" in removal.stan
    total_abund_per_sample <- unname(apply(rem, 1, function(x) sum(x, na.rm = TRUE)))
    
    # Design matrix, just NULL model for now
    X <- matrix(data = 1, nrow = n, ncol = 1)
    
    #' Corresponds with "abund_per_band" in removal.stan
    abundance_per_band <- rem
    abundance_per_band[is.na(abundance_per_band)] <- 0
    
    #' Corresponds with "max_time" in removal.stan
    time_design <- sim_data$time_design
    time_design[is.na(time_design)] <- 0
    
    stan_data_rem <- list(n_samples = n,
                          n_covariates = ncol(X),
                          max_intervals = ncol(rem),
                          abund_per_band = abundance_per_band,
                          abund_per_sample = total_abund_per_sample,
                          bands_per_sample = time_bands_per_sample,
                          max_time = time_design,
                          X = X)
    
    removal_bayes <- sampling(rem_stan,
                              data = stan_data_rem,
                              verbose = TRUE,
                              chains = 4,
                              iter = 2000,
                              warmup = 1000,
                              cores = 4,
                              pars = c("gamma"),
                              control = list(adapt_delta = 0.8,
                                             max_treedepth = 15))
    save(removal_bayes, file = file_name_rem)
    
    ####### Distance Modelling loop #########################
    
    dis <- sim_data$dis[, 3:ncol(sim_data$dis)]
    
    #' Corresponds with "bands_per_sample" in distance.stan
    dist_bands_per_sample <- unname(apply(dis,
                                          1,
                                          function (x) sum(!is.na(x))))
    
    #' Total abundance per sampling event.
    #' I.e., this is the sum of Y_ik over k
    #' Corresponds with "abund_per_sample" in distance.stan
    total_abund_per_sample <- unname(apply(dis, 1, function(x) sum(x, na.rm = TRUE)))
    
    # Design matrix, just NULL model for now
    X <- matrix(data = 1, nrow = n, ncol = 1)
    
    #' Corresponds with "abund_per_band" in distance.stan
    abundance_per_band <- dis
    abundance_per_band[is.na(abundance_per_band)] <- 0
    
    #' Corresponds with "max_dist" in distance.stan
    dist_design <- sim_data$dist_design
    dist_design[is.na(dist_design)] <- 0
    
    stan_data_dis <- list(n_samples = n,
                          n_covariates = ncol(X),
                          max_intervals = ncol(dis),
                          abund_per_band = abundance_per_band,
                          abund_per_sample = total_abund_per_sample,
                          bands_per_sample = dist_bands_per_sample,
                          max_dist = dist_design,
                          X = X)
    
    distance_bayes <- sampling(dis_stan,
                               data = stan_data_dis,
                               verbose = TRUE,
                               chains = 4,
                               iter = 2000,
                               warmup = 1000,
                               cores = 4,
                               pars = c("theta"),
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15))
    save(distance_bayes, file = file_name_dis)
    
  }
}

