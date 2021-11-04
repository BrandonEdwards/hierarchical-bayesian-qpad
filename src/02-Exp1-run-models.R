####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 02-Exp1-run-models.R
# Created October 2021
# Last Updated November 2021

####### Import Libraries and External Files #######

library(detect)
library(doParallel)
library(foreach)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


####### Set Constants #############################

n_cores <- 30
n_obs <- c(75,500, 1000)
n_sim <- 400

# phi <- 0.4
# tau <- 100
# den <- 10
# n_protocols <- 4
# n_time_bins <- c(3, 4, 10, 7)
# n_dist_bins <- c(3, 4, 7, 11)
# max_times <- matrix(data = c(3, 4, 5, NA, NA, NA, NA, NA, NA, NA, 
#                              2, 4, 5, 6, NA, NA, NA, NA, NA, NA,
#                              1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
#                              2, 3, 4, 5, 6, 7, 8, NA, NA, NA),
#                     nrow = 4, byrow = TRUE)
# max_dist <- matrix(data = c(50, 100, 400, NA, NA, NA, NA, NA, NA, NA, NA,
#                             25, 50, 100, 150, NA, NA, NA, NA, NA, NA, NA, 
#                             25, 50, 75, 100, 125, 150, 400, NA, NA, NA, NA, 
#                             10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150),
#                    nrow = 4, byrow = TRUE)


####### Maximum Likelihood ########################

cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)

for (n in n_obs)
{
  foreach (i = 1:n_sim, .packages = 'detect') %dopar%
    {
      file_name_sim <- paste0("output/exp1/sim_data/n",
                              n,
                              "-",
                              i,
                              ".rda")
      file_name_dis <- paste0("output/exp1/dis_mods/n",
                              n,
                              "-",
                              i,
                              "_dis_mle.rda")
      file_name_rem <- paste0("output/exp1/rem_mods/n",
                              n,
                              "-",
                              i,
                              "_rem_mle.rda")
      x <- load(file_name_sim)
      
      removal_mle <- cmulti(as.matrix(x$rem[, 3:ncol(x$rem)]) | as.matrix(x$time_design) ~ 1, type = "rem")
      save(removal_mle, file = file_name_rem)
      
      distance_mle <- cmulti(as.matrix(x$dis[, 3:ncol(x$dis)]) | as.matrix(x$dist_design) ~ 1, type = "dis")
      save(distance_mle, file = file_name_dis)
      
    }
}

stopCluster(cluster)

####### Bayesian Model (Stan) #####################


rem_stan <- stan_model(file = "stan/removal.stan")
dis_stan <- stan_model(file = "stan/distance.stan")

cluster <- makeCluster(10, type = "PSOCK")
registerDoParallel(cluster)

for (n in n_obs)
{
  foreach (i = 1:n_sim, .packages = 'stan') %dopar%
  {
    file_name_sim <- paste0("output/exp1/sim_data/n",
                            n,
                            "-",
                            i,
                            ".rda")
    file_name_dis <- paste0("output/exp1/dis_mods/n",
                            n,
                            "-",
                            i,
                            "_dis_bayes.rda")
    file_name_rem <- paste0("output/exp1/rem_mods/n",
                            n,
                            "-",
                            i,
                            "_rem_bayes.rda")
    x <- load(file_name_sim)
    
    ############ Removal Modelling loop ################
    rem <- x$rem[, 3:ncol(x$rem)]
    
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
    time_design <- x$time_design
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
                              chains = 3,
                              iter = 2000,
                              warmup = 1000,
                              cores = 3,
                              pars = c("gamma"),
                              control = list(adapt_delta = 0.8,
                                             max_treedepth = 15))
    save(removal_bayes, file = file_name_rem)
    
    ####### Distance Modelling loop #########################
    
    dis <- x$dis[, 3:ncol(x$dis)]
    
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
    dist_design <- x$dist_design
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
                               chains = 3,
                               iter = 2000,
                               warmup = 1000,
                               cores = 3,
                               pars = c("theta"),
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15))
    save(distance_bayes, file = file_name_dis)
    
  }
}

stopCluster(cluster)
