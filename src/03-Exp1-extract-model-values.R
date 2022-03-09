####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 03-Exp1-extract-model-values.R
# Created March 2022
# Last Updated March 2022

####### Import Libraries and External Files #######

library(detect)
library(rstan)

####### Set Constants #############################

n_obs <- c(75,500, 1000)
n_sim <- 400
models <- c("mle", "bayes")

####### Extract Data ##############################

# Create empty data frame
ext_data <- data.frame(Sample_Size = rep(n_obs, each = n_sim * length(models)),
                       Sim_Number = rep(seq(1,n_sim), times = length(n_obs) * length(models)),
                       Model = rep(rep(models, each = n_sim), times = length(n_obs)),
                       Phi_Est = NA,
                       Phi_SD = NA,
                       Tau_Est = NA,
                       Tau_SD = NA)

for (n in n_obs)
{
  for (sim in 1:n_sim)
  {
    # Removal MLE model
    load(paste0("output/exp1/rem_mods/n",
                n,
                "-",
                sim,
                "_rem_mle.rda"))
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "mle"),
             "Phi_Est"] <- exp(coef(removal_mle)[1])
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "mle"),
             "Phi_SD"] <- sqrt(vcov(removal_mle))
    
    # Distance MLE Model
    load(paste0("output/exp1/dis_mods/n",
                n,
                "-",
                sim,
                "_dis_mle.rda"))
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "mle"),
             "Tau_Est"] <- exp(coef(distance_mle)[1])
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "mle"),
             "Tau_SD"] <- sqrt(vcov(distance_mle))
    
    # Removal Bayes model
    load(paste0("output/exp1/rem_mods/n",
                n,
                "-",
                sim,
                "_rem_bayes.rda"))
    gamma_list <- extract(removal_bayes)$gamma
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "bayes"),
             "Phi_Est"] <- exp(mean(gamma_list))
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "bayes"),
             "Phi_SD"] <- sqrt(var(gamma_list))
    
    # Distance Bayes model
    load(paste0("output/exp1/dis_mods/n",
                n,
                "-",
                sim,
                "_dis_bayes.rda"))
    theta_list <- extract(distance_bayes)$theta
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "bayes"),
             "Tau_Est"] <- exp(mean(theta_list))
    ext_data[which(ext_data$Sample_Size == n &
                     ext_data$Sim_Number == sim &
                     ext_data$Model == "bayes"),
             "Tau_SD"] <- sqrt(var(theta_list))
  }
}

write.table(ext_data,
            file = "output/exp1/extracted_data.csv",
            sep = ",",
            row.names = FALSE)