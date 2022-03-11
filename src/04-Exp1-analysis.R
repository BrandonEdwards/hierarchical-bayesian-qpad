####### Script Information ########################
# Brandon P.M. Edwards
# Hierarchical QPAD Simulation Study
# 04-Exp1-analysis.R
# Created March 2022
# Last Updated March 2022

####### Import Libraries and External Files #######

library(ggplot2)
library(ggpubr)
library(viridis)
theme_set(theme_pubclean())

####### Set Constants #############################

phi_true <- 0.4
tau_true <- 100

####### Read Data #################################

ext_data <- read.csv("output/exp1/extracted_data.csv")

####### Clean Data ################################

# Remove erronous data but save to data frame for further inspection
error_data <- ext_data[which(ext_data$Phi_Est < 0.001), ]
error_data <- rbind(error_data, ext_data[which(ext_data$Tau_Est > 1000), ])

ext_data <- ext_data[-which(ext_data$Phi_Est < 0.001), ]
ext_data <- ext_data[-which(ext_data$Tau_Est > 1000), ]

# Factor the sample sizes
ext_data$Sample_Size <- factor(ext_data$Sample_Size,
                               levels = c("75", "500", "1000"))

####### Exploratory Plots #########################

  ## Univariate Plots ############

phi_boxplot <- ggplot(data = ext_data,
                      aes(x = Sample_Size, y = Phi_Est)) +
  geom_boxplot(aes(fill = Model)) +
  geom_hline(yintercept = phi_true, linetype="dashed", color = "red") +
  xlab("Sample Size") + 
  ylab("Cue Rate Estimate") +
  scale_fill_viridis(discrete = TRUE) +
  NULL

tau_boxplot <- ggplot(data = ext_data,
                      aes(x = Sample_Size, y = Tau_Est)) +
  geom_boxplot(aes(fill = Model)) +
  geom_hline(yintercept = tau_true, linetype="dashed", color = "red") +
  xlab("Sample Size") + 
  ylab("EDR Estimate") +
  scale_fill_viridis(discrete = TRUE) +
  NULL

 ### Bivariate plots #######



####### Calculate Performance Measures ############

# Bias



