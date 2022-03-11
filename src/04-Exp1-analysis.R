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

  ### Univariate Plots ############

phi_boxplot <- ggplot(data = ext_data,
                      aes(x = Sample_Size, y = Phi_Est)) +
  geom_boxplot(aes(fill = Model)) +
  geom_hline(yintercept = phi_true, linetype="dashed", color = "red") +
  xlab("Sample Size") + 
  ylab("Cue Rate Estimate") +
  scale_fill_viridis(discrete = TRUE) +
  NULL

phi_sd_boxplot <- ggplot(data = ext_data,
						 aes(x = Sample_Size, y = Phi_SD)) +
	geom_boxplot(aes(fill = Model)) +
	xlab("Sample Size") +
	ylab("Cue Rate Standard Deviation") +
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

tau_sd_boxplot <- ggplot(data = ext_data,
						 aes(x = Sample_Size, y = Tau_SD)) +
	geom_boxplot(aes(fill = Model)) +
	xlab("Sample Size") +
	ylab("EDR Standard Deviation") +
	scale_fill_viridis(discrete = TRUE) +
	NULL

  ### Bivariate Estimates vs Standard Deviation plots #######

phi_vs_sd_list <- vector(mode = "list", length = 6)
tau_vs_sd_list <- vector(mode = "list", length = 6)

i <- 1
for (n in unique(ext_data$Sample_Size))
{
	for (m in unique(ext_data$Model))
	{
		to_plot <- ext_data[which(ext_data$Sample_Size == n &
									ext_data$Model == m), ]

		phi_vs_sd_list[[i]] <- ggplot(data = to_plot, aes(x = Phi_Est, y = Phi_SD)) +
									  geom_point() +
									  xlab("Cue Rate Estimate") +
									  ylab("Cue Rate Standard Deviation") +
									  #xlim(min_phi - 0.1, max_phi + 0.1) +
									  #ylim(min_phi_sd - 0.005, max_phi_sd + 0.005) +
									  NULL

		tau_vs_sd_list[[i]] <- ggplot(data = to_plot, aes(x = Tau_Est, y = Tau_SD)) +
							  			geom_point() +
							 			xlab("EDR Estimate") +
							  			ylab("EDR Standard Deviation") +
							  			NULL
		i <- i + 1
	}
}

	### Bivariate Estimates vs Estimate plots #######

	# Incomplete, come back to this.

phi_plot_list <- vector(mode = "list", length = 3)
tau_plot_list <- vector(mode = "list", length = 3)

i <- 1
for (n in unique(ext_data$Sample_Size))
{
	to_plot <- ext_data[which(ext_data$Sample_Size == n), ]

}


####### Calculate Performance Measures ############

# Bias

####### Write Results #############################

png(file = "output/exp1/plots/univariate.png", width = 8, height = 8, res = 300, units = "in")
ggarrange(phi_boxplot, phi_sd_boxplot, tau_boxplot, tau_sd_boxplot,
		  labels = c("A", "B", "C", "D"),
		  common.legend = TRUE,
		  legend = "top")
dev.off()

png(file = "output/exp1/plots/bivariate_phi_vs_sd.png", width = 8, height = 8, res = 300, units = "in")
ggarrange(plotlist = phi_vs_sd_list, nrow = 3, ncol = 2)
dev.off()

png(file = "output/exp1/plots/bivariate_tau_vs_sd.png", width = 8, height = 8, res = 300, units = "in")
ggarrange(plotlist = tau_vs_sd_list, nrow = 3, ncol = 2)
dev.off()
