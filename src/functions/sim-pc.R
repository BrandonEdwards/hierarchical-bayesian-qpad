#' Simulate point count data with removal and distance sampling
#' 
#' @param n_obs How many sampling events in this point count data set?
#' @param phi Singing/cue rate
#' @param tau Effective detection radius
#' @param den Density for bSims
#' @param n_protocols Number of point count protocols to simulate
#' @param n_time_bins Number of time bins per protocol (vector)
#' @param n_dis_bins Number of distance bins per protocol (vector)
#' @param max_times Maximum time for each protocol
#' @param max_dist Maximum distance for each protocol
#' @param n_cores 2D vector of cores for each protocol and each replicate
#' 

sim_pc <- function(n_obs = 1000,
                   phi = NULL,
                   tau = NULL,
                   den = NULL,
                   n_protocols = NULL,
                   n_time_bins = NULL,
                   n_dist_bins = NULL,
                   max_times = NULL,
                   max_dist = NULL,
                   n_cores = NULL)
{
  # Randomly allocate sample size per protocol
  protocol_prop <- round(n_obs * rdirichlet(n = 1, alpha = rep(1, n_protocols)))
  protocol_num <- NULL
  for (i in 1:n_protocols)
  {
    protocol_num <- c(protocol_num, rep(i, protocol_prop[i]))
  }

  # Due to rounding, ensure that the total sample size remains == n_obs
  if (length(protocol_num) > n_obs)
  {
    protocol_num <- protocol_num[1:n_obs]
  }else if (length(protocol_num) < n_obs)
  {
    remainder <- n_obs - length(protocol_num)
    to_add <- rep(4, remainder)
    protocol_num <- c(protocol_num, to_add)
  }
  
  rem_names <- c("N", "Protocol",
                 paste0(rep("Int", times = max(n_time_bins)), 1:max(n_time_bins)))
  rem_df_list <- vector(mode = "list", length = n_protocols)
  
  dis_names <- c("N", "Protocol",
                 paste0(rep("Int", times = max(n_dist_bins)), 1:max(n_dist_bins)))
  dis_df_list <- vector(mode = "list", length = n_protocols)
  
  for (p in 1:n_protocols)
  {
    n_obs_p <- unname(table(protocol_num)[p])
    
    rem_df_list[[p]] <- data.frame(as.data.frame(matrix(data = NA,
                                                        nrow = n_obs_p,
                                                        ncol = length(rem_names))))
    names(rem_df_list[[p]]) <- rem_names
    
    dis_df_list[[p]] <- as.data.frame(matrix(data = NA,
                                             nrow = n_obs_p,
                                             ncol = length(dis_names)))
    names(dis_df_list[[p]]) <- dis_names
    
    rem_df_list[[p]]$Protocol <- dis_df_list[[p]]$Protocol <- p
  }

  sim_landscape <- bsims_all(density = den, vocal_rate = phi, tau = tau/100)
  
  cluster_protocols <- makeCluster(n_cores[1], type = "PSOCK")
  registerDoParallel(cluster_protocols)
  
  foreach(p = 1:n_protocols, .packages = 'bSims') %dopar%
  {
    n_replicates <- nrow(rem_df_list[[p]])
    
    cluster_reps <- makeCluster(n_cores[2], type = "PSOCK")
    sim_reps <- sim_landscape$replicate(n_replicates, cl = cluster_reps)
    stopCluster(cluster_reps)
    
    tint <- max_times[p, 1:n_time_bins[p]]
    rint <- max_dist[p, 1:n_dist_bins[p]] / 100 
    
    for (n in 1:n_replicates)
    {
      tr <- bsims_transcribe(sim_reps[[n]], tint=tint, rint=rint)
      tally <- tr$removal
      dis_df[n, 2 + c(1:n_dist_bins[p])] <- unname(rowSums(tally))
      rem_df[n, 2 + c(1:n_time_bins[p])] <- unname(colSums(tally))      
    }
  }
  stopCluster(cluster_protocols)
  
  rem_df <- do.call(rbind, rem_df_list)
  dis_df <- do.call(rbind, dis_df_list)
  
  time_design <- matrix(data = NA,
                        nrow = n_obs,
                        ncol = ncol(max_times))
  time_design <- data.frame(max_times[protocol_num, ])
  names(time_design) <- paste0(rep("Int", times = ncol(max_times)), 1:ncol(max_times))
  
  dist_design <- matrix(data = NA,
                        nrow = n_obs,
                        ncol = ncol(max_dist))
  dist_design <- data.frame(max_dist[protocol_num, ])
  names(dist_design) <- paste0(rep("Int", times = ncol(max_dist)), 1:ncol(max_dist))
  
  return(list(rem = rem_df,
              dis = dis_df,
              dist_design = dist_design,
              time_design = time_design))
}
#'   
#'   
#'   #' Main algorithm
#'   #' First,
#'   n <- 1
#'   chances_counter <- 0
#'   while (n <= n_obs)
#'   {
#'     N <- rem_df[n, "N"]
#'     p <- rem_df[n, "Protocol"]
#'     max_times_p <- max_times[p, ]
#'     max_dist_p <- max_dist[p, ]
#'     
#'     # First, generate
#'     available <- vector(mode = "numeric", length = n_time_bins[p])
#'     dist_total <- vector(mode = "numeric", length = n_dist_bins[p])
#'     for (j in 1:n_time_bins[p])
#'     {
#'       if (j == 1)
#'       {
#'         available[j] <- round(N * (1 - exp(-max_times_p[j] * phi)))
#'       }else
#'       {
#'         available[j] <- round(N * (1 - exp(-max_times_p[j] * phi))) - sum(available)
#'       }
#'       
#'       if (available[j] > 0)
#'       {
#'         bird_dists <- runif(n = available[j], min = 0, max = max(max_dist_p, na.rm = TRUE))
#'         binned_birds <- as.data.frame(table(cut(bird_dists, breaks = c(-1, max_dist_p[1:n_dist_bins[p]]))))$Freq
#'         
#'         u_r <- 1 - exp(-(max_dist_p[1:n_dist_bins[p]] / tau) ^ 2)
#'         recorded <- round(binned_birds * ((tau^2 / max_dist_p[1:n_dist_bins[p]] ^ 2) * (u_r)))
#'         rem_df[n, paste0("Int", j)] <- sum(recorded)       
#'       }else{
#'         recorded <- rep(0, n_dist_bins[p])
#'         rem_df[n, paste0("Int", j)] <- sum(recorded)       
#'       }
#'       dist_total <- dist_total + recorded
#'     }
#'     
#'     dis_df[n, 2 + c(1:n_dist_bins[p])] <- dist_total
#'     
#'     # Check for 0 counts for a survey. If so, throw that away and try again
#'     if (sum(dist_total) == 0)
#'     {
#'       if (chances_counter == 10)
#'       {
#'         N <- N + 1
#'         rem_df[n, "N"] <- N
#'         dis_df[n, "N"] <- N
#'         chances_counter <- 0
#'         n <- n - 1
#'       }else{
#'         chances_counter <- chances_counter + 1
#'         n <- n - 1
#'       }
#'     }else{
#'       chances_counter <- 0
#'     }
#' 
#'     n <- n + 1
#'   }
#'   
#'   time_design <- matrix(data = NA,
#'                         nrow = n_obs,
#'                         ncol = ncol(max_times))
#'   time_design <- data.frame(max_times[protocol_num, ])
#'   names(time_design) <- paste0(rep("Int", times = ncol(max_times)), 1:ncol(max_times))
#'   
#'   dist_design <- matrix(data = NA,
#'                         nrow = n_obs,
#'                         ncol = ncol(max_dist))
#'   dist_design <- data.frame(max_dist[protocol_num, ])
#'   names(dist_design) <- paste0(rep("Int", times = ncol(max_dist)), 1:ncol(max_dist))
#'                                
#'   return(list(rem = rem_df,
#'               dis = dis_df,
#'               dist_design = dist_design,
#'               time_design = time_design))
#' }