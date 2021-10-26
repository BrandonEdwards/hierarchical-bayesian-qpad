#' Simulate point count data with removal and distance sampling
#' 
#' @param n_obs How many sampling events in this point count data set?
#' @param phi Singing/cue rate
#' @param tau Effective detection radius
#' @param n_protocols Number of point count protocols to simulate
#' @param n_time_bins Number of time bins per protocol (vector)
#' @param n_dis_bins Number of distance bins per protocol (vector)
#' @param max_times Maximum time for each protocol
#' @param max_dist Maximum distance for each protocol
#' @param poisson_lambda Lambda parameter for generating true abundance
#' 

sim_pc <- function(n_obs = 1000,
                   phi = NULL,
                   tau = NULL,
                   n_protocols = NULL,
                   n_time_bins = NULL,
                   n_dist_bins = NULL,
                   max_times = NULL,
                   max_dist = NULL,
                   poisson_lambda = 10)
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
  rem_df <- as.data.frame(matrix(data = NA,
                                 nrow = n_obs,
                                 ncol = length(rem_names)))
  names(rem_df) <- rem_names
  
  dis_names <- c("N", "Protocol",
                 paste0(rep("Int", times = max(n_dist_bins)), 1:max(n_dist_bins)))
  dis_df <- as.data.frame(matrix(data = NA,
                                 nrow = n_obs,
                                 ncol = length(dis_names)))
  names(dis_df) <- dis_names
  
  rem_df$N <- dis_df$N <- round(rpois(n = n_obs, lambda = poisson_lambda)) + 1
  rem_df$Protocol <- dis_df$Protocol <- protocol_num
  
  
  #' Main algorithm
  #' First,
  n <- 1
  chances_counter <- 0
  while (n <= 1000)
  {
    N <- rem_df[n, "N"]
    p <- rem_df[n, "Protocol"]
    max_times_p <- max_times[p, ]
    max_dist_p <- max_dist[p, ]
    
    # First, generate
    available <- vector(mode = "numeric", length = n_time_bins[p])
    dist_total <- vector(mode = "numeric", length = n_dist_bins[p])
    for (j in 1:n_time_bins[p])
    {
      if (j == 1)
      {
        available[j] <- round(N * (1 - exp(-max_times_p[j] * phi)))
      }else
      {
        available[j] <- round(N * (1 - exp(-max_times_p[j] * phi))) - sum(available)
      }
      
      if (available[j] > 0)
      {
        bird_dists <- runif(n = available[j], min = 0, max = max(max_dist_p, na.rm = TRUE))
        binned_birds <- as.data.frame(table(cut(bird_dists, breaks = c(-1, max_dist_p[1:n_dist_bins[p]]))))$Freq
        
        u_r <- 1 - exp(-(max_dist_p[1:n_dist_bins[p]] ^ 2 / tau^2))
        recorded <- round(binned_birds * ((tau^2 / max_dist_p[1:n_dist_bins[p]] ^ 2) * (u_r)))
        rem_df[n, paste0("Int", j)] <- sum(recorded)       
      }else{
        recorded <- rep(0, n_dist_bins[p])
        rem_df[n, paste0("Int", j)] <- sum(recorded)       
      }
      dist_total <- dist_total + recorded
    }
    
    dis_df[n, 2 + c(1:n_dist_bins[p])] <- dist_total
    
    # Check for 0 counts for a survey. If so, throw that away and try again
    if (sum(dist_total) == 0)
    {
      if (chances_counter == 10)
      {
        N <- N + 1
        rem_df[n, "N"] <- N
        dis_df[n, "N"] <- N
        chances_counter <- 0
        n <- n - 1
      }else{
        chances_counter <- chances_counter + 1
        n <- n - 1
      }
    }else{
      chances_counter <- 0
    }

    n <- n + 1
  }
  
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