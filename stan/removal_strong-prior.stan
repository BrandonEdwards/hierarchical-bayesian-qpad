data {
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 1> n_covariates;        // total number of covariates
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int abund_per_band[n_samples, max_intervals];// abundance in time band j for sample i
  vector[n_samples] abund_per_sample; // total abundnace for sample i
  int bands_per_sample[n_samples]; // number of time bands for sample i
  matrix[n_samples, max_intervals] max_time; // max time duration for time band j
  matrix[n_samples, n_covariates] X;  // matrix of covariates
}

parameters {
  vector[n_covariates] gamma;
}

model {
  vector[n_samples] log_phi;             // singing rate
  matrix[n_samples, max_intervals] Pi;   // probabilities
  
  gamma ~ normal(-1, 0.5);
  log_phi = X * gamma;
  
  Pi = rep_matrix(0, n_samples, max_intervals);
  
  for (i in 1:n_samples)
  {
    for (j in 2:bands_per_sample[i])
    {
      Pi[i,j] = (exp(-max_time[i,j-1] * exp(log_phi[i])) - 
                 exp(-max_time[i,j] * exp(log_phi[i]))) / 
                (1 - exp(-max_time[i,bands_per_sample[i]] * exp(log_phi[i])));
    }
    Pi[i,1] = 1 - sum(Pi[i,]);
    
    abund_per_band[i,] ~ multinomial(to_vector(Pi[i,]));
  }

}
