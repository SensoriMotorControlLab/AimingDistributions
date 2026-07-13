
StepWiseSimulations <- function() {
  
  step_size_distributions <- read.csv('data/step_size_multimodal_parameters.csv', stringsAsFactors = FALSE)
  step_time_distributions <- read.csv('data/step_time_gamma_parameters.csv', stringsAsFactors = FALSE)
  step_SD_distributions   <- read.csv('data/step_SD_gamma_parameters.csv', stringsAsFactors = FALSE)
  
  for (rotation in c(20,30,40,50,60)) {
    
    step_size_distr <- step_size_distributions[step_size_distributions$rotation==rotation,]
    step_time_distr <- step_time_distributions[step_time_distributions$rotation==rotation,]
    step_SD_distr   <- step_SD_distributions[step_SD_distributions$rotation==rotation,]
    
    bootstrapStepWiseModel(step_size_distr, step_time_distr, step_SD_distr, n_simulations = 15000)
    
  }
  
}

bootstrapStepWiseModel <- function(step_size_dist, step_time_dist, step_SD_dist, n_simulations = 15000) {
  
  trials <- 120
  # Create a matrix to store the results
  results <- matrix(NA, nrow = n_simulations, ncol = trials)
  
  mode <- as.integer( runif(n_simulations) > step_size_dist$w[2] ) + 1
  
  step_sizes <- rep(0, n_simulations)
  step_sizes[mode == 1] <- 0
  step_sizes[mode == 2] <- rnorm(sum(mode == 2), mean = step_size_dist$m[2], sd = step_size_dist$s[2])
  
  step_times <- rep(NA, n_simulations)
  step_times[mode == 2] <- rgamma(n_simulations, shape = step_time_dist$shape, scale = step_time_dist$scale)

}