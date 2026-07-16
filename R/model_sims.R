
# STEP FUNCTION simulations -----

generateSimpleStepfunctionDistribution <- function(signal) {
  
  if (signal %in% c('aiming','adapt')) {
  } else {
    cat('signal must be either "aiming" or "adapt"\n')
  }
  
  step_size_distr_name <- list('aiming' = 'multimodal', 'adapt' = 'normal')[[signal]]
  
  step_size_distributions <- read.csv(sprintf('data/distributions/%s_step_size_%s_parameters.csv', signal, step_size_distr_name), stringsAsFactors = FALSE)
  step_time_distributions <- read.csv(sprintf('data/distributions/%s_step_time_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
  step_SD_distributions   <- read.csv(sprintf('data/distributions/%s_step_SD_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
  step_SD_distributions   <- step_SD_distributions[is.na(step_SD_distributions$makestep),]
  

  # in order to have 1000 runs even in the 7.5% participants without a strategy in the 60 degree rotation
  # we need to run 15000 simulations in total (1000/0.075 = 13333.33... round up to be sure)
  
  # in the exponential model of aiming responses, it's only 5.4% though, which means 18k
  # simulations, which we round up to 20k
  n_simulations <- 20000
  simulations <- list()
  
  for (rotation in c(20,30,40,50,60)) {
    
    if (signal == 'adapt') {
      step_size_distr <- step_size_distributions[step_size_distributions$rotation==rotation,]
      step_size_distr <- data.frame('m'=c(step_size_distr$mean, 0), 's'=c(step_size_distr$sd,1), 'w'=c(1,0))
    } else {
      step_size_distr <- step_size_distributions[step_size_distributions$rotation==rotation,c('m','s','w')]
    }
    step_time_distr <- step_time_distributions[step_time_distributions$rotation==rotation,c('shape','rate')]
    step_SD_distr   <- step_SD_distributions[step_SD_distributions$rotation==rotation,c('shape','rate')]
    
    AR <- bootstrapStepWiseModel(step_size_distr = step_size_distr, 
                                 step_time_distr = step_time_distr, 
                                 step_SD_distr   = step_SD_distr, 
                                 n_simulations   = n_simulations)
    
    simulations[[as.character(rotation)]] <- AR
    
  }
  
  saveRDS(simulations, file = sprintf('data/simulations/%s_stepfunction_simulations.rds', signal))
  
}

bootstrapStepWiseModel <- function(step_size_distr, step_time_distr, step_SD_distr, n_simulations = 20000) {
  
  trials <- 120
  
  # Create a matrix to store the results
  results <- matrix(0, nrow = n_simulations, ncol = trials)
  
  # step or no step?
  mode <- as.integer( runif(n_simulations) > step_size_distr$w[2] ) + 1
  
  # step sizes (0 for no step):
  step_sizes <- rep(0, n_simulations)
  step_sizes[mode == 1] <- 0
  step_sizes[mode == 2] <- rnorm(sum(mode == 2), mean = step_size_distr$m[2], sd = step_size_distr$s[2])
  
  # step times (NA for no step - should not be used later on, which will throw errors):
  step_times <- rep(NA, n_simulations)
  step_times[mode == 2] <- ceiling(rgamma( n=sum(mode == 2), 
                                           shape = step_time_distr$shape, 
                                           rate = step_time_distr$rate ) )
  
  step_times[which(step_times > trials)] <- trials
  
  # simple model has just one level of noise throughout:
  step_SD <- rgamma(n=n_simulations, shape=step_SD_distr$shape, rate=step_SD_distr$rate)
  
  noise <- matrix(rnorm(n=trials*n_simulations,
                        mean=0,
                        sd=rep(step_SD,each=trials)), # same SD for all trials in a simulated participant 
                  nrow=n_simulations,
                  ncol=trials,
                  byrow = TRUE) 
  
  # steps are added in a loop... can't think of a better way right now
  for (idx in which(mode == 2)) {
    # cat(sprintf('idx: %d, step_time: %d, step_size: %.2f\n', idx, step_times[idx], step_sizes[idx]))
    results[idx,c(max(1,step_times[idx]):trials)] <- step_sizes[idx]
  }
  
  # noise is added in one go:
  aiming_responses <- results + noise
  
  return(aiming_responses)
  
}

# EXPONENTIAL FUNCTION simulations -----

generateExponentialFunctionDistribution <- function(signal) {
  
  if (signal %in% c('aiming','adapt')) {
  } else {
    cat('signal must be either "aiming" or "adapt"\n')
  }
  
  exp_asymptote_distr_name <- list('aiming' = 'multimodal', 'adapt' = 'normal')[[signal]]
  
  exp_asymptote_distributions <- read.csv(sprintf('data/distributions/%s_step_size_%s_parameters.csv', signal, step_size_distr_name), stringsAsFactors = FALSE)
  step_time_distributions <- read.csv(sprintf('data/distributions/%s_step_time_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
  step_SD_distributions   <- read.csv(sprintf('data/distributions/%s_step_SD_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
  step_SD_distributions   <- step_SD_distributions[is.na(step_SD_distributions$makestep),]
  
  
  # in order to have 1000 runs even in the 7.5% participants without a strategy in the 60 degree rotation
  # we need to run 15000 simulations in total (1000/0.075 = 13333.33... round up to be sure)
  n_simulations <- 20000
  simulations <- list()
  
  for (rotation in c(20,30,40,50,60)) {
    
    if (signal == 'adapt') {
      step_size_distr <- step_size_distributions[step_size_distributions$rotation==rotation,]
      step_size_distr <- data.frame('m'=c(step_size_distr$mean, 0), 's'=c(step_size_distr$sd,1), 'w'=c(1,0))
    } else {
      step_size_distr <- step_size_distributions[step_size_distributions$rotation==rotation,c('m','s','w')]
    }
    step_time_distr <- step_time_distributions[step_time_distributions$rotation==rotation,c('shape','rate')]
    step_SD_distr   <- step_SD_distributions[step_SD_distributions$rotation==rotation,c('shape','rate')]
    
    AR <- bootstrapExponentialModel(step_size_distr = step_size_distr, 
                                 step_time_distr = step_time_distr, 
                                 step_SD_distr   = step_SD_distr, 
                                 n_simulations   = n_simulations)
    
    simulations[[as.character(rotation)]] <- AR
    
  }
  
  saveRDS(simulations, file = sprintf('data/simulations/%s_stepfunction_simulations.rds', signal))
  
}

bootstrapExponentialModel <- function(step_size_distr, step_time_distr, step_SD_distr, n_simulations = 20000) {
  
  trials <- 120
  
  # Create a matrix to store the results
  results <- matrix(0, nrow = n_simulations, ncol = trials)
  
  # step or no step?
  mode <- as.integer( runif(n_simulations) > step_size_distr$w[2] ) + 1
  
  # step sizes (0 for no step):
  step_sizes <- rep(0, n_simulations)
  step_sizes[mode == 1] <- 0
  step_sizes[mode == 2] <- rnorm(sum(mode == 2), mean = step_size_distr$m[2], sd = step_size_distr$s[2])
  
  # step times (NA for no step - should not be used later on, which will throw errors):
  step_times <- rep(NA, n_simulations)
  step_times[mode == 2] <- ceiling(rgamma( n=sum(mode == 2), 
                                           shape = step_time_distr$shape, 
                                           rate = step_time_distr$rate ) )
  
  step_times[which(step_times > trials)] <- trials
  
  # simple model has just one level of noise throughout:
  step_SD <- rgamma(n=n_simulations, shape=step_SD_distr$shape, rate=step_SD_distr$rate)
  
  noise <- matrix(rnorm(n=trials*n_simulations,
                        mean=0,
                        sd=rep(step_SD,each=trials)), # same SD for all trials in a simulated participant 
                  nrow=n_simulations,
                  ncol=trials,
                  byrow = TRUE) 
  
  # steps are added in a loop... can't think of a better way right now
  for (idx in which(mode == 2)) {
    # cat(sprintf('idx: %d, step_time: %d, step_size: %.2f\n', idx, step_times[idx], step_sizes[idx]))
    results[idx,c(max(1,step_times[idx]):trials)] <- step_sizes[idx]
  }
  
  # noise is added in one go:
  aiming_responses <- results + noise
  
  return(aiming_responses)
  
}