
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
      step_size_distr <- data.frame('m'=c(0, step_size_distr$mean), 's'=c(1,step_size_distr$sd), 'w'=c(0,1))
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
  if (step_size_distr$w[1] == 0) {
    mode <- rep(2, n_simulations) # for now, just use the above-0 asymptote mode)
  } else {
    mode <- as.integer( runif(n_simulations) > step_size_distr$w[1] ) + 1
  }

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
  
  asymptote_distr_name <- list('aiming' = 'multimodal', 'adapt' = 'normal')[[signal]]
  
  asymptote_distributions <- read.csv(sprintf('data/distributions/%s_exp_asymptote_%s_parameters.csv', signal, asymptote_distr_name), stringsAsFactors = FALSE)
  changerate_distributions <- read.csv(sprintf('data/distributions/%s_exp_changerate_exponential_parameter.csv', signal), stringsAsFactors = FALSE)
  exponential_SD_distributions   <- read.csv(sprintf('data/distributions/%s_exp_sd_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
  # step_SD_distributions   <- step_SD_distributions[is.na(step_SD_distributions$makestep),]
  
  # print(asymptote_distributions)
  # print(changerate_distributions)
  # print(step_SD_distributions)
  
  n_simulations <- 20000
  simulations <- list()
  
  for (rotation in c(20,30,40,50,60)) {
    
    if (signal == 'adapt') {
      asymp_distr <- asymptote_distributions[asymptote_distributions$rotation==rotation,]
      asymp_distr <- data.frame('m'=c(0, asymp_distr$mean), 's'=c(1, asymp_distr$sd), 'w'=c(0,1))
    } else {
      asymp_distr <- asymptote_distributions[asymptote_distributions$rotation==rotation,c('m','s','w')]
    }
    # print(asymp_distr)
    
    changerate_distr <- changerate_distributions[changerate_distributions$rotation==rotation,]
    # print(changerate_distr)
    
    exp_SD_distr   <- exponential_SD_distributions[exponential_SD_distributions$rotation==rotation,c('shape','rate')]
    # print(exp_SD_distr)

    AR <- bootstrapExponentialModel( asymp_distr      = asymp_distr,
                                     changerate_distr = changerate_distr,
                                     exp_SD_distr     = exp_SD_distr,
                                     n_simulations    = n_simulations)
    # print(str(AR))
    
    simulations[[as.character(rotation)]] <- AR

  }

  saveRDS(simulations, file = sprintf('data/simulations/%s_exponential_simulations.rds', signal))
  # return(AR)
  
}

bootstrapExponentialModel <- function(asymp_distr, 
                                      changerate_distr, 
                                      exp_SD_distr, 
                                      n_simulations = 20000) {
  
  trials <- 120
  
  # Create a matrix to store the results
  results <- matrix(0, nrow = n_simulations, ncol = trials)
  
  # ~0 asymptote (1) or above-0 asymptote (2)?
  if (asymp_distr$w[1] == 0) {
    mode <- rep(2, n_simulations) # for now, just use the above-0 asymptote mode)
  } else {
    mode <- as.integer( runif(n_simulations) > asymp_distr$w[1] ) + 1
  }

  # asymptotes:
  asymptotes <- rep(0, n_simulations)
  asymptotes[mode == 1] <- rnorm(sum(mode == 1), mean = asymp_distr$m[1], sd = asymp_distr$s[1])
  asymptotes[mode == 2] <- rnorm(sum(mode == 2), mean = asymp_distr$m[2], sd = asymp_distr$s[2])
  
  # print(length(which(asymptotes <= 0)))
  
  # step sizes (0 for no step):
  change_rates <- rexp(n = n_simulations, rate = changerate_distr$rate)
  # cat(sprintf("percentage change rates > 1: %0.1f\n", 100*sum(change_rates > 1)/length(change_rates))) 
  change_rates[which(change_rates > 1)] <- 1 # cap at 1, otherwise the exponential function will overshoot the asymptote
  
  change_rates[which(change_rates < 0.001)] <- 0.001 # cap at 0, otherwise the exponential function will undershoot the asymptote]
  print(sort(change_rates)[1:100])
  
  rellevels <- mapply(function(r,t) {r^t}, change_rates, matrix(rep(c(0:(trials-1)), each=length(change_rates)), ncol=trials,nrow=length(change_rates)))
  curves <- 1 - matrix(rellevels, nrow=length(change_rates))
  
  # multiply each row of curves with the asymptote of that simulated participant
  curves <- curves * matrix(rep(asymptotes, each=trials), nrow=n_simulations, ncol=trials)
  
  
  # figure out the noise:
  exp_SD <- rgamma(n=n_simulations, shape=exp_SD_distr$shape, rate=exp_SD_distr$rate)
  
  noise <- matrix(rnorm(n=trials*n_simulations,
                        mean=0,
                        sd=rep(step_SD,each=trials)), # same SD for all trials in a simulated participant
                  nrow=n_simulations,
                  ncol=trials,
                  byrow = TRUE)
  
  # add noise to the curves:
  responses <- curves + noise
  
  return(responses)
  
}



# plot and compare fits ----

readData <- function(signal) {
  
  if (signal %in% c('aiming','adapt')) {
  } else {
    cat('signal must be either "aiming" or "adapt"\n')
  }
  
  behavior <- list()
  
  for (rot in c(20,30,40,50,60)) {
    
    rotfiles <- list.files(path = "data/summaries", 
                           pattern = sprintf("SUMMARY_aiming%d", rot))
    
    responses <- matrix(NA, nrow = length(rotfiles), ncol = 120)
    
    for (rotfile_no in c(1:length(rotfiles))) {
      
      rotfile <- rotfiles[[rotfile_no]]
      
      # ppno <- ppno + 1
      
      ppid <- substr(strsplit(rotfile, "_")[[1]][3], 1, 6)
      
      # cat(sprintf('working on participant %d (%s, %d° rotation)\n', ppno, ppid, rot))
      
      # read participant data
      data <- read.csv(sprintf("data/summaries/%s", rotfile))
      
      if (signal == 'aiming') { 
        
        responses[rotfile_no,] <- data$aimdeviation_deg[which(data$rotation_deg == -1 * rot)]
        
      } else if (signal == 'adapt') {
        
        responses[rotfile_no,] <- data$reachdeviation_deg[which(data$rotation_deg == -1 * rot)]
        
      }
      
    }
    
    behavior[[as.character(rot)]] <- responses
    
  }
  
  return(behavior)

}

plotDataAndFits <- function(signal='aiming') {
  
  if (signal %in% c('aiming','adapt')) {
  } else {
    cat('signal must be either "aiming" or "adapt"\n')
  }
  
  behavior <- readData(signal)
  stepfunction <- readRDS(sprintf('data/simulations/%s_stepfunction_simulations.rds', signal))
  exponential <- readRDS(sprintf('data/simulations/%s_exponential_simulations.rds', signal))
  
  layout(mat=matrix(c(1:(3*5)), nrow=5, ncol=3, byrow=TRUE))
  par(mar=c(2,3,1.5,.01))
  
  x <- seq(.5,120.5, by=1)
  y <- seq(-10,70, by=4) 
  yd <- seq(-10,70,by=.5)
  
  bluePal <- colorRampPalette(c("white","blue"))
  blueramp <- bluePal(200)
  
  for (rotation in c(20,30,40,50,60)) {
    
    # plot the data
    plot(NULL,NULL,
         xlim=c(1,120), ylim=c(-10,70), main='',
         xlab='', ylab='', xaxt='n', yaxt='n',
         bty='n'
         )
    if (rotation == 20) {
      title(main=sprintf('%s data', signal), x=60, y=65, cex=1.5)
    }
    title(ylab=sprintf('%s° condition', rotation),line=2)
    axis(side=2, at=seq(0,60,by=20), labels=seq(0,60,by=20), las=2)
    if (rotation == 60) {
      axis(side=1, at=c(1,30,60,90,120))
    } else {
      axis(side=1, at=c(1,30,60,90,120), labels=rep('',5))
    }

    # image(x=x, y=y,
    #       z=t(get2Dcounts(matrix=behavior[[as.character(rotation)]])), 
    #       add=TRUE, axes=FALSE )
    
    image(x=x, y=yd,
          z=t(get2Ddensity(matrix=behavior[[as.character(rotation)]])), 
          add=TRUE, axes=FALSE, col=blueramp)
    
    
    # plot the exponential simulations
    plot(NULL,NULL,
         xlim=c(1,120), ylim=c(-10,70), main='',
         xlab='', ylab='', xaxt='n', yaxt='n',
         bty='n')
    if (rotation == 20) {
      title(main=sprintf('exponential distribution'), x=60, y=65, cex=1.5)
    }
    if (rotation == 60) {
      axis(side=1, at=c(1,30,60,90,120))
    } else {
      axis(side=1, at=c(1,30,60,90,120), labels=rep('',5))
    }
    
    axis(side=2, at=seq(0,60,by=20), labels=rep('',4), las=2)
    
    # image(x=x, y=y,
    #       z=t(get2Dcounts(matrix=exponential[[as.character(rotation)]])), 
    #       add=TRUE, axes=FALSE )
    image(x=x, y=yd,
          z=t(get2Ddensity(matrix=exponential[[as.character(rotation)]])), 
          add=TRUE, axes=FALSE, col=blueramp )
    
    # plot the step-function simulations
    plot(NULL,NULL,
         xlim=c(1,120), ylim=c(-10,70), main='',
         xlab='', ylab='', xaxt='n', yaxt='n',
         bty='n')
    if (rotation == 20) {
      title(main=sprintf('step-function distribution'), x=60, y=65, cex=1.5)
    }
    if (rotation == 60) {
      axis(side=1, at=c(1,30,60,90,120))
    } else {
      axis(side=1, at=c(1,30,60,90,120), labels=rep('',5))
    }
    
    axis(side=2, at=seq(0,60,by=20), labels=rep('',4), las=2)
    
    # image(x=x, y=y,
    #       z=t(get2Dcounts(matrix=stepfunction[[as.character(rotation)]])), 
    #       axes=FALSE, add=TRUE)
    image(x=x, y=yd,
          z=t(get2Ddensity(matrix=stepfunction[[as.character(rotation)]])), 
          axes=FALSE, add=TRUE, col=blueramp)

    
  }
  
}

get2Dcounts <- function(matrix, from=-10, to=70, n=20) {
  
  breaks <- seq(from,to,(to - from)/n)
  
  output <- matrix(NA, nrow=n, ncol=ncol(matrix)) 
  
  for (col_idx in 1:ncol(matrix)) {
    data <- matrix[,col_idx]
    data <- data[data >= from & data < to]
    # print(data)
    # print(hist(data, breaks=breaks, plot=FALSE)$counts)
    output[,col_idx] <- hist(data, breaks=breaks, plot=FALSE)$counts
  }
  
  return(output)
  
}


get2Ddensity <- function(matrix, from=-10, to=70, n=160) {
  
  breaks <- seq(from,to,(to - from)/n)
  
  output <- matrix(NA, nrow=n, ncol=ncol(matrix)) 
  
  for (col_idx in 1:ncol(matrix)) {
    data <- matrix[,col_idx]
    data <- data[data >= from & data < to]
    # output[,col_idx] <- hist(data, breaks=breaks, plot=FALSE)$counts
    output[,col_idx] <- density(data, from=from, to=to, n=n)$y
  }
  
  return(output^.5)
  
}