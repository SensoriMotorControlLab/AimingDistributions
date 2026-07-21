
# initial parameters -----

## OUTDATED functions: -----
### [STEP FUNCTION simulations] -----

generateSimpleStepfunctionDistribution <- function(signal) {
  
  # cat('start simulations...\n')
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
  # cat('simulations done, saving to file...\n')
  saveRDS(simulations, file = sprintf('data/simulations/%s_stepfunction_simulations.rds', signal))
  # cat('done!\n')
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

### [EXPONENTIAL FUNCTION simulations] -----

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
  
  change_rates[which(change_rates < 0.01)] <- 0.01 # cap at 0, otherwise the exponential function will undershoot the asymptote]
  # print(sort(change_rates)[1:100])
  
  rellevels <- mapply(function(r,t) {r^t}, change_rates, matrix(rep(c(0:(trials-1)), each=length(change_rates)), ncol=trials,nrow=length(change_rates)))
  curves <- 1 - matrix(rellevels, nrow=length(change_rates))
  
  # multiply each row of curves with the asymptote of that simulated participant
  curves <- curves * matrix(rep(asymptotes, each=trials), nrow=n_simulations, ncol=trials)
  
  
  # figure out the noise:
  exp_SD <- rgamma(n=n_simulations, shape=exp_SD_distr$shape, rate=exp_SD_distr$rate)
  
  noise <- matrix(rnorm(n=trials*n_simulations,
                        mean=0,
                        sd=rep(exp_SD,each=trials)), # same SD for all trials in a simulated participant
                  nrow=n_simulations,
                  ncol=trials,
                  byrow = TRUE)
  
  # add noise to the curves:
  responses <- curves + noise
  
  return(responses)
  
}



## plot and compare fits ----

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
  
  simulations <- readRDS('data/simulations/initial_model_simulations.rds')
  
  # stepfunction <- readRDS(sprintf('data/simulations/%s_stepfunction_simulations.rds', signal))
  stepfunction <- simulations[[signal]][['stepfunction']][['simulations']]
  # exponential <- readRDS(sprintf('data/simulations/%s_exponential_simulations.rds', signal))
  exponential <- simulations[[signal]][['exponential']][['simulations']]
  
  step_d <- c()
  exp_d  <- c()
  
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


get2Ddensity <- function(matrix, from=-10, to=70, n=160, bw=1) {
  
  breaks <- seq(from,to,(to - from)/n)
  
  output <- matrix(NA, nrow=n, ncol=ncol(matrix)) 
  
  # print(str(matrix))
  
  for (col_idx in 1:ncol(matrix)) {
    data <- matrix[,col_idx]
    # print(summary(data))
    data <- data[data >= from & data < to]
    # output[,col_idx] <- hist(data, breaks=breaks, plot=FALSE)$counts
    output[,col_idx] <- density(data, from=from, to=to, n=n, bw=bw)$y
  }
  
  return(output^.5)
  
}


saveSimulations <- function() {
  
  
  simulations <- list()
  
  for (signal in c('aiming', 'adapt')) {
    
    simulations[[signal]] <- list()
    
    for (model in c('stepfunction', 'exponential')) {
      
      simulations[[signal]][[model]] <- list()
      
      startPars <- getStartingParameters(signal, model)
      
      if (model == 'stepfunction') {
        out <- getStepfunctionSimulations( par   = startPars$pars,
                                           fixed = startPars$fixed,
                                           n_simulations = 20000)
      } else if (model == 'exponential') {
        out <- getExponentialSimulations( par   = startPars$pars,
                                          fixed = startPars$fixed,
                                          n_simulations = 20000)
      }
      
      simulations[[signal]][[model]][['simulations']] <- out
      simulations[[signal]][[model]][['parameters']] <- startPars
      
    }
    
    
  }
  
  saveRDS(simulations, file = 'data/simulations/initial_model_simulations.rds')
  
}


getStepfunctionSimulations <- function(par, fixed=NULL, n_simulations=20000) {
  
  # combine the fixed and free parameters into one list
  if (!is.null(fixed)) {
    all_pars <- c(par, fixed)
  } else {
    all_pars <- par
  }
  
  gen_idx <- grep('all_', names(all_pars))
  gen_pars <- all_pars[gen_idx]
  names(gen_pars) <- substr(names(gen_pars), 5, nchar(names(gen_pars)))
  
  # the simulations use lots of random numbers, which makes the model
  # a bit, well "random"... to make it more stable by setting the rng seed here 
  # that way, the same "random" numbers are used in each run of the model
  set.seed(37331)
  
  out <- list()
  
  for (rotation in c(20,30,40,50,60)) {
    
    rot_idx <- grep(sprintf('r%d_', rotation), names(all_pars))
    rot_pars <- all_pars[rot_idx]
    # rot_pars <- all_pars[grep(sprintf('r%d_', rotation), names(all_pars))]
    
    # cat('read from here:\n')
    names(rot_pars) <- substr(names(rot_pars), 5, nchar(names(rot_pars)))
    # print(rot_pars)
    
    model <- simulateStepfunctionModel( par           = c(gen_pars, rot_pars),
                                        n_simulations = n_simulations)
    
    print(str(model))
    
    out[[as.character(rotation)]] <- model
    
  }
  
  cat('done?\n')
  return(out)
  
}

getExponentialSimulations <- function(par, fixed=NULL, n_simulations=5000) {
  
  # combine the fixed and free parameters into one list
  if (!is.null(fixed)) {
    all_pars <- c(par, fixed)
  } else {
    all_pars <- par
  }
  
  gen_idx <- grep('all_', names(all_pars))
  gen_pars <- all_pars[gen_idx]
  names(gen_pars) <- substr(names(gen_pars), 5, nchar(names(gen_pars)))
  
  # the simulations use lots of random numbers, which makes the model
  # a bit, well "random"... to make it more stable by setting the rng seed here 
  # that way, the same "random" numbers are used in each run of the model
  set.seed(37331)
  
  out <- list()
  
  for (rotation in c(20,30,40,50,60)) {
    
    rot_idx <- grep(sprintf('r%d_', rotation), names(all_pars))
    rot_pars <- all_pars[rot_idx]
    # rot_pars <- all_pars[grep(sprintf('r%d_', rotation), names(all_pars))]
    
    # cat('read from here:\n')
    names(rot_pars) <- substr(names(rot_pars), 5, nchar(names(rot_pars)))
    # print(rot_pars)
    
    model <- simulateExponentialModel( par           = c(gen_pars, rot_pars),
                                        n_simulations = n_simulations)
    
    print(str(model))
    
    out[[as.character(rotation)]] <- model
    
  }
  
  cat('done?\n')
  
  return(out)
  
}

# statistics ----

nll <- function(d) {
  
  # remove values resulting in errors or Infs?
  d <- d[which(!is.na(d))]
  d[which(d < .Machine$double.eps)] <- .Machine$double.eps
  
  # d[which((d-1) == 0)] <- .Machine$double.eps # maybe this should be made closer to 1?
  
  nll <- -1 * sum(log(d))
  
  if (!is.finite(nll)) {
    cat('non-finite nll!\n')
    nll <- -1 * sum(log(rep(.Machine$double.eps, length(d)))) # minimum probability
  }
  
  return(nll)
  
}


calculateAICs <- function() {
  
  for (signal in c('aiming', 'adapt')) {
    
    behavior <- readData(signal)
    # cat('loading data...\n')
    stepfunction <- readRDS(sprintf('data/simulations/%s_stepfunction_simulations.rds', signal))
    exponential <- readRDS(sprintf('data/simulations/%s_exponential_simulations.rds', signal))
    
    step_d <- c()
    exp_d  <- c()
    # cat('calculating likelihoods...\n')
    for (rotation in c(60)) {
      
      data <- behavior[[as.character(rotation)]]
      step_model <- stepfunction[[as.character(rotation)]]
      exp_model  <- exponential[[as.character(rotation)]]
      
      step_densities <- getProbabilityDensities(data, step_model)
      exp_densities  <- getProbabilityDensities(data, exp_model)
      
      step_d <- c(step_d, step_densities)
      exp_d  <- c(exp_d, exp_densities)
      
    }
    # cat('calculating AICs...\n')
    step_k <- 5 * (5 + 2 + 2)
    exp_k  <- 5 * (5 + 1 + 2)
    
    step_AIC <- Reach::AIC(logLik = -1*Reach::nll(step_d), k = step_k, N=length(step_d))
    exp_AIC  <- Reach::AIC(logLik = -1*Reach::nll(exp_d ), k = exp_k,  N=length(exp_d ))
    
    cat(sprintf('%s: step-function AIC = %.2f, exponential AIC = %.2f\n', signal, step_AIC, exp_AIC))
    
  }
  
}

getProbabilityDensities <- function(data, model, bw=.5) {
  
  if (ncol(data) != ncol(model)) {
    cat('data and model must have the same number of columns\n')
    return(NA)
  }
  
  # probd <- c()
  # for (col_idx in 1:ncol(data)) {
  #   data_col <- data[,col_idx]
  #   model_col <- model[,col_idx]
  #   
  #   dens <- density(model_col, 
  #                   from=min(model_col), 
  #                   to=max(model_col), 
  #                   n=250) # does n matter here?
  #   # function(xs, t, h = bw.nrd0(xs)) mean(dnorm(t, mean = xs, sd = h))
  #   probd <- c(probd, dnorm(data_col, mean = dens$x, sd = dens$bw))
  # }
  
  col_idx <- c(1:ncol(data))
  
  probd <- lapply(col_idx, function(idx) {
    # print(idx)
    # print(as.double(data[,idx]))
    data_col <- data[,idx]
    model_col <- model[,idx]
    # print(summary(model_col))
    # print(as.double(model_col))
    
    data_col <- data_col[!is.na(data_col)]
    model_col <- model_col[!is.na(model_col)]
    
    # if (!is.null(noise)) {
    #   bw = noise
    #   # mx <- seq(min(model_col), max(model_col), length.out=1000)
    #   
    # } else {
    #   dens <- density(model_col,
    #                   from=min(model_col),
    #                   to=max(model_col),
    #                   n=1000) # does n matter here?
    #   bw = dens$bw
    #   mx = dens$x
    # }
    
    # cat('got density\n')
    # function(xs, t, h = bw.nrd0(xs)) mean(dnorm(t, mean = xs, sd = h))
    # return(dnorm(data_col, mean = mx, sd = bw))
    
    # return(dnorm(data_col, model_col, sd = bw))
    # a <- lapply(data_col, function(x) mean(dnorm(x, mean = model_col, sd = bw)))
    
    # return(rowMeans( matrix(unlist(lapply(data_col, function(x) mean(dnorm(x, mean = model_col, sd = bw))) ), 
    #                         nrow=length(data_col), 
    #                         ncol=length(model_col), 
    #                         byrow=TRUE) ))
    
    # changed to be the probability of the model given the data
    # with a fixed bandwidth (bw) for the kernel density estimation
    return( colMeans( dnorm( matrix(data_col,  ncol=length(data_col), nrow=length(model_col), byrow=TRUE),
                             matrix(model_col, ncol=length(data_col), nrow=length(model_col), byrow=FALSE),
                             sd = bw) ))
    # with using means of probability of each data point given all model simulations,
    # the nll should become independent of the number of simulations used
    
  })
  
  # print(probd)
  
  return(unlist(probd))
  
}

# FIT the models? -----

initialFits <- function() {
  
  signal <- c()
  model  <- c()
  nll    <- c()
  AIC    <- c()
  
  for (sig in c('aiming', 'adapt')) {
    behavior <- readData(sig)
    for (mod in c('stepfunction', 'exponential')) {
      startPar <- getStartingParameters(sig,mod)
      
      if (mod == 'stepfunction') {
        out <- NLLstepfunctionModel(data  = behavior, 
                                    par   = startPar$pars, 
                                    fixed = startPar$fixed,
                                    n_simulations=20000)
      } else if (mod == 'exponential') {
        out <- NLLexponentialModel(data  = behavior, 
                                   par   = startPar$pars, 
                                   fixed = startPar$fixed,
                                   n_simulations=20000)
      }
      signal <- c(signal, sig)
      model  <- c(model, mod)
      nll    <- c(nll, out)
      AIC    <- c(AIC, Reach::AIC(logLik = -1*out, k = length(startPar$pars), N=120*200))
    }
  }
  
  write.csv(data.frame(signal=signal, model=model, nll=nll),
            file='data/fits/initial_model_fits.csv', row.names=FALSE)
  
}


# optimize the model parameters? -----
# we could still reduce some of the parameters:

fitModels <- function() {
  
  for (signal in c('aiming', 'adapt')) {
    for (model in c('stepfunction', 'exponential')) {
      out <- fitDistributionModel(signal, model)
      cat(sprintf('fitted %s model to %s data\n', model, signal))
      print(out)
      saveRDS(out, file = sprintf('data/fits/%s_%s_fit.rds', signal, model))
    }
  }
  
}

fitDistributionModel <- function(signal, model) {
  
  if (signal %in% c('aiming','adapt')) {
  } else {
    cat('signal must be either "aiming" or "adapt"\n')
  }
  
  if (model %in% c('stepfunction','exponential')) {
  } else {
    cat('model must be either "stepfunction" or "exponential"\n')
  }
  
  behavior <- readData(signal)
  
  startPar <- getStartingParameters(signal,model)
  
  if (model == 'stepfunction') {
    fitpar <- fitStepfunctionModel(data  = behavior, 
                                   par   = startPar$pars, 
                                   fixed = startPar$fixed,
                                   lower = startPar$lower,
                                   upper = startPar$upper)
  } else if (model == 'exponential') {
    fitpar <- fitExponentialModel(data  = behavior, 
                                  par   = startPar$pars, 
                                  fixed = startPar$fixed,
                                  lower = startPar$lower,
                                  upper = startPar$upper)
  }
  
  return(fitpar)

}

getStartingParameters <- function(signal, model) {
  
  if (signal %in% c('aiming','adapt')) {
  } else {
    cat('signal must be either "aiming" or "adapt"\n')
  }
  
  if (model %in% c('stepfunction','exponential')) {
  } else {
    cat('model must be either "stepfunction" or "exponential"\n')
  }
  
  pars <- list()
  lower <- c()
  upper <- c()
  fixed <- list()
  
  if (model == 'stepfunction') {
    
    # single parameters:
    #  - step time (gamma)
    
    # per rotation parameters:
    #  - step size (normal or multimodal)
    #  - noise SD (mean of the fitted gamma distributions)
    
    step_time_distributions <- read.csv(sprintf('data/distributions/%s_step_time_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
    step_time_distr <- step_time_distributions[step_time_distributions$rotation==20,c('shape','rate')]
    
    pars['all_steptime_rate'] = step_time_distr$rate
    lower <- c(lower, .0001)
    upper <- c(upper, Inf)
    pars['all_steptime_shape'] = step_time_distr$shape
    lower <- c(lower, .1)
    upper <- c(upper, Inf)
    
    
    for (rotation in c(20,30,40,50,60)) {
      
      step_size_distr_name <- list('aiming' = 'multimodal', 'adapt' = 'normal')[[signal]]
      step_size_distributions <- read.csv(sprintf('data/distributions/%s_step_size_%s_parameters.csv', signal, step_size_distr_name), stringsAsFactors = FALSE)
      step_SD_distributions   <- read.csv(sprintf('data/distributions/%s_step_SD_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
      
      if (signal == 'adapt') {
        step_size_distr <- step_size_distributions[step_size_distributions$rotation==rotation,]
        asymp_distr <- data.frame('m'=c(0, step_size_distr$mean), 's'=c(1,step_size_distr$sd), 'w'=c(0,1))
      } else {
        asymp_distr <- step_size_distributions[step_size_distributions$rotation==rotation,c('m','s','w')]
      }
      step_SD_distr   <- step_SD_distributions[which(step_SD_distributions$rotation==rotation & is.na(step_SD_distributions$makestep)),c('shape','rate')]
      
      if (signal == 'adapt') {
        
        fixed[sprintf('r%d_asymp_m0',       rotation)] = asymp_distr$m[1]
        fixed[sprintf('r%d_asymp_s0',       rotation)] = asymp_distr$s[1]
        fixed[sprintf('r%d_asymp_w0',       rotation)] = asymp_distr$w[1]
        
        pars[sprintf('r%d_asymp_m1',       rotation)] = asymp_distr$m[2]
        lower <- c(lower, 0)
        upper <- c(upper, rotation+10)
        pars[sprintf('r%d_asymp_s1',       rotation)] = asymp_distr$s[2]
        lower <- c(lower, .0001)
        upper <- c(upper, Inf)
        
        fixed[sprintf('r%d_asymp_w1',       rotation)] = asymp_distr$w[2]
        
      } else {
        fixed[sprintf('r%d_asymp_m0',       rotation)] = asymp_distr$m[1]
        
        pars[ sprintf('r%d_asymp_s0',       rotation)] = asymp_distr$s[1]
        lower <- c(lower, .0001)
        upper <- c(upper, Inf)
        
        fixed[sprintf('r%d_asymp_w0',       rotation)] = asymp_distr$w[1]
        fixed[sprintf('r%d_asymp_m1',       rotation)] = asymp_distr$m[2]
        
        pars[ sprintf('r%d_asymp_s1',       rotation)] = asymp_distr$s[2]
        lower <- c(lower, .0001)
        upper <- c(upper, Inf)
        
        fixed[sprintf('r%d_asymp_w1',       rotation)] = asymp_distr$w[2]
      }
      
      pars[sprintf('r%d_noise',          rotation)] = step_SD_distr$shape/step_SD_distr$rate
      lower <- c(lower, .1)
      upper <- c(upper, Inf)
    }
    

  } else if (model == 'exponential') {
    
    # single parameters:
    #  - change rate (gamma NOT exponential)
    
    # rotation specific parameters:
    # - asymptote (normal or multimodal)
    # - noise SD (mean of the fitted gamma distributions) 
    
    changerate_distributions <- read.csv(sprintf('data/distributions/%s_exp_changerate_gamma_parameter.csv', signal), stringsAsFactors = FALSE)
    
    changerate_distr <- changerate_distributions[changerate_distributions$rotation==20,]
    
    pars['all_roc_rate'] = changerate_distr$rate
    lower <- c(lower, .001)
    upper <- c(upper, Inf)
    pars['all_roc_shape'] = changerate_distr$shape
    lower <- c(lower, 1.001)
    upper <- c(upper, Inf)
    
    
    for (rotation in c(20,30,40,50,60)) {
      
      asymptote_distr_name <- list('aiming' = 'multimodal', 'adapt' = 'normal')[[signal]]
      asymptote_distributions <- read.csv(sprintf('data/distributions/%s_exp_asymptote_%s_parameters.csv', signal, asymptote_distr_name), stringsAsFactors = FALSE)
      exponential_SD_distributions   <- read.csv(sprintf('data/distributions/%s_exp_sd_gamma_parameters.csv', signal), stringsAsFactors = FALSE)
      
      if (signal == 'adapt') {
        asymp_distr <- asymptote_distributions[asymptote_distributions$rotation==rotation,]
        asymp_distr <- data.frame('m'=c(0, asymp_distr$mean), 's'=c(1, asymp_distr$sd), 'w'=c(0,1))
      } else {
        asymp_distr <- asymptote_distributions[asymptote_distributions$rotation==rotation,c('m','s','w')]
      }
      
      exp_SD_distr   <- exponential_SD_distributions[exponential_SD_distributions$rotation==rotation,c('shape','rate')]
      
      # startPar <- list(asymp_distr=asymp_distr,
      #                  changerate_distr=changerate_distr,
      #                  exp_SD_distr=exp_SD_distr)
      

      # fixed <- list( fixed,
      #                sprintf('r%d_asymp_m0', rotation) = asymp_distr$m[1],
      #                sprintf('r%d_asymp_s0', rotation) = asymp_distr$s[1],
      #                sprintf('r%d_asymp_w0', rotation) = asymp_distr$w[1])
      # pars <- list(  pars, 
      #                sprintf('r%d_asymp_m1', rotation) = asymp_distr$m[2],
      #                sprintf('r%d_asymp_s1', rotation) = asymp_distr$s[2],
      #                sprintf('r%d_asymp_w1', rotation) = asymp_distr$w[2],
      #                sprintf('r%d_roc_rate', rotation) = changerate_distr$rate,
      #                sprintf('r%d_noise', rotation)    = exp_SD_distr$shape/exp_SD_distr$rate)
      
      
      if (signal == 'adapt') {
        fixed[sprintf('r%d_asymp_m0', rotation)] = asymp_distr$m[1]
        fixed[sprintf('r%d_asymp_s0', rotation)] = asymp_distr$s[1]
        fixed[sprintf('r%d_asymp_w0', rotation)] = asymp_distr$w[1]
        pars[sprintf('r%d_asymp_m1', rotation)] = asymp_distr$m[2]
        lower <- c(lower, 0)
        upper <- c(upper, rotation+10)
        pars[sprintf('r%d_asymp_s1', rotation)] = asymp_distr$s[2]
        lower <- c(lower, .0001)
        upper <- c(upper, Inf)
        fixed[sprintf('r%d_asymp_w1', rotation)] = asymp_distr$w[2]
        
      } else {
        fixed[sprintf('r%d_asymp_m0', rotation)] = asymp_distr$m[1]
        pars[ sprintf('r%d_asymp_s0', rotation)] = asymp_distr$s[1]
        lower <- c(lower, .0001)
        upper <- c(upper, Inf)
        fixed[sprintf('r%d_asymp_w0', rotation)] = asymp_distr$w[1]
        fixed[sprintf('r%d_asymp_m1', rotation)] = asymp_distr$m[2]
        pars[ sprintf('r%d_asymp_s1', rotation)] = asymp_distr$s[2]
        lower <- c(lower, .0001)
        upper <- c(upper, Inf)
        fixed[sprintf('r%d_asymp_w1', rotation)] = asymp_distr$w[2]
        
      }
      
      pars[sprintf('r%d_noise', rotation)]    = exp_SD_distr$shape/exp_SD_distr$rate
      lower <- c(lower, .1)
      upper <- c(upper, Inf)
      
    }
    
  }
  
  return(list( 'pars'  = unlist(pars), 
               'fixed' = unlist(fixed), 
               'lower' = lower, 
               'upper' = upper))
  
}

### stepfunction specific ----

fitStepfunctionModel <- function(data, par, fixed, lower=NULL, upper=NULL) {
  
  
  # nlm     WORKS (no constraints)
  # nlminb  WORKS (box constraints)
  # spg     doesn't work (has 'simple' constraints)
  # ucminf  WORKS (no constraints) comes from separate package
  # newuoa  is for unconstrained fitting
  # bobyqa  WORKS (and has constraints)
  # nmkb    doesn't work (package not installed)
  # hjkb    doesn't work (package not installed)
  # Rcgmin  WORKS... or not? (supposedly has constraints, from separate package)
  # Rvmmin  WORKS? maybe (not sure if it has constraints... probably from separate package... not available for this version of R)
  
  # methods <- c("L-BFGS-B", "nlminb", "bobyqa") #, "Rcgmin", "Rvmmin")
  
  methods <- "L-BFGS-B"
  
  if (!is.null(lower) & !is.null(upper)) {
    model <- optimx::optimx(  par = par, 
                              NLLstepfunctionModel,
                              method = methods,
                              lower = lower,
                              upper = upper,
                              data = data,
                              fixed = fixed
    )
  } else {
    model <- optim( par = par, 
                    NLLstepfunctionModel,
                    method = "BFGS",
                    data = data,
                    fixed = fixed
    )
  }
  
  
  
  return(list('model'=model, 'fixedpar'=fixed))
  
}

NLLstepfunctionModel <- function(par, data, fixed=NULL, n_simulations=5000) {
  
  # combine the fixed and free parameters into one list
  if (!is.null(fixed)) {
    all_pars <- c(par, fixed)
  } else {
    all_pars <- par
  }

  probdens <- c()
  
  # print(str(data))
  
  gen_idx <- grep('all_', names(all_pars))
  gen_pars <- all_pars[gen_idx]
  names(gen_pars) <- substr(names(gen_pars), 5, nchar(names(gen_pars)))
  
  # the simulations use lots of random numbers, which makes the model
  # a bit, well "random"... to make it more stable by setting the rng seed here 
  # that way, the same "random" numbers are used in each run of the model
  set.seed(37331)
  
  for (rotation in c(20,30,40,50,60)) {
    
    rot_idx <- grep(sprintf('r%d_', rotation), names(all_pars))
    rot_pars <- all_pars[rot_idx]
    # rot_pars <- all_pars[grep(sprintf('r%d_', rotation), names(all_pars))]
    
    # cat('read from here:\n')
    names(rot_pars) <- substr(names(rot_pars), 5, nchar(names(rot_pars)))
    # print(rot_pars)
    
    model <- simulateStepfunctionModel( par           = c(gen_pars, rot_pars),
                                        n_simulations = n_simulations)
    
    # noise <- c(gen_pars, rot_pars)['noise']
    
    # print(str(model))
    probdens <- c(probdens, 
                  getProbabilityDensities(data[[sprintf('%d',rotation)]], 
                                          model))
    
  }
  
  # negLogLik <- Reach::nll(probdens)
  negLogLik <- nll(probdens)
  cat(sprintf('negLogLik: %.2f\n', negLogLik))
  return(negLogLik)
  
}

simulateStepfunctionModel <- function( par, n_simulations = 20000) {
  
  # unpack the parameters:
  asymp_m1 <- par['asymp_m0']
  asymp_s1 <- par['asymp_s0']
  asymp_w1 <- par['asymp_w0']
  
  asymp_m2 <- par['asymp_m1']
  asymp_s2 <- par['asymp_s1']
  asymp_w2 <- par['asymp_w1']
  
  steptime_rate  <- par['steptime_rate']
  steptime_shape <- par['steptime_shape']
  
  noise          <- par['noise']
  
  

  trials <- 120 # should this be a parameter as well?
  
  # Create a matrix to store the results
  results <- matrix(0, nrow = n_simulations, ncol = trials)
  
  # step or no step?
  if (asymp_w1 == 0) {
    mode <- rep(2, n_simulations) # for now, just use the above-0 asymptote mode)
  } else {
    mode <- as.integer( runif(n_simulations) > asymp_w1 ) + 1
  }
  
  # step sizes (0 for no step):
  step_sizes <- rep(0, n_simulations)
  if (any(mode == 1)) {
    step_sizes[mode == 1] <- rnorm(sum(mode == 1), mean = asymp_m1, sd = asymp_s1)
  }
  # step_sizes[mode == 1] <- 0
  step_sizes[mode == 2] <- rnorm(sum(mode == 2), mean = asymp_m2, sd = asymp_s2)
  
  # step times (NA for no step - should not be used later on, which will throw errors):
  step_times <- rep(NA, n_simulations)
  
  # this was to catch a few NAs from rgamma, but apparently,
  # there are either none, or all of them are NA
  # N_times <- length(step_times[mode == 2])
  # m2st    <- rep(NA, N_times) 
  # while (any(is.na(m2st))) {
  #   m2st[is.na(m2st)] <- ceiling(rgamma( n = sum(is.na(m2st)), 
  #                                        shape = steptime_shape, 
  #                                        rate = steptime_rate))
  # }
  # 
  # step_times[mode == 2] <- m2st
  
  step_times[mode == 2] <- ceiling(rgamma( n=sum(mode == 2),
                                           shape = steptime_shape,
                                           rate  = steptime_rate ) )
  
  # if gamma returns NAs, we make the step time, the latest possible
  # and increase the noise, so that the model can still fit the dat
  # but returns low likelihoods, so that the optimizer can find a better solution
  if (any(is.na(step_times[mode == 2]))) {
    step_times[mode == 2] <- trials
    noise <- 999999
  }
  
  step_times[which(step_times > trials)] <- trials
  
  # # simple model has just one level of noise throughout:
  # step_SD <- rgamma(n=n_simulations, shape=step_SD_distr$shape, rate=step_SD_distr$rate)
  # 
  # rand_noise <- matrix( rnorm(n=trials*n_simulations,
  #                             mean=0,
  #                             sd=rep(step_SD,each=trials)), # same SD for all trials in a simulated participant 
  #                       nrow=n_simulations,
  #                       ncol=trials,
  #                       byrow = TRUE) 
  
  
  rand_noise <- matrix( rnorm(n=trials*n_simulations,
                              mean=0,
                              sd=noise), # same SD for all trials in a simulated participant
                        nrow=n_simulations,
                        ncol=trials,
                        byrow = TRUE)
  
  
  # steps are added in a loop... can't think of a better way right now
  if (any(mode == 2)){
    for (idx in which(mode == 2)) {
      # print(step_times[idx])
      # cat(sprintf('idx: %d, step_time: %d, step_size: %.2f\n', idx, step_times[idx], step_sizes[idx]))
      results[idx,c(max(1, step_times[idx]):trials)] <- step_sizes[idx]
    }
  }
  # noise is added in one go:
  responses <- results + rand_noise
  
  return(responses)
  
}

### exponential specific -----

fitExponentialModel <- function(data, par, fixed, lower=NULL, upper=NULL) {
  
  # methods <- c("nlm", "nlminb", "ucminf", "newuoa", "bobyqa")
  methods <- "L-BFGS-B"
  
  if (!is.null(lower) & !is.null(upper)) {
    out <- optim( par = par, 
                  NLLexponentialModel,
                  method = "L-BFGS-B",
                  lower = lower,
                  upper = upper,
                  data = data,
                  fixed = fixed
    )
  } else {
    out <- optim( par = par, 
                  NLLexponentialModel,
                  method = "BFGS",
                  data = data,
                  fixed = fixed
    )
  }
  
  out$fixed = fixed
  
  return(out)
  
}

NLLexponentialModel <- function(par, data, fixed=NULL, n_simulations=5000) {
  
  # combine the fixed and free parameters into one list
  if (!is.null(fixed)) {
    all_pars <- c(par, fixed)
  } else {
    all_pars <- par
  }
  
  probdens <- c()
  
  gen_idx <- grep('all_', names(all_pars))
  gen_pars <- all_pars[gen_idx]
  names(gen_pars) <- substr(names(gen_pars), 5, nchar(names(gen_pars)))
  
  set.seed(37331)
  
  for (rotation in c(20,30,40,50,60)) {
    
    rot_pars <- all_pars[grep(sprintf('r%d_', rotation), names(all_pars))]
    
    # cat('read from here:\n')
    names(rot_pars) <- substr(names(rot_pars), 5, nchar(names(rot_pars)))
    
    model <- simulateExponentialModel( par           = c(gen_pars, rot_pars),
                                       n_simulations = n_simulations)
    
    probdens <- c(probdens, 
                  getProbabilityDensities(data[[sprintf('%d',rotation)]], 
                                          model))
    
  }
  
  # negLogLik <- Reach::nll(probdens)
  negLogLik <- nll(probdens)
  cat(sprintf('negLogLik: %.2f\n', negLogLik))
  return(negLogLik)
  
}

simulateExponentialModel <- function( par, n_simulations = 20000) {
 
  # unpack the parameters:
  asymp_m1  <- par['asymp_m0']
  asymp_s1  <- par['asymp_s0']
  asymp_w1  <- par['asymp_w0']
  
  asymp_m2  <- par['asymp_m1']
  asymp_s2  <- par['asymp_s1']
  asymp_w2  <- par['asymp_w1']
  
  roc_rate  <- par['roc_rate']
  roc_shape <- par['roc_shape']

  noise     <- par['noise']
  
  trials <- 120
  
  # Create a matrix to store the results
  results <- matrix(0, nrow = n_simulations, ncol = trials)
  
  # ~0 asymptote (1) or above-0 asymptote (2)?
  if (asymp_w1 == 0) {
    mode <- rep(2, n_simulations) # for now, just use the above-0 asymptote mode)
  } else {
    mode <- as.integer( runif(n_simulations) > asymp_w1 ) + 1
  }
  
  # asymptotes:
  asymptotes <- rep(0, n_simulations)
  asymptotes[mode == 1] <- rnorm(sum(mode == 1), mean = asymp_m1, sd = asymp_s1)
  asymptotes[mode == 2] <- rnorm(sum(mode == 2), mean = asymp_m2, sd = asymp_s2)
  
  # print(length(which(asymptotes <= 0)))
  
  # step sizes (0 for no step):
  # change_rates <- rexp(n = n_simulations, rate = roc_rate)
  change_rates <- rgamma(n = n_simulations, shape = roc_shape, rate = roc_rate)
  # cat(sprintf("percentage change rates > 1: %0.1f\n", 100*sum(change_rates > 1)/length(change_rates))) 
  change_rates[which(change_rates > 1)] <- 1 # cap at 1, otherwise the exponential function will overshoot the asymptote
  
  change_rates[which(change_rates < 0.01)] <- 0.01 # min at 0.01, otherwise the exponential function will undershoot the asymptote]
  # print(sort(change_rates)[1:100])
  
  rellevels <- mapply(function(r,t) {r^t}, change_rates, matrix(rep(c(0:(trials-1)), each=length(change_rates)), ncol=trials,nrow=length(change_rates)))
  curves <- 1 - matrix(rellevels, nrow=length(change_rates))
  
  # multiply each row of curves with the asymptote of that simulated participant
  curves <- curves * matrix(rep(asymptotes, each=trials), nrow=n_simulations, ncol=trials)
  
  
  # figure out the noise:
  # exp_SD <- rgamma(n=n_simulations, shape=exp_SD_distr$shape, rate=exp_SD_distr$rate)
  # 
  # noise <- matrix(rnorm(n=trials*n_simulations,
  #                       mean=0,
  #                       sd=rep(exp_SD,each=trials)), # same SD for all trials in a simulated participant
  #                 nrow=n_simulations,
  #                 ncol=trials,
  #                 byrow = TRUE)
  
  rand_noise <- matrix( rnorm(n    = trials*n_simulations,
                              mean = 0,
                              sd   = noise), # same SD for all participants and trials
                        nrow=n_simulations,
                        ncol=trials,
                        byrow = TRUE)
  
  # add noise to the curves:
  responses <- curves + rand_noise
  
  return(responses)
   
}