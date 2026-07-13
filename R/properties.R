
# generic / all properties ----

extractEmpiricalProperties <- function() {
  
  
  ppno <- 0

  participant    <- c()
  rotation       <- c()
  
  aiming_exp_asymptote  <- c()
  aiming_exp_changerate <- c()
  aiming_exp_sd         <- c()
  
  aiming_step_time      <- c()
  aiming_step_size      <- c()
  aiming_prestep_sd     <- c()
  aiming_poststep_sd    <- c()
  
  adapt_exp_asymptote  <- c()
  adapt_exp_changerate <- c()
  adapt_exp_sd         <- c()
  
  adapt_step_time      <- c()
  adapt_step_size      <- c()
  adapt_prestep_sd     <- c()
  adapt_poststep_sd    <- c()
  
  
  aiming_final_strategy <- c()
  aiming_stratdev_onset <- c()
  aiming_strat_stable   <- c()
  aiming_devel_duration <- c()
  aiming_predev_sd      <- c()
  aiming_devel_sd       <- c()
  aiming_stable_sd      <- c()
  
  for (rot in c(20,30,40,50,60)) {
    
    rotfiles <- list.files(path = "data/summaries", 
                           pattern = sprintf("SUMMARY_aiming%d", rot))
    
    for (rotfile in rotfiles) {
      
      ppid <- substr(strsplit(rotfile, "_")[[1]][3], 1, 6)
      # take first 6 characters of string:
      # ppid <- substr(fnend, 1, 6)
      # ppid <- strsplit(strsplit(rotfile, ".")[[1]][1], "_")[[1]][3]
      # print (ppid)
      
      # read participant data
      data <- read.csv(sprintf("data/summaries/%s", rotfile))
      
      # aiming
      ARtimecourse <- data$aimdeviation_deg[which(data$rotation_deg == -1 * rot)]
      
      # plot(ARtimecourse, main=ppid)
      # 
      # "d8706f" maybe a mirror strategy? or something different for targets on the left and right?
      # "74cdcd"
      # "6a544a"
      # "acd861" (but only up to trial 75 or so)
      # "ab3b79" (very short double strategy)
      # "83456b" (short double strategy)
      # "1e2a6b" (short double strategy)
      
      # "9812d3" 20-30 degrees strategy followed by no strategy at all... what's going on?
      
      # "096819" exploration followed by no strategy at all
      # "5a3568"
      # "57481d" (maybe a little)
      # "514daa"
      # "1ab537" one trial at -8... could be counted as non-strategy - but let's make sure?
      # "43c6f3"
      
      # "3e3a73" 30 degree strategy declines to 0 then a 15-20 degree strategy?
      # "3a78d0" similar but: 0 - 10 - 0 - 8 - 0...
      # "8d426d" relapse to zero around trial 60
      
      # "657fba" 6 ramps of strategy increases
      
      # "5e13b6" radnom responses in a fairly large range (-20 to +15 degrees)
      # "e50c54" no strategy, ends in exploration (-5 to +10 degrees)
      
      # "9dd5aa" negative strategy, followed by no strategy, followed by slow ramp up to almost 100 degrees, followed by double strategy
      
      ntrials <- 8
      
      participant    <- c(participant, ppid)
      rotation       <- c(rotation, rot)
      
      # EXPONENTIAL AIMING TIMECOURSE PROPERTIES
      
      expfit <- Reach::exponentialFit(signal=ARtimecourse,
                                      gridpoints=9, gridfits=5,
                                      asymptoteRange=c(-10,rotation+10))
      
      # print(expfit)
      
      aiming_exp_asymptote  <- c(aiming_exp_asymptote, expfit['N0'])
      aiming_exp_changerate <- c(aiming_exp_changerate, expfit['lambda'])
      expmodel <- Reach::exponentialModel(par=expfit, timepoints=length(ARtimecourse))
      # print(str(expmodel))
      aiming_exp_sd         <- c(aiming_exp_sd, sd(ARtimecourse - expmodel$output, na.rm=TRUE))
      
      
      
      # STEPWISE AIMING TIMECOURSE PROPERTIES
      
      step_df <- data.frame('trial'=c(1:length(ARtimecourse)), 'deviation'=ARtimecourse)
      step_par <- stepFit(data=step_df, gridpoints=6, gridfits=4)
      # print(step_par)
      
      if (step_par['s']>=5) {
        
        aiming_step_time   <- c(aiming_step_time,   step_par['t'])
        aiming_step_size   <- c(aiming_step_size,   step_par['s']) 
        
        # use max, so that we do not have an index lower than 0, and min to ensure no more than 8 trials are used?
        aiming_prestep_sd  <- c(aiming_prestep_sd,  sd(ARtimecourse[1:min(8,max(1, floor(step_par['t'])))]))
        
        # use min, so that we do not have an index higher than the length of the timecourse
        # print( min((length(ARtimecourse)-1),max(1,ceiling(step_par['t']))) )
        aiming_poststep_sd <- c(aiming_poststep_sd, sd(ARtimecourse[min((length(ARtimecourse)-1),max(1,ceiling(step_par['t']))) :  length(ARtimecourse)]))
      } else {
        aiming_step_time   <- c(aiming_step_time,   NA)
        aiming_step_size   <- c(aiming_step_size,   step_par['s'])
        aiming_prestep_sd  <- c(aiming_prestep_sd,  sd(ARtimecourse[1:8]))
        aiming_poststep_sd <- c(aiming_poststep_sd, NA)
      }
      
      # EXPANDED STEPWISE AIMING PROPERTIES
      
      final_strat    <- median(ARtimecourse[c((length(ARtimecourse)-ntrials+1):length(ARtimecourse))])
      aiming_final_strategy <- c(aiming_final_strategy, final_strat)
      
      onset <- which(ARtimecourse > 5)[1]
      if (final_strat < 5) {
        onset <- NA
      }
      aiming_stratdev_onset <- c(aiming_stratdev_onset, onset)
      
      if (is.na(onset)) {
        stab_trial <- NA
      } else {
        stab_tc <- ARtimecourse[onset:length(ARtimecourse)]
        stab_trial <- Reach::findStabilizationTrial(stab_tc)+onset
      }
      aiming_strat_stable <- c(aiming_strat_stable, stab_trial)
      
      aiming_devel_duration <- c(aiming_devel_duration, stab_trial - onset)
      
      if (is.na(onset)) {
        aiming_predev_sd      <- c(aiming_predev_sd, sd(ARtimecourse)) # whole timecourse?
        aiming_devel_sd       <- c(aiming_devel_sd,  NA)
        aiming_stable_sd      <- c(aiming_stable_sd, NA)
      } else if (is.na(stab_trial)) {
        aiming_predev_sd      <- c(aiming_predev_sd, sd(ARtimecourse[1:onset])) # only until onset? when does this ever happen?
        aiming_devel_sd       <- c(aiming_devel_sd,  NA)
        aiming_stable_sd      <- c(aiming_stable_sd, NA)
      } else {
        aiming_predev_sd      <- c(aiming_predev_sd, sd(ARtimecourse[1:onset]))
        devdur <- (stab_trial - onset) + 1
        mdevstrat <- (c(1:devdur)/devdur) * final_strat
        aiming_devel_sd       <- c(aiming_devel_sd, sd(ARtimecourse[onset:stab_trial]-mdevstrat))
        aiming_stable_sd      <- c(aiming_stable_sd, sd(ARtimecourse[stab_trial:length(ARtimecourse)]))
      }
      
      # only get steps for the people with a final strategy larger than 5 deg?
      # no, let's have it depend on stepsize... should strongly correlate though


    }
 
  }
  
  
  aiming_exp_prop_sd <- data.frame(
                                    participant,
                                    rotation,
                                    aiming_exp_asymptote,
                                    aiming_exp_changerate,
                                    aiming_exp_sd
  )
  
  write.csv(aiming_exp_prop_sd,
            file='data/aiming_exp_properties.csv',
            row.names=FALSE, quote=TRUE)

  aiming_stepwise_prop_df <- data.frame(
                                        participant           = participant,
                                        rotation              = rotation,
                                        aiming_step_time      = aiming_step_time,
                                        aiming_step_size      = aiming_step_size,
                                        aiming_prestep_sd     = aiming_prestep_sd,
                                        aiming_poststep_sd    = aiming_poststep_sd    
  )
  
  write.csv(aiming_stepwise_prop_df, 
            file='data/aiming_stepwise_properties.csv', 
            row.names=FALSE, quote=TRUE)
    
  aiming_expanded_prop_df <- data.frame(
                                        participant    = participant,
                                        rotation       = rotation,
                                        aiming_final_strategy = aiming_final_strategy,
                                        aiming_stratdev_onset = aiming_stratdev_onset,
                                        aiming_strat_stable   = aiming_strat_stable,
                                        aiming_devel_duration = aiming_devel_duration,
                                        aiming_predev_sd      = aiming_predev_sd,
                                        aiming_devel_sd       = aiming_devel_sd,
                                        aiming_stable_sd      = aiming_stable_sd
  )
  
  write.csv(aiming_expanded_prop_df, 
            file='data/aiming_expanded_properties.csv', 
            row.names=FALSE, quote=TRUE)
  
}

getProperties <- function() {
  
  if (file.exists('data/empirical_properties.csv')) {
    prop_df <- read.csv('data/empirical_properties.csv')
  } else {
    extractEmpiricalProperties()
    prop_df <- read.csv('data/empirical_properties.csv')
  }
  
  return(prop_df)
  
}

# moved to Reach package:

# # this is the robust mean method found in this paper:
# # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006501
# 
# findStabilizationTrial <- function(timecourse, minn=10) {

plotPropertyDistributionsByRotation <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- extractEmpiricalProperties()
  }
  
  colnames <- c('final_strategy', 'stratdev_onset', 'strat_stable', 'devel_duration', 'predev_sd', 'devel_sd', 'stable_sd', 'step_time', 'step_size', 'prestep_sd', 'poststep_sd')
  
  par(mar=c(4,3,.2,.2))
  ncols <- 3
  layout(matrix(1:(ncols*(ceiling(length(colnames)/ncols))), nrow=ceiling(length(colnames)/ncols), ncol=ncols, byrow=TRUE))
  
  for (colname in colnames) {
    valrange <- range(properties[, colname], na.rm=TRUE)
    plot(NA, 
         xlim=valrange, ylim=c(0.5,6.5),
         xlab=colname, ylab='density', 
         main='',axes=FALSE,bty='n')
    X <- seq(valrange[1], valrange[2], length.out=201)
    for (rot_idx in c(1,2,3,4,5)) {
      rotation <- c(20,30,40,50,60)[rot_idx]
      propvals <- properties[which(properties$rotation == rotation), colname]
      pvd <- density(propvals, na.rm=TRUE, 
                     n = length(X), from=min(X), to=max(X))
      lines(pvd$x, (pvd$y/max(pvd$y))+rot_idx, col=rot_idx)
      points(propvals, rep(rot_idx, length(propvals)), col=rot_idx, pch=20, cex=0.5)
      
    }
    axis(side=1,at=pretty(valrange),labels=pretty(valrange))
    axis(side=2,at=c(1,2,3,4,5),labels=c(20,30,40,50,60))
  }
} 


# step model ----

stepFunction <- function(par, trials) {
  
  if (length(trials) == 1) {
    trials <- c(0:(trials-1))
  }
  
  predictions <- rep(0, length(trials))
  
  predictions[which(trials >= par['t'])] <- par['s']
  
  # print(predictions)
  
  return(predictions)
  
}

stepMSE <- function(par, data) {
  
  # trials <- unique(data$trial)
  # 
  # predictions <- stepFunction(par    = par,
  #                             trials = trials)
  # 
  # errors <- data$deviation - predictions
  # 
  # MSE <- mean(errors^2, na.rm=TRUE)
  # 
  # return(MSE)
  
  return(mean((stepFunction(par = par, trials = unique(data$trial)) - data$deviation)^2, na.rm=TRUE))
  # return(mean((stepFunction(par = par, trials = data$trial) - data$deviation)^2, na.rm=TRUE))
  
}

require('Reach')

stepFit <- function(data, gridpoints=9, gridfits=5) {
  
  # set the search grid:
  parvals <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  # stepsizerange <- diff(range(data$deviation, na.rm=TRUE)))
  stepsizerange <- c(-10, 70)
  stepsizes <- parvals * (diff(stepsizerange) - min(stepsizerange))
  
  # steptimemax <- max(data$trial, na.rm=TRUE)
  steptimemax <- 100
  steptimes <- parvals * steptimemax
  
  searchgrid <- expand.grid('t' = steptimes,
                            's' = stepsizes)
  
  MSE <- apply(searchgrid, FUN=stepMSE, MARGIN=c(1), data=data)
  # print(MSE)
  
  # print(data$deviation)
  
  # lo <- c(0,min(data$deviation, na.rm=TRUE))
  # hi <- c(max(data$trial, na.rm=TRUE), max(data$deviation, na.rm=TRUE))
  
  
  lo <- c(1, -10)
  hi <- c(100, 70)
  
  # print(lo)
  # print(hi)
  
  # print(data.frame(searchgrid[order(MSE)[1:gridfits],]))
  
  # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( data.frame(searchgrid[order(MSE)[1:gridfits],]),
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=stepMSE,
                            
                            # do we use the bounded approach?
                            method     = 'L-BFGS-B',
                            lower      = lo,
                            upper      = hi,
                            
                            data       = data) )
  
  # pick the best fit:
  win <- allfits[order(allfits$value)[1],]
  winpar <- unlist(win[1:2])
  
  # return the best parameters:
  return(winpar)
  
}

# step-size and step-time -----


plotStepPars <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  cat('percentage of people who do not have a strategy:\n(FALSE == strategy, TRUE == no strategy)\n')
  strat_by_rot_table <- table(properties$rotation, is.na(properties$step_time))/40
  print(strat_by_rot_table)
  cat('\n(table not used, just confirms patterns)\n')
  
  par(mfrow=c(1,2))
  
  plot(y = NULL, x = NULL,
       ylab = 'rotation size / density',
       xlab = 'step time (trials)',
       xlim=c(0, 120), ylim=c(0.5,5.5),
       bty='n', axes=FALSE)
  
  
  propvals <- round(properties[which( properties$step_size > 5 & 
                                        properties$step_time >= 0 ), 'step_time'])
  all_gamma_fit <- MASS::fitdistr(propvals, densfun = "gamma")
  
  # cat('overall gamma distribution for step time:\n')
  # print(all_gamma_fit)
  
  rotation <- c(20,30,40,50,60)
  shape <- rep(all_gamma_fit$estimate['shape'], 5)
  rate  <- rep(all_gamma_fit$estimate['rate'], 5)
  
  step_time_gamma_par <- data.frame('rotation'=rotation, 'shape'=shape, 'rate'=rate)
  write.csv(step_time_gamma_par, file='data/step_time_gamma_parameters.csv', row.names=FALSE)
  
  # all_poisson_fit <- MASS::fitdistr(propvals, densfun = "poisson")
  # print(all_poisson_fit)
  
  rot_gamma_d <- c()
  
  X <- seq(.25, 120, length.out=480)
  gY <- dgamma(X, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
  # pY <- dpois(X, lambda=all_poisson_fit$estimate['lambda'])
  # print(pY)
  
  for (rot_idx in c(1,2,3,4,5)) {
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation
                                   & properties$step_size > 5
                                   & properties$step_time >= 0 
                                   # & properties$step_time >= 0 
                                   # & properties$step_time < 100
                                   & !is.na(properties$step_time)
                                 ), 'step_time']
    
    pvd <- density(propvals, na.rm=TRUE, bw=1.6,
                   n = 241, from=0, to=120)
    
    lines(pvd$x, (pvd$y/max(pvd$y))+rot_idx-0.5, col=rot_idx)
    points(propvals, rep(rot_idx-0.5, length(propvals)), col=rot_idx, pch=20, cex=0.5)
    
    # print(propvals)
    # gamma_fit <- MASS::fitdistr(propvals, densfun = "gamma", start=list(shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate']))
    gamma_fit <- MASS::fitdistr(propvals, densfun = "gamma")
    Y <- dgamma(X, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'])
    # print(Y)
    lines(X, (Y/max(Y))+rot_idx-0.5, col=rot_idx, lw=1, lty=2)
    
    lines(X, (gY/max(gY))+rot_idx-0.5, col='purple', lw=1, lty=2)
    # lines(X, (pY/max(pY))+rot_idx-0.5, col='orange', lw=1, lty=2)
    
    rot_gamma_d <- c(rot_gamma_d, dgamma(propvals, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'] ))
    
    # norm_fit <- MASS::fitdistr(propvals, densfun = "normal")
    # Y <- dnorm(X, mean=norm_fit$estimate['mean'], sd=norm_fit$estimate['sd'])
    # lines(X, (Y/max(Y))+rot_idx-0.5, col=rot_idx, lw=1, lty=2)
    # print(norm_fit)
    
    
  }
  
  axis(side=1, at=c(0,30,60,90,120))
  axis(side=2, at=c(1,2,3,4,5), labels=c(20,30,40,50,60))
  
  rot_gamma_nll <- Reach::nll(d = rot_gamma_d)
  rotgAIC <- Reach::AIC(logLik = -1*rot_gamma_nll, k=10, N=length(rot_gamma_d))
  cat('\nStep Time:\n')
  cat(sprintf('one gamma AIC: %0.1f, rot. spec. gamma AIC: %0.1f\n\n', stats::AIC(all_gamma_fit, k=2), rotgAIC))
  
  plot(y = NULL, x = NULL,
       ylab = 'rotation size / density',
       xlab = 'step size (deg)',
       xlim=c(-10, 70), ylim=c(0.5,5.5),
       bty='n', axes=FALSE)
  
  
  X <- seq(-10, 70, length.out=241)
  
  fixed <- data.frame('m'=c(0, NA), 's'=c(NA,NA), 'w'=c(NA,NA))
  
  allpar <- NA
  
  for (rot_idx in c(1,2,3,4,5)) {
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation 
                                   # & !is.na(properties$step_size)
                                   # & properties$step_size >= -10
                                   # & properties$step_size <=  rotation + 10
                                   
                                   # & properties$step_time < 100
                                   # & !is.na(properties$stratdev_onset)
                                 ), 'step_size']
    
    # print(propvals)
    
    pvd <- density(propvals, na.rm=TRUE, bw=1.6,
                   n = 161, from=-10, to=70)

    lines(pvd$x, 0.9*(pvd$y/max(pvd$y))+rot_idx-0.45, col=rot_idx)
    points(propvals, rep(rot_idx-0.5, length(propvals)), col=rot_idx, pch=20, cex=0.5)
    
    fixed$m[2] <- rotation/2
    # fixed$w <- unname(strat_by_rot_table[rot_idx,c(1,2)])
    
    fitpar <- Reach::multiModalFit(x=propvals, n=2, points=6, best=4, fixed=fixed)

    fitpar$rotation <- rotation

    if (is.data.frame(allpar)) {
      allpar <- rbind(allpar, fitpar)
    } else {
      allpar <- fitpar
    }
    
    Reach::multiModalModel(x=X, par=fitpar) -> yvals
    lines(X, 0.9*(yvals/(max(yvals)))+rot_idx-0.45, col=rot_idx, lw=2)
    points(x = c(rotation, rotation/2),
           rep(rot_idx+0.4,2),
           col=rot_idx, pch=6, cex=1)
    
  }
  
  print(allpar)
  
  write.csv(allpar, file='data/step_size_multimodal_parameters.csv', row.names=FALSE)
  
  axis(side=1, at=c(0,20,40,60))
  axis(side=2, at=c(1,2,3,4,5), labels=c(20,30,40,50,60))
  
}

plotStepSD <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  par(mfrow=c(1,3))
  
  steptrue <- c(FALSE,        TRUE,         TRUE         )
  depvars   <- c('prestep_sd', 'prestep_sd', 'poststep_sd')
  
  X <- seq(.25, 40, length.out=160)
  
  
  rotation <- c()
  makestep <- c()
  phase    <- c()
  shape    <- c()
  rate     <- c()
  
  for (situation in c(1,2,3)) {
    if (steptrue[situation]) {
      sitprop <- properties[which(!is.na(properties$step_time)),]
    } else {
      sitprop <- properties[which(is.na(properties$step_time)),]
    }
    depvar <- depvars[situation]
    print(depvar)
    
    
    
    plot(y = NULL, x = NULL,
         ylab = 'density (by rotation size)',
         xlab = sprintf('%s (deg)', depvar ),
         xlim=c(0, max(X)), ylim=c(0.5,5.5),
         bty='n', axes=FALSE)
    
    
    all_propval <- sitprop[, depvar]
    # print(range(all_propval))
    # print(range(all_propval, na.rm=TRUE))
    all_propval <- all_propval[which(!is.na(all_propval))]
    all_propval[which(all_propval <= 0)] <- .Machine$double.eps
    all_gamma_fit <- MASS::fitdistr(all_propval, densfun = "gamma", lower=c(1.001, 0.001), upper=c(1000,1000))
    agY <- dgamma(X, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
    
    allAIC <- stats::AIC(all_gamma_fit, k=2)
    
    # all_norm_fit <- MASS::fitdistr(all_propval, densfun = "normal")
    # anY <- dnorm(X, mean=all_norm_fit$estimate['mean'], sd=all_norm_fit$estimate['sd'])
    
    rot_d <- c()
    
    for (rot_idx in c(1,2,3,4,5)) {
      rot <- c(20,30,40,50,60)[rot_idx]
      propvals <- sitprop[which(sitprop$rotation == rot
                                
      ), depvar]
      # print(propvals)
      propvals <- propvals[which(!is.na(propvals))]
      propvals[which(propvals <= 0)] <- .Machine$double.eps
      
      points(propvals, rep(rot_idx-0.55, length(propvals)), col=rot_idx, pch=20, cex=0.5)
      
      pvd <- density(propvals, na.rm=TRUE, 
                     n = length(X), from=min(X), to=max(X))
      
      lines(x = c(0,max(X)),
            y = rep(0,2)+rot_idx-0.5,
            col='#ddd')
      lines(x   = pvd$x, 
            y   = 0.9*(pvd$y/max(pvd$y))+rot_idx-0.5, 
            col = rot_idx)
      
      lines(x   = X, 
            y   = 0.9*(agY/max(agY))+rot_idx-0.5, 
            col = 'purple', lw=1, lty=2)
      # lines(x   = X, 
      #       y   = 0.9*(anY/max(anY))+rot_idx-0.5, 
      #       col = 'orange', lw=1, lty=2)
      
      
      gamma_fit <- MASS::fitdistr(propvals, densfun = "gamma", lower=c(1.001, 0.001), upper=c(1000,1000))
      Y <- dgamma(X, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'])
      lines(X, 0.9*(Y/max(Y))+rot_idx-0.5, col=rot_idx, lw=1, lty=2)
      
      # norm_fit <- MASS::fitdistr(propvals, densfun = "normal")
      # Y <- dnorm(X, mean=norm_fit$estimate['mean'], sd=norm_fit$estimate['sd'])
      # lines(X, 0.9*(Y/max(Y))+rot_idx-0.5, col=rot_idx, lw=1, lty=3)
      
      rot_d <- c(rot_d, dgamma(propvals, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'] ))
      
      # gamma_nll <- nll(d = dgamma(propvals, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'] ))
      # agAIC <- Reach::AIC(logLik = -1*gamma_nll, k=2, N=length(propvals))
      # cat(sprintf('gamma AIC: %0.1f, all gamma AIC: %0.1f\n',stats::AIC(gamma_fit, k=2), agAIC))
      # print(gamma_fit$estimate)
      
      # cat(sprintf('gamma AIC: %0.1f, normal AIC: %0.1f\n',stats::AIC(gamma_fit, k=2), stats::AIC(norm_fit, k=2)))
      
      # norm_nll  <- nll(d = dnorm(propvals,  mean=all_norm_fit$estimate['mean'], sd=all_norm_fit$estimate['sd']))
      # anAIC <- Reach::AIC(logLik = -1*norm_nll,  k=2, N=length(propvals))
      # cat(sprintf('all gamma AIC: %0.1f, all normal AIC: %0.1f\n',agAIC, anAIC))
      
      rotation <- c(rotation, rot)
      makestep <- c(makestep, steptrue[situation])
      phase    <- c(phase, depvar)
      shape    <- c(shape, gamma_fit$estimate['shape'])
      rate     <- c(rate, gamma_fit$estimate['rate'])
      
    }
    
    rot_gamma_nll <- nll(d = rot_d)
    rotgAIC <- Reach::AIC(logLik = -1*rot_gamma_nll, k=10, N=length(rot_d))
    cat(sprintf('one gamma AIC: %0.1f, 5rot gamma AIC: %0.1f\n', allAIC, rotgAIC))
    
    axis(side=1, at=c(0,20,40))
    axis(side=2, at=c(1,2,3,4,5), labels=c(20,30,40,50,60))
    
  }
  
  step_SD_gamma_distr <- data.frame('rotation'=rotation, 'makestep'=makestep, 'phase'=phase, 'shape'=shape, 'rate'=rate)
  write.csv(step_SD_gamma_distr, file='data/step_SD_gamma_parameters.csv', row.names=FALSE)
  
}

# bi-modal final strategy fits -----

fitFinalStrategyDistributions <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }

  fixed <- data.frame('m'=c(0, NA), 's'=c(NA,NA), 'w'=c(NA,NA))
  
  colname <- 'final_strategy'
  
  allpar <- NA
  
  for (rot_idx in c(1,2,3,4,5)) {
    
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation), colname]
    
    propvals <- propvals[which(propvals > -10)] # remove negative outliers
    
    # distribution_fit <- MASS::fitdistr(propvals, densfun = "normal")
    print(sprintf('variable %s, rotation %d', colname, rotation))
    fixed$m[2] <- rotation/2
    fitpar <- Reach::multiModalFit(x=propvals, n=2, points=6, best=4, fixed=fixed)
    
    fitpar$rotation <- rotation
    
    if (is.data.frame(allpar)) {
      allpar <- rbind(allpar, fitpar)
    } else {
      allpar <- fitpar
    }
  }
  
  write.csv(allpar, file='data/final_strategy_multimodal_parameters.csv', row.names=FALSE)
  return(allpar)
}


# fitNormalDistributions <- function() {
#   
#   properties <- extractEmpiricalProperties()
#   
#   colnames <- c('stratdev_onset', 'strat_stable', 'devel_duration', 'predev_sd', 'devel_sd', 'stable_sd')
#   
#   fixed <- data.frame('m'=c(0, NA), 's'=c(0.75,NA), 'w'=c(NA,NA))
#   
#   for (colname in colnames) {
#     
#     valrange <- range(properties[, colname], na.rm=TRUE)
#     
#     for (rot_idx in c(1,2,3,4,5)) {
#       
#       rotation <- c(20,30,40,50,60)[rot_idx]
#       propvals <- properties[which(properties$rotation == rotation), colname]
#       
#       # distribution_fit <- MASS::fitdistr(propvals, densfun = "normal")
# 
#       print(sprintf('variable %s, rotation %d', colname, c(20,30,40,50,60)[rot_idx]))
#       print(Reach::multiModalFit(x=propvals, n=1, points=20, best=10))
#     }
#     
#   }
#   
# }

plotFinalStratModeWeights <- function() {
  
  read.csv('data/final_strategy_multimodal_parameters.csv') -> allpar
  
  plot(x=NULL,y=NULL,
       xlim=c(15,65), ylim=c(0,1),
       xlab='rotation', ylab='proportion observations',
       bty='n', axes=FALSE)
  
  for (rotation in unique(allpar$rotation)) {
    rotpars <- allpar[which(allpar$rotation == rotation), ]
    weights <- rotpars$w
    polygon( x=c(rotation-4, rotation-4, rotation+4, rotation+4),
             y=c(0, weights[2], weights[2], 0),
             col='blue', border='black')
    polygon( x=c(rotation-4, rotation-4, rotation+4, rotation+4),
             y=c(weights[2], 1, 1, weights[2]),
             col='white', border='black')
  }
  
  axis(side=1, at=c(20,30,40,50,60))
  axis(side=2, at=c(0,0.25,0.5,0.75,1))
  
}

plotFinalStratDistributions <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  read.csv('data/final_strategy_multimodal_parameters.csv') -> allpar
  
  # properties <- extractEmpiricalProperties()
  colname <- 'final_strategy'
  
  valrange <- c(-10, 70)
  plot(NA, 
       xlim=valrange, ylim=c(0.5,5.5),
       xlab='final strategy magnitude', ylab='density', 
       main='',axes=FALSE,bty='n')
  X <- seq(valrange[1], valrange[2], length.out=321)
  
  for (rot_idx in c(1,2,3,4,5)) {
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation), colname]
    propvals <- propvals[which(propvals > -10)] # remove negative outliers
    
    pvd <- density(propvals, na.rm=TRUE, bw=1.6,
                   n = length(X), from=min(X), to=max(X))
    lines(pvd$x, 0.8*(pvd$y/(max(pvd$y)))+rot_idx-0.4, col=rot_idx, lw=2, lty=3)
    points(propvals, rep(rot_idx-0.4, length(propvals)), col=rot_idx, pch=20, cex=0.5)
    
    
    rotpars <- allpar[which(allpar$rotation == rotation), ]
    rotpars$s[1] <- 1.5 # fix the standard deviation of the first mode to 0.75
    Reach::multiModalModel(x=X, par=rotpars) -> yvals
    lines(X, 0.8*(yvals/(max(yvals)))+rot_idx-0.4, col=rot_idx, lw=2)
    points(x = c(rotation, rotation/2), 
           rep(rot_idx+0.4,2), 
           col=rot_idx, pch=6, cex=1)
    
  }
  
  axis(side=1,at=c(0,20,30,40,50,60),labels=c(0,20,30,40,50,60))
  axis(side=2,at=c(1,2,3,4,5),labels=c(20,30,40,50,60))
    
}

# strategy development onset - alpha distribution -----

fitOnsetGammaDistributions <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  colname <- 'stratdev_onset'
  
  valrange <- range(properties[, colname], na.rm=TRUE)
  
  propvals <- properties[,colname]
  
  
  all_gamma_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "gamma")
  cat('\nOnset Trial ALL conditions:\n')
  # print(all_gamma_fit)
  
  normal_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "normal")
  poisson_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "poisson")
  
  cat(sprintf(' gamma AIC: %0.1f\n',stats::AIC(all_gamma_fit, k=2)))
  cat(sprintf(' normal AIC: %0.1f\n',stats::AIC(normal_fit, k=2)))
  cat(sprintf(' poisson AIC: %0.1f\n\n',stats::AIC(poisson_fit, k=2)))
  
  # dgamma <- dgamma(propvals, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
  # gamma_AIC_m <- Reach::AIC(logLik=-1*Reach::nll(dgamma),k=2,N=length(propvals))
  # cat(sprintf(' mgamma AIC: %0.1f\n',gamma_AIC_m))
  
  rot_gamma_d <- c()
  
  
  for (rot_idx in c(1,2,3,4,5)) {
    
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation), colname]

    gamma_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "gamma")
    # cat(sprintf('\nOnset Trial for %d deg:\n', rotation))
    # print(gamma_fit)
    
    normal_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "normal")
    # cat('normal:\n')
    # print(normal_fit)

    poisson_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "poisson")
    
    # cat(sprintf(' gamma AIC: %0.1f\n',stats::AIC(gamma_fit, k=10)))
    # cat(sprintf(' normal AIC: %0.1f\n',stats::AIC(normal_fit, k=10)))
    # cat(sprintf(' poisson AIC: %0.1f\n',stats::AIC(poisson_fit, k=10)))
    
    dgamma <- dgamma(propvals, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
    rot_gamma_d <- c(rot_gamma_d, dgamma)
    # gamma_AIC_m <- Reach::AIC(logLik=-1*Reach::nll(dgamma),k=2,N=length(propvals))
    # cat(sprintf(' all gamma AIC: %0.1f\n',gamma_AIC_m))
    
    
  }
  
  rot_gamma_AIC_m <- Reach::AIC(logLik=-1*Reach::nll(rot_gamma_d),k=10,N=length(rot_gamma_d))
  cat(sprintf(' rot. spec. gammas AIC: %0.1f\n',rot_gamma_AIC_m))
  
  
  cat('\n')

  rotation <- c(20,30,40,50,60)
  shape <- rep(all_gamma_fit$estimate['shape'], length(rotation)) 
  rate <- rep(all_gamma_fit$estimate['rate'], length(rotation))
  
  par_df <- data.frame('rotation'=rotation, 'shape'=shape, 'rate'=rate)
  write.csv(par_df, file='data/stratdev_onset_gamma_parameters.csv', row.names=FALSE)

}

plotOnsetGamma <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  colname <- 'stratdev_onset'
  
  propvals <- properties[,colname]
  gamma_fit_all <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "gamma")
  
  valrange <- range(properties[, colname], na.rm=TRUE)
  valrange[1] <- 0
  X <- seq(valrange[1], valrange[2], length.out=321)

  plot(NULL,NULL,
       xlim=valrange, ylim=c(0.5,5.5),
       xlab='strategy development onset trial', ylab='density', 
       main='',axes=FALSE,bty='n')
  
  for (rot_idx in c(1,2,3,4,5)) {
    
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation), colname]
    propvals <- propvals[which(propvals < 80)]
    
    gamma_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "gamma")
    
    pvd <- density(propvals, na.rm=TRUE, bw=1.6,
                   n = length(X), from=min(X), to=max(X))
    lines(pvd$x, 0.8*(pvd$y/(max(pvd$y)))+rot_idx-0.4, col=rot_idx, lw=2, lty=3)
    points(propvals, rep(rot_idx-0.4, length(propvals)), col=rot_idx, pch=20, cex=0.5)
    
    yvals <- dgamma(X, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'])
    lines(X, 0.8*(yvals/(max(yvals)))+rot_idx-0.4, col=rot_idx, lw=0.5)
    
    yvals_all <- dgamma(X, shape=gamma_fit_all$estimate['shape'], rate=gamma_fit_all$estimate['rate'])
    lines(X, 0.8*(yvals_all/(max(yvals_all)))+rot_idx-0.4, col='purple', lw=2, lty=2)
    
  }
  
  axis(side=2,at=c(1,2,3,4,5),labels=c(20,30,40,50,60))
  axis(side=1,at=pretty(valrange),labels=pretty(valrange))
  # dev.off()
  
}

# strategy development duration -----

checkStratDevDuration <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  properties[which(!is.na(properties$devel_duration)),] -> properties
  
  par(mar=c(5,5,5,5))
  plot(y = properties$stratdev_onset, 
       x = properties$devel_duration,
       xlim=c(0,120), ylim=c(0,120),
       ylab = 'strategy development onset trial',
       xlab = 'strategy development duration (trials)',
       pch=20, col=properties$rotation,
       asp=1, bty='n', axes=FALSE)
  lines(x = c(0, 120), y = c(120, 0), col='#999999', lw=1, lty=1)
  axis(side=1, at=pretty(c(0,120)), labels=pretty(c(0,120)))
  axis(side=2, at=pretty(c(0,120)), labels=pretty(c(0,120)))
  legend(90, 120, legend=c(20,30,40,50,60), col=c(1,2,3,4,5), pch=20, bty='n')
  
  onset_density <- density(properties$stratdev_onset[which(!is.na(properties$stratdev_onset))], na.rm=TRUE, from=0, to=120, n=241)
  devdur_density <- density(properties$devel_duration[which(!is.na(properties$devel_duration))], na.rm=TRUE, from=0, to=120, n=241)
  
  polygon(y = c(0,onset_density$x,120), 
          x = c(0,((onset_density$y/max(onset_density$y))*25),0)+125, 
          col='#00009966', border=NA, xpd=TRUE)
  polygon(y = c(0,((devdur_density$y/max(devdur_density$y))*15),0)+125, 
          x = c(0,devdur_density$x, 120), 
          col='#00009966', border=NA, xpd=TRUE)
  
  # is there a correlation with onset?
  print(sprintf('correlation onset vs duration:'))
  print(cor.test(properties$stratdev_onset, properties$devel_duration, na.rm=TRUE))
  
  # is it bimodal?
  devdur <- properties$devel_duration[which(!is.na(properties$devel_duration))]
  fitpar <- Reach::multiModalFit(x=devdur, n=2, points=6, best=4)
  print(sprintf('strategy development duration bimodal fit:'))
  print(fitpar)
  
  rotdurbimodpar <- NA
  for (rotation in c(20,30,40,50,60)) {
    fitpar$rotation <- rotation
    if (is.data.frame(rotdurbimodpar)) {
      rotdurbimodpar <- rbind(rotdurbimodpar, fitpar)
    } else {
      rotdurbimodpar <- fitpar
    }
  }
  # print(rotdurbimodpar)
  write.csv(rotdurbimodpar, file='data/stratdev_duration_multimodal_parameters.csv', row.names=FALSE)
  
  X <- devdur_density$x
  
  Reach::multiModalModel(x=X, par=fitpar) -> yvals
  lines(x=X, y=(18*(yvals/(max(yvals))))+125, col='red',   lw=1, lty=2, xpd=TRUE)
  
  gamma_fit <- MASS::fitdistr(devdur, densfun = "gamma")
  yvals <- dgamma(X, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'])
  lines(x=X, y=(18*(yvals/(max(yvals))))+125, col='green', lw=1, lty=2, xpd=TRUE)
  
  gamma_d <- dgamma(devdur, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'])
  gamma_nll <- Reach::nll(gamma_d)

  norm_d <- (fitpar$w[1] * dnorm(devdur, mean=fitpar$m[1], sd=fitpar$s[1])) + (fitpar$w[2] * dnorm(devdur, mean=fitpar$m[2], sd=fitpar$s[2]))
  norm_nll <- Reach::nll(norm_d)
  
  # print(log(gamma_d))
  # print(log(norm_d))
  print(sprintf('gamma NLL: %0.2f', gamma_nll))
  print(sprintf('bimodal NLL: %0.2f', norm_nll))
  
  logLik <- c('gamma'=gamma_nll, 'bimodal'=norm_nll) * -1
  print(Reach::AIC(logLik=logLik, k=c(2,5), N=length(devdur)))  
}

# standard deviations -----

correlateSDs <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  # one of them is too high:
  properties <- properties[which(properties$stable_sd < 40),]
  
  sd_vars <- c('predev_sd','devel_sd','stable_sd')
  # cat('outer loop:\n')
  # print(c(1,length(sd_vars)-1))
  for (first in c(1,length(sd_vars)-1)) {
    # cat('  inner loop:\n  ')
    # print(unique(c(first+1,length(sd_vars))))
    for (second in unique(c(first+1,length(sd_vars)))) {
      
      cat(sprintf('%s ~ %s:\n', sd_vars[second], sd_vars[first]))
      x = properties[,sd_vars[first]]
      y = properties[,sd_vars[second]]
      print(cor.test(x, y))
      
    }
  }
  # cat('done\n')
  # return()
  
}

plotSDbyRotation <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- getProperties()
  }
  
  # one of them is too high:
  properties <- properties[which(properties$stable_sd < 40),]
  
  layout(mat=matrix(c(1,2,3),nrow=1,byrow=TRUE))
  
  phase = c()
  rotation = c()
  shape = c()
  rate = c()
  
  sd_vars <- c('predev_sd','devel_sd','stable_sd')
  for (col_idx in c(1:length(sd_vars))) {
    colname <- sd_vars[col_idx]
    valrange <- c(0, c(90,60,15)[col_idx])
    
    #  valrange <- c(-10, 70)
    plot(NA, 
         xlim=valrange, ylim=c(0.5,5.5),
         xlab=colname, ylab='density by rotation size', 
         main='',axes=FALSE,bty='n')
    X <- seq(valrange[1], valrange[2], length.out=321)
    
    propvals <- properties[, colname]
    propvals <- propvals[which(!is.na(propvals))]
    all_gamma_fit <- MASS::fitdistr(propvals, densfun = "gamma")
    cat(sprintf('\n -= %s =-\n\n', toupper(colname)))
    cat(sprintf('one gamma AIC: %0.1f\n',stats::AIC(all_gamma_fit, k=2)))
    # print(all_gamma_fit$estimate)
    
    all_norm_fit <- MASS::fitdistr(propvals, densfun = "normal")
    cat(sprintf('one normal AIC: %0.1f\n',stats::AIC(all_norm_fit, k=2)))
    # print(all_norm_fit$estimate)
    
    rot_gamma_d <- c()
    
    phase <- c(phase, colname)
    rotation <- c(rotation, 200)
    shape <- c(shape, all_gamma_fit$estimate['shape'])
    rate <- c(rate, all_gamma_fit$estimate['rate'])

    for (rot_idx in c(1,2,3,4,5)) {
      rot <- c(20,30,40,50,60)[rot_idx]
      propvals <- properties[which(properties$rotation == rot), colname]
      propvals <- propvals[which(propvals > -10)] # remove negative outliers
      
      pvd <- density(propvals, na.rm=TRUE, bw=1.6,
                     n = length(X), from=min(X), to=max(X))
      lines(pvd$x, 0.8*(pvd$y/(max(pvd$y)))+rot_idx-0.4, col=rot_idx, lw=1, lty=1)
      points(propvals, rep(rot_idx-0.5, length(propvals)), col=rot_idx, pch=20, cex=0.5)
      
      propvals <- propvals[which(!is.na(propvals))]
      gamma_fit <- MASS::fitdistr(propvals, densfun = "gamma")
      normal_fit <- MASS::fitdistr(propvals, densfun = "normal")
      # poisson_fit <- MASS::fitdistr(propvals, densfun = "poisson")
      
      phase <- c(phase, colname)
      rotation <- c(rotation, rot)
      shape <- c(shape, gamma_fit$estimate['shape'])
      rate <- c(rate, gamma_fit$estimate['rate'])
      
      # dgamma <- dgamma(propvals, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
      # all_gamma_AIC <- Reach::AIC(logLik=-1*Reach::nll(dgamma),k=2,N=length(propvals))[1]
      # # print(all_gamma_AIC)
      # 
      # cat(sprintf('%s, %d -- gamma AIC: %0.1f,  normal AIC: %0.1f\n', colname, rot ,stats::AIC(gamma_fit, k=10),stats::AIC(normal_fit, k=10)))
      # cat(sprintf("all gamma AIC: %0.1f\n", all_gamma_AIC))
      
      rot_gamma_d <- dgamma(propvals, shape=gamma_fit$estimate['shape'], rate=gamma_fit$estimate['rate'])

    }
    
    axis(side=1,at=seq(0,max(valrange),max(valrange)/3))
    axis(side=2,at=c(1,2,3,4,5),labels=c(20,30,40,50,60))
    
    rot_gamma_AIC_m <- Reach::AIC(logLik=-1*Reach::nll(rot_gamma_d),k=10,N=length(rot_gamma_d))
    cat(sprintf('rot. spec. gammas AIC: %0.1f\n', rot_gamma_AIC_m))
  
  }
  
  gamma_fit_par <- data.frame(phase, rotation, shape, rate)
  
  write.csv(gamma_fit_par, 'data/SD_gamma_par.csv', row.names = FALSE)
  
  # return(gamma_fit_par)
  
}


plotSDfitPars <- function() {
  
  df <- read.csv('data/SD_gamma_par.csv', stringsAsFactors = FALSE)
  
  layout(mat=matrix(c(1:6),ncol=2,byrow=TRUE))
  
  for (phase in unique(df$phase)) {
    
    for (gamma_par in c('shape','rate')) {
      
      subdf <- df[which(df$phase == phase & df$rotation < 200),]
      
      plot(x=subdf$rotation, y=subdf[,gamma_par],
           xlim=c(0,65),ylim=c(0,max(subdf[,gamma_par])*1.05),
           main=sprintf('%s', phase),
           xlab='rotation', ylab=gamma_par,
           bty='n')
      
    }
  }
  
}