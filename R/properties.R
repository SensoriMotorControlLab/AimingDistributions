
# generic / all properties ----

extractEmpiricalProperties <- function() {
  
  
  ppno <- 0
  
  participant    <- c()
  rotation       <- c()
  final_strategy <- c()
  stratdev_onset <- c()
  strat_stable   <- c()
  devel_duration <- c()
  predev_sd      <- c()
  devel_sd       <- c()
  stable_sd      <- c()
  
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
      
      final_strat    <- median(ARtimecourse[c((length(ARtimecourse)-ntrials+1):length(ARtimecourse))])
      final_strategy <- c(final_strategy, final_strat)
      
      onset <- which(ARtimecourse > 5)[1]
      if (final_strat < 5) {
        onset <- NA
      }
      stratdev_onset <- c(stratdev_onset, onset)
      
      if (is.na(onset)) {
        stab_trial <- NA
      } else {
        stab_tc <- ARtimecourse[onset:length(ARtimecourse)]
        stab_trial <- Reach::findStabilizationTrial(stab_tc)+onset
      }
      strat_stable <- c(strat_stable, stab_trial)
      
      devel_duration <- c(devel_duration, stab_trial - onset)
      
      if (is.na(onset)) {
        predev_sd <- c(predev_sd, sd(ARtimecourse))
        devel_sd  <- c(devel_sd, NA)
        stable_sd <- c(stable_sd, NA)
      } else if (is.na(stab_trial)) {
        predev_sd      <- c(predev_sd, sd(ARtimecourse[1:onset]))
        devel_sd       <- c(devel_sd, NA)
        stable_sd      <- c(stable_sd, NA)
      } else {
        predev_sd      <- c(predev_sd, sd(ARtimecourse[1:onset]))
        devel_sd       <- c(devel_sd, sd(ARtimecourse[onset:stab_trial]))
        stable_sd      <- c(stable_sd, sd(ARtimecourse[stab_trial:length(ARtimecourse)]))
      }

    }
   
  }
  
  return( data.frame (
                        participant    = participant,
                        rotation       = rotation,
                        final_strategy = final_strategy,
                        stratdev_onset = stratdev_onset,
                        strat_stable   = strat_stable,
                        devel_duration = devel_duration,
                        predev_sd      = predev_sd,
                        devel_sd       = devel_sd,
                        stable_sd      = stable_sd
  ) )
  
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
  
  colnames <- c('final_strategy', 'stratdev_onset', 'strat_stable', 'devel_duration', 'predev_sd', 'devel_sd', 'stable_sd')
  
  par(mar=c(4,3,.2,.2))
  layout(matrix(1:(2*(ceiling(length(colnames)/2))), nrow=ceiling(length(colnames)/2), ncol=2, byrow=TRUE))
  
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

# bi-modal final strategy fits -----

fitFinalStrategyDistributions <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- extractEmpiricalProperties()
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
    properties <- extractEmpiricalProperties()
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
    properties <- extractEmpiricalProperties()
  }
  
  colname <- 'stratdev_onset'
  
  valrange <- range(properties[, colname], na.rm=TRUE)
  
  propvals <- properties[,colname]
  
  
  all_gamma_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "gamma")
  cat('\nOnset Trial ALL conditions:\n')
  # print(gamma_fit)
  
  normal_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "normal")
  poisson_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "poisson")
  
  cat(sprintf(' gamma AIC: %0.1f\n',stats::AIC(all_gamma_fit, k=2)))
  cat(sprintf(' normal AIC: %0.1f\n',stats::AIC(normal_fit, k=2)))
  cat(sprintf(' poisson AIC: %0.1f\n\n',stats::AIC(poisson_fit, k=2)))
  
  dgamma <- dgamma(propvals, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
  gamma_AIC_m <- Reach::AIC(logLik=-1*Reach::nll(dgamma),k=2,N=length(propvals))
  cat(sprintf(' mgamma AIC: %0.1f\n',gamma_AIC_m))
  
  for (rot_idx in c(1,2,3,4,5)) {
    
    rotation <- c(20,30,40,50,60)[rot_idx]
    propvals <- properties[which(properties$rotation == rotation), colname]

    gamma_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "gamma")
    cat(sprintf('\nOnset Trial for %d deg:\n', rotation))
    # print(gamma_fit)
    
    normal_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "normal")
    # cat('normal:\n')
    # print(normal_fit)

    poisson_fit <- MASS::fitdistr(propvals[which(!is.na(propvals))], densfun = "poisson")
    
    cat(sprintf(' gamma AIC: %0.1f\n',stats::AIC(gamma_fit, k=10)))
    cat(sprintf(' normal AIC: %0.1f\n',stats::AIC(normal_fit, k=10)))
    cat(sprintf(' poisson AIC: %0.1f\n',stats::AIC(poisson_fit, k=10)))
    
    dgamma <- dgamma(propvals, shape=all_gamma_fit$estimate['shape'], rate=all_gamma_fit$estimate['rate'])
    gamma_AIC_m <- Reach::AIC(logLik=-1*Reach::nll(dgamma),k=2,N=length(propvals))
    cat(sprintf(' all gamma AIC: %0.1f\n',gamma_AIC_m))
    
    
  }
  
  cat('\n')

  return(all_gamma_fit)
  
}

plotOnsetGamma <- function(properties=NULL) {
  
  if (is.null(properties)) {
    properties <- extractEmpiricalProperties()
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
    properties <- extractEmpiricalProperties()
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



