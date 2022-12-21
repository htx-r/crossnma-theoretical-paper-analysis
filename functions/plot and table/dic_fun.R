# CHECK:  others (BUGSnet and Phillipo) extract resdev to compute "dev" not dev as I did. 
# But resdev is computed per study not per arm, so that the pmdev while the pmdev_fitted computed per arm
# how will we subtract the pmdev_fitted from pmdev if they don't have the same length (when we compute the leverage)
# computing r and n doesnt work that way, they are matrices which could have NA's when we have different number of arms.
dic.fun <- function(x){
  # To compute DIC, we need to calculate 1. D_bar_model 2. D_at_rtilde 
  jagssamples <- as.mcmc(x) # I NEED R2jags and coda packages to run that
  
  samples = do.call(rbind, x$samples)%>% data.frame()
  
  #** For IPD
  if(x$model$data$ns.ipd!=0){
  # 1. calculate D_bar_model (the posterior mean of residual deviance)
  totresdev.ipd <- samples$totresdev.ipd %>% mean()
  
  # 2. calculate D_at_rtilde
  dev.ipd <- samples %>% select(., starts_with("dev.ipd")) # the samples of deviance 
  pmdev.ipd <- colMeans(dev.ipd) # the posterior mean deviance at each study
  rhat.ipd <- samples %>% select(., starts_with("rhat.ipd")) # samples of the expected event at each study arm
  rtilde.ipd <- rhat.ipd %>%colMeans()
  n.ipd00 <- x$model$data$n.ipd
  n.ipd0 <- n.ipd00[!is.na(n.ipd00)]
  n.ipd <- rep(n.ipd0, times=n.ipd0)     # get the sample size
  
  r.ipd00 <- x$model$data$r.ipd
  r.ipd0 <- r.ipd00[!is.na(r.ipd00)]
  r.ipd <- rep(r.ipd0, times=r.ipd0)
  
  
  pmdev_fitted.ipd <- -2*(ifelse(rtilde.ipd>.Machine$double.eps,r.ipd*log(rtilde.ipd/n.ipd),0)+ifelse(rtilde.ipd-1>.Machine$double.eps,(n.ipd-r.ipd)*log((1-(rtilde.ipd/n.ipd))),0))
  #2 * (ifelse(pasi75_r > .Machine$double.eps, pasi75_r * log(pasi75_r / fitted), 0) +
   #      ifelse(pasi75_n - pasi75_r > .Machine$double.eps, (pasi75_n - pasi75_r) * log((pasi75_n - pasi75_r) / (pasi75_n - fitted)), 0))
  # FINAL: claculate DIC
  leverage.ipd = pmdev.ipd-pmdev_fitted.ipd
  DIC.ipd = totresdev.ipd + sum(leverage.ipd,na.rm = T)
  DIC.ipd
  } else{
    DIC.ipd <- 0
    totresdev.ipd <- 0
  }
  
  if(x$model$data$ns.ad!=0){
  #** For AD
  # 1. calculate D_bar_model (the posterior mean of residual deviance)
  totresdev.ad <- samples$totresdev.ad %>% mean()
  
  # 2. calculate D_at_rtilde
  dev.ad <- samples %>% select(., starts_with("resdev.ad")) # the samples of deviance 
  pmdev.ad <- colMeans(dev.ad) # the posterior mean deviance at each study
  rhat.ad <- samples %>% select(., starts_with("rhat.ad")) # samples of the expected event at each study arm
  rtilde.ad <- rhat.ad %>%colMeans()
  r.ad <- c(x$model$data$r)     # get the observed number of events
  r.ad <- r.ad[!is.na(r.ad)]
  n.ad <- c(x$model$data$n)     # get the sample size
  n.ad <- n.ad[!is.na(n.ad)]
  pmdev_fitted.ad <- 2*(r.ad*log(r.ad/rtilde.ad)+(n.ad-r.ad)*log((n.ad-r.ad)/(n.ad-rtilde.ad)))
  
  # FINAL: claculate DIC
  leverage.ad = pmdev.ad-pmdev_fitted.ad
  DIC.ad = totresdev.ad + sum(leverage.ad,na.rm = T)
  DIC.ad
  } else{
    DIC.ad <- 0
  }
  cat("The Deviance Information Criterion (DIC) for IPDs is", DIC.ipd,
      "\n The Deviance Information Criterion (DIC) for ADs is", DIC.ad,
      "\n The residual deviance for IPDs is",totresdev.ipd,
      "\n The residual deviance for ADs is",totresdev.ad)
}

