model_with_dic <- "model {
  #*** IPD part
  for (i in 1:np) {  # loop through individuals
    y[i]~dbern(p[i]) # binomial likelihood
    logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]] # logistic transformation - to estimate Odds Ratio (OR)
    # to compute the IPD deviance
    for (k in 1:na.ipd[study[i]]) {
      rhat.ipd[i,k] <- p[i]*n.ipd[study[i],k]
      dev.ipd[i,k] <- -2*(r.ipd[study[i],k]*log(p[i])+(n.ipd[study[i],k]-r.ipd[study[i],k])*log(1-p[i]))
    }
    resdev.ipd[i] <- sum(dev.ipd[i,1:na.ipd[study[i]]]) # IPD residual deviance for each indvidual (data point)
    }
  totresdev.ipd <- sum(resdev.ipd[]) # total IPD residual deviance 
  for(j in 1:(ns.ipd)){      # loop through IPD studies
    w[j,1]<- 0              # multi-arm correction is zero for reference arm
    theta[j,t.ipd[j,1]]<- 0 # treatment effect is zero for reference arm
    
    for (k in 2:na.ipd[j]){ # loop through non-referent IPD arms
      
      # Synthesize relative treatment effects
      
      theta[j,t.ipd[j,k]] ~ dnorm(md[j,t.ipd[j,k]],precd[j,t.ipd[j,k]])
      # multi-arm correction
      md[j,t.ipd[j,k]]<- mean[j,k] + sw[j,k]
      w[j,k]<- (theta[j,t.ipd[j,k]]  - mean[j,k])
      sw[j,k]<- sum(w[j,1:(k-1)])/(k-1)
      precd[j,t.ipd[j,k]]<- prec *2*(k-1)/k
      
      # consistency equation
      mean[j,k] <-d[t.ipd[j,k]] - d[t.ipd[j,1]]
      
      
    }
  }
  
  #*** AD part
  for (j in 1:ns.ad) {            # loop through AD studies
    w.ad[j,1]<- 0                 # multi-arm correction is zero for referent arm
    theta[j+ns.ipd,t.ad[j,1]]<- 0 # treatment effect is zero for referent arm
    
    for (k in 1:na.ad[j]) {                 # loop through AD arms
      r[j,k] ~ dbin(pa[j,t.ad[j,k]],n[j,k]) # binomial likelihood of number of events
      ## to compute the AD deviance
      rhat.ad[j,k] <- pa[j,t.ad[j,k]]*n[j,k]
      dev.ad[j,k] <- 2*(r[j,k] * (log(r[j,k]-log(rhat.ad[j,k]))) + (n[j,k]-r[j,k])*(log(n[j,k]-r[j,k])-log(n[j,k]-rhat.ad[j,k])))
    }
    resdev.ad[j] <- sum(dev.ad[j,1:na.ad[j]]) # AD residual deviance for each study (data point)
    logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm
    
    for (k in 2:na.ad[j]){                  # loop through non-referent AD arms
      logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]]) # logistic transformation Log Odds Ratio
      
      # Synthesize relative treatment effects
      theta[j+ns.ipd,t.ad[j,k]] ~ dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
      # multi-arm correction
      md.ad[j,t.ad[j,k]]<- mean.ad[j,k] + sw.ad[j,k]
      w.ad[j,k]<- (theta[j+ns.ipd,t.ad[j,k]]  - mean.ad[j,k])
      sw.ad[j,k]<- sum(w.ad[j,1:(k-1)])/(k-1)
      precd.ad[j,t.ad[j,k]]<- prec *2*(k-1)/k
      # consistency equations
      mean.ad[j,k] <-d[t.ad[j,k]] - d[t.ad[j,1]]
    }
  }
  totresdev.ad <- sum(resdev.ad[]) # total AD residual deviance 
  
  totresdev <- totresdev.ipd + totresdev.ad # IPD and AD final residual deviance 
  #*** PRIORS
  for (j in 1:(ns.ipd+ns.ad)) {u[j] ~ dnorm(0,.01)} # log-odds in referent arm
  # heterogeneity between theta's
  tau ~dunif(0,2)
  prec<- pow(tau,-2)
  d[1] <- 0 # treatment effect is zero for reference treatment
  for(k in 2:nt) {d[k] ~ dnorm(0,.01)}
  # Compute OR for each comparison
  for(i in 1:(nt-1)) {
    for (j in (i+1):nt) {
      OR[j,i]<- exp(d[j] - d[i])
      LOR[j,i]<- d[j] - d[i]}
  }
  
 

}

"
