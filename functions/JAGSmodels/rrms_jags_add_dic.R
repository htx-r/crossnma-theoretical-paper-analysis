rrms_jags_unadj_dic <- 
"
model {
  #*** IPD part
  for (i in 1:np) {  # loop through individuals
    y[i]~dbern(p[i]) # binomial likelihood
    logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]] # logistic transformation - to estimate Odds Ratio (OR)
    # calculate residual deviance per indvidual
    #for (k in 1:na.ipd[study[i]]) {
      #rhat.ipd[i] <- p[i]*n.ipd[study[i],k]
      dev.ipd[i] <- -2*(stnd[study[i],arm[i]]+r.ipd[study[i],arm[i]]*log(p[i])+(n.ipd[study[i],arm[i]]-r.ipd[study[i],arm[i]])*log(1-p[i])) 
    #}
    #resdev.ipd[i] <- sum(dev.ipd[i,1:na.ipd[study[i]]]) 
  }
  totresdev.ipd <- sum(dev.ipd[]) # total IPD residual deviance 

  for(j in 1:(ns.ipd)){      # loop through IPD studies
     w[j,1]<- 0              # multi-arm correction is zero for reference arm
     theta[j,t.ipd[j,1]]<- 0 # treatment effect is zero for reference arm
     
     for (k in 2:na.ipd[j]){ # loop through non-referent IPD arms

      # Synthesize relative treatment effects
      
    theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]
    md[j,t.ipd[j,k]]<- mean[j,k]

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
      dev.ad[j,k] <- 2*(r[j,k] * (log(r[j,k])-log(rhat.ad[j,k])) + (n[j,k]-r[j,k])*(log(n[j,k]-r[j,k])-log(n[j,k]-rhat.ad[j,k])))
    }
    resdev.ad[j] <- sum(dev.ad[j,1:na.ad[j]]) # AD residual deviance for each study (data point)
    logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm

    for (k in 2:na.ad[j]){                  # loop through non-referent AD arms
      logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]]) # logistic transformation Log Odds Ratio

      # Synthesize relative treatment effects
      
    theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]
    md.ad[j,t.ad[j,k]]<- mean.ad[j,k]
      # consistency equations
      mean.ad[j,k] <-d[t.ad[j,k]] - d[t.ad[j,1]]
      
    }
  }
totresdev.ad <- sum(resdev.ad[]) # total AD residual deviance 

  #*** PRIORS
  for (j in 1:(ns.ipd+ns.ad)) {u[j] ~ dnorm(0,.01)} # log-odds in referent arm
  
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
rrms_jags_prior_dic <- "model {
  #*** IPD part
  for (i in 1:np) {  # loop through individuals
    y[i]~dbern(p[i]) # binomial likelihood
    logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]] # logistic transformation - to estimate Odds Ratio (OR)
  # calculate residual deviance per indvidual
    #for (k in 1:na.ipd[study[i]]) {
      #rhat.ipd[i] <- p[i]*n.ipd[study[i],k]
      dev.ipd[i] <- -2*(stnd[study[i],arm[i]]+r.ipd[study[i],arm[i]]*log(p[i])+(n.ipd[study[i],arm[i]]-r.ipd[study[i],arm[i]])*log(1-p[i])) 
    #}
    #resdev.ipd[i] <- sum(dev.ipd[i,1:na.ipd[study[i]]]) 
  }
  totresdev.ipd <- sum(dev.ipd[]) # total IPD residual deviance 

  for(j in 1:(ns.ipd)){      # loop through IPD studies
     w[j,1]<- 0              # multi-arm correction is zero for reference arm
     theta[j,t.ipd[j,1]]<- 0 # treatment effect is zero for reference arm
     
     for (k in 2:na.ipd[j]){ # loop through non-referent IPD arms

      # Synthesize relative treatment effects
      
    theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]
    md[j,t.ipd[j,k]]<- mean[j,k]

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
      dev.ad[j,k] <- 2*(r[j,k] * (log(r[j,k])-log(rhat.ad[j,k])) + (n[j,k]-r[j,k])*(log(n[j,k]-r[j,k])-log(n[j,k]-rhat.ad[j,k])))
    }
    resdev.ad[j] <- sum(dev.ad[j,1:na.ad[j]]) # AD residual deviance for each study (data point)
    logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm

    for (k in 2:na.ad[j]){                  # loop through non-referent AD arms
      logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]]) # logistic transformation Log Odds Ratio

      # Synthesize relative treatment effects
      
    theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]
    md.ad[j,t.ad[j,k]]<- mean.ad[j,k]
      # consistency equations
      mean.ad[j,k] <-d[t.ad[j,k]] - d[t.ad[j,1]]
      
    }
  }
totresdev.ad <- sum(resdev.ad[]) # total AD residual deviance 

  #*** PRIORS
  for (j in 1:(ns.ipd+ns.ad)) {u[j] ~ dnorm(0,.01)} # log-odds in referent arm
  
  d[1] <- 0 # treatment effect is zero for reference treatment
  d[2]~dnorm(0,0.01)
                    d[3]~dnorm(-0.00157505304230088,4.67727342972019)
                    d[4]~dnorm(0,0.01)
                    
  # Compute OR for each comparison
  for(i in 1:(nt-1)) {
    for (j in (i+1):nt) {
      OR[j,i]<- exp(d[j] - d[i])
      LOR[j,i]<- d[j] - d[i]}
      }
  
  }"
rrms_jags_adj1_dic <- 
  "
model {
  #*** IPD part
  for (i in 1:np) {  # loop through individuals
    y[i]~dbern(p[i]) # binomial likelihood
    logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]]+R[study[i]]*gamma[study[i]] # logistic transformation - to estimate Odds Ratio (OR)
  # calculate residual deviance per indvidual
    #for (k in 1:na.ipd[study[i]]) {
      #rhat.ipd[i] <- p[i]*n.ipd[study[i],k]
      dev.ipd[i] <- -2*(stnd[study[i],arm[i]]+r.ipd[study[i],arm[i]]*log(p[i])+(n.ipd[study[i],arm[i]]-r.ipd[study[i],arm[i]])*log(1-p[i])) 
    #}
    #resdev.ipd[i] <- sum(dev.ipd[i,1:na.ipd[study[i]]]) 
  }
  totresdev.ipd <- sum(dev.ipd[]) # total IPD residual deviance 

  for(j in 1:(ns.ipd)){      # loop through IPD studies
     w[j,1]<- 0              # multi-arm correction is zero for reference arm
     theta[j,t.ipd[j,1]]<- 0 # treatment effect is zero for reference arm
     
     for (k in 2:na.ipd[j]){ # loop through non-referent IPD arms

      # Synthesize relative treatment effects
      
    theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]
    md[j,t.ipd[j,k]]<- mean[j,k]

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
      dev.ad[j,k] <- 2*(r[j,k] * (log(r[j,k])-log(rhat.ad[j,k])) + (n[j,k]-r[j,k])*(log(n[j,k]-r[j,k])-log(n[j,k]-rhat.ad[j,k])))
    }
    resdev.ad[j] <- sum(dev.ad[j,1:na.ad[j]]) # AD residual deviance for each study (data point)
    logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm

    for (k in 2:na.ad[j]){                  # loop through non-referent AD arms
      logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]])+R[j+ns.ipd]*gamma[(j+ns.ipd)] # logistic transformation Log Odds Ratio

      # Synthesize relative treatment effects
      
    theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]
    md.ad[j,t.ad[j,k]]<- mean.ad[j,k]
      # consistency equations
      mean.ad[j,k] <-d[t.ad[j,k]] - d[t.ad[j,1]]
      
    }
  }
totresdev.ad <- sum(resdev.ad[]) # total AD residual deviance 

  #*** PRIORS
  for (j in 1:(ns.ipd+ns.ad)) {u[j] ~ dnorm(0,.01)} # log-odds in referent arm
  
  d[1] <- 0 # treatment effect is zero for reference treatment
  for(k in 2:nt) {d[k] ~ dnorm(0,.01)}
  # Compute OR for each comparison
  for(i in 1:(nt-1)) {
    for (j in (i+1):nt) {
      OR[j,i]<- exp(d[j] - d[i])
      LOR[j,i]<- d[j] - d[i]}
      }
  
  
  
  
  # Common effect for gamma (bias effect)
for (j in std.in) {gamma[j]<-g}
g~dnorm(0, 0.01)
prec.gamma <- 0
                      # bias adjustment
                      for (j in 1:(ns.ipd+ns.ad)) {R[j]~dbern(pi[bias_index[j]])}
                      pi[1]~dbeta(100,1)
                      pi[2]~dbeta(1,100)
                      pi[3]~dbeta(30,1)
                      pi[4]~dbeta(1,30)
                      pi[5]~dbeta(1,1)
  
  }
"
rrms_jags_adj2_dic <- 
"model {
  #*** IPD part
  for (i in 1:np) {  # loop through individuals
    y[i]~dbern(p[i]) # binomial likelihood
    logit(p[i]) <- u[study[i]]+theta[study[i],trt[i]] # logistic transformation - to estimate Odds Ratio (OR)
 # calculate residual deviance per indvidual
    #for (k in 1:na.ipd[study[i]]) {
      #rhat.ipd[i] <- p[i]*n.ipd[study[i],k]
      dev.ipd[i] <- -2*(stnd[study[i],arm[i]]+r.ipd[study[i],arm[i]]*log(p[i])+(n.ipd[study[i],arm[i]]-r.ipd[study[i],arm[i]])*log(1-p[i])) 
    #}
    #resdev.ipd[i] <- sum(dev.ipd[i,1:na.ipd[study[i]]]) 
  }
  totresdev.ipd <- sum(dev.ipd[]) # total IPD residual deviance 

  for(j in 1:(ns.ipd)){      # loop through IPD studies
     w[j,1]<- 0              # multi-arm correction is zero for reference arm
     theta[j,t.ipd[j,1]]<- 0 # treatment effect is zero for reference arm
     
     for (k in 2:na.ipd[j]){ # loop through non-referent IPD arms

      # Synthesize relative treatment effects
      theta[j,t.ipd[j,k]] <- md[j,t.ipd[j,k]]+(pi[bias_index[j]]*gamma[j])
        md[j,t.ipd[j,k]]<- mean[j,k]

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
      dev.ad[j,k] <- 2*(r[j,k] * (log(r[j,k])-log(rhat.ad[j,k])) + (n[j,k]-r[j,k])*(log(n[j,k]-r[j,k])-log(n[j,k]-rhat.ad[j,k])))
    }
    resdev.ad[j] <- sum(dev.ad[j,1:na.ad[j]]) # AD residual deviance for each study (data point)

    logit(pa[j,t.ad[j,1]]) <- u[j+ns.ipd]   # Log odds at referent arm

    for (k in 2:na.ad[j]){                  # loop through non-referent AD arms
      logit(pa[j,t.ad[j,k]]) <- u[j+ns.ipd]+(theta[j+ns.ipd,t.ad[j,k]]) # logistic transformation Log Odds Ratio

      # Synthesize relative treatment effects
      theta[j+ns.ipd,t.ad[j,k]] <- md.ad[j,t.ad[j,k]]+(pi[bias_index[j]]*gamma[j+ns.ipd])
        md.ad[j,t.ad[j,k]]<- mean.ad[j,k]
      # consistency equations
      mean.ad[j,k] <-d[t.ad[j,k]] - d[t.ad[j,1]]
    }
  }
totresdev.ad <- sum(resdev.ad[]) # total AD residual deviance 

  #*** PRIORS
  for (j in 1:(ns.ipd+ns.ad)) {u[j] ~ dnorm(0,.01)} # log-odds in referent arm
   
  d[1] <- 0 # treatment effect is zero for reference treatment
  for(k in 2:nt) {d[k] ~ dnorm(0,.01)}
  # Compute OR for each comparison
  for(i in 1:(nt-1)) {
    for (j in (i+1):nt) {
      OR[j,i]<- exp(d[j] - d[i])
      LOR[j,i]<- d[j] - d[i]}
      }
  
  
  
  
  # Common effect for gamma (bias effect)
for (j in std.in) {gamma[j]<-g}
g~dnorm(0, 0.01)
prec.gamma <- 0
                      # bias adjustment
                      pi[1]~dbeta(100,1)
                      pi[2]~dbeta(1,100)
                      pi[3]~dbeta(30,1)
                      pi[4]~dbeta(1,30)
                      pi[5]~ dbeta(1,1)
  
  }"