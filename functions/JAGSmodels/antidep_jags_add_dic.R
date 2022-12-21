## 22.08.2022 I created JAGS model using the crossnma.model(), 
# then I copy/paste it here and add the part to calculate the deviance 

antidep_jags_unadj_dic <- "
model {
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

antidep_jags_adj1_dic <- "
model {
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
  
  
  
  
  # Random effect for gamma (bias effect)
for (j in std.in) {gamma[j]~dnorm(g,prec.gamma)}
for (j in std.act.no) {gamma[j]~dnorm(0,prec.gamma)}
for (j in std.act.yes) {gamma[j]~dnorm(g.act,prec.gamma)}
g~dnorm(0, 0.01)
g.act~dnorm(0, 0.01)
tau.gamma~dunif(0,2)
 prec.gamma <- pow(tau.gamma,-2)
                      # bias adjustment
                      for (j in 1:(ns.ipd+ns.ad)) {R[j]~dbern(pi[bias_index[j]])}
                      pi[1]~dbeta(20,1)
                      pi[2]~dbeta(1,20)
                      pi[3]~dbeta(30,1)
                      pi[4]~dbeta(1,30)
                      pi[5]~dbeta(1,1)
  
  }
"

antidep_jags_adj2_dic <-"
model {
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
      
        theta1[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]],precd.ad[j,t.ad[j,k]])
      theta2[j+ns.ipd,t.ad[j,k]]~dnorm(md.ad[j,t.ad[j,k]]+gamma[j+ns.ipd],precd.ad[j,t.ad[j,k]]+(prec.gamma*2*(k-1)/k))
      theta[j+ns.ipd,t.ad[j,k]] <- (1-pi[bias_index[j]])*theta1[j+ns.ipd,t.ad[j,k]]+pi[bias_index[j]]*theta2[j+ns.ipd,t.ad[j,k]]
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
  
  
  
  
  # Random effect for gamma (bias effect)
for (j in std.in) {gamma[j]~dnorm(g,prec.gamma)}
for (j in std.act.no) {gamma[j]~dnorm(0,prec.gamma)}
for (j in std.act.yes) {gamma[j]~dnorm(g.act,prec.gamma)}
g~dnorm(0, 0.01)
g.act~dnorm(0, 0.01)
tau.gamma~dunif(0,2)
 prec.gamma <- pow(tau.gamma,-2)
                      # bias adjustment
                      pi[1]~dbeta(20,1)
                      pi[2]~dbeta(1,20)
                      pi[3]~dbeta(30,1)
                      pi[4]~dbeta(1,30)
                      pi[5]~ dbeta(1,1)
  
  }
"
