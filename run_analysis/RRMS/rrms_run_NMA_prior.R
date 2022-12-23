#  Cross-Network meta-analysis of relapsing-remitting multiple sclerosis RRMS ** Use NRS as a prior **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file apply the prior method to RRMS: 

# we conduct NMA with only the RCTs. The priors in the RCT model are 
 # the posteriors of the relative treatment effects estimated from the model with NRS data
# This model is fitted repeatedly with different downweight values (w) of NRS posteriors to the RCT model

#-------- load the libray --------#
#install.packages("crossnma")
library(crossnma)

#-------- data --------#
myprt.data <- read.csv("data/RRMS/final data/rrms_final_IPD") 
mystd.data <- read.csv("data/RRMS/final data/rrms_final_AD") 

#-------- MCMC settings --------#
n.adapt = 2000
n.iter=100000
n.burnin = 40000
thin=1
n.chains=2
#-------- Use NRS as a prior with different weights w --------#

#--------  prior NMR with different w (inflation factor of variance) --------#
w <- seq(0,1,l=11)[-1]
jagsfit.prior <- list()
for (i in 1:length(w)) {
  # jags model: code+data
  mod.prior <- crossnma.model(trt=treat,study=study,
                              outcome=r,n=n,design=design,
                              prt.data=myprt.data,
                              std.data=mystd.data,
                              reference='Natalizumab',
                              trt.effect='common',
                              #---------- bias adjustment ----------
                              method.bias='prior',
                              run.nrs=list(var.infl=w[i],
                                           mean.shift=0,
                                           trt.effect="common",
                                           n.adapt = n.adapt,
                                           n.iter=n.iter,
                                           n.burnin = n.burnin,
                                           thin=thin,
                                           n.chains=n.chains)
                              )
  
  jagsfit.prior[[i]] <- crossnma(mod.prior,
                                     n.adapt = n.adapt,
                                     n.iter=n.iter,
                                     n.burnin = n.burnin,
                                     thin=thin,
                                     n.chains=n.chains,
                                     monitor='LOR')
}


#== SAVE output
save(jagsfit.prior,file = "output/RRMS/JAGS/jagsfit_RRMS_prior.RData")





