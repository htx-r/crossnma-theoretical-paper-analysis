#  Revision: Cross-Network meta-analysis of relapsing-remitting multiple sclerosis RRMS 
# ** Main analysis with deviance and OR for each study **
#  Authors: Tasnim Hamza and Georgia Salanti

# 1. unadjusted 
# 2. use NRS as a prior
# 3. bias-adjust1 
# 4. bias-adjust2 

#-------- load the libray --------#
devtools::install_github("htx-r/crossnma",force = TRUE)
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

#-------- NMA  --------#
# load all JAGS models (include the DIC part)
source("functions/JAGSmodels/rrms_jags_add_dic.R")
#** 1. unadjusted 
# jags model: code+data
mod_unadjust <- crossnma.model(trt=treat,study=study,
                               outcome=r,n=n,design=design,
                               prt.data=myprt.data,
                               std.data=mystd.data,
                               reference='Placebo',
                               trt.effect='common',
                               method.bias = 'naive'
                               )

mod_unadjust$jagsmodel <- rrms_jags_unadj_dic #   DIC part added to JAGS model
# r.ipd and n.ipd needed to calculate DIC and deviance 
source("run_analysis/RRMS/R1_rrms_run_NMA_prep.R")
extra_data <- find_r_n_mat_ipd(mod_unadjust)
mod_unadjust$data$r.ipd <-extra_data$r.ipd
mod_unadjust$data$n.ipd <-extra_data$n.ipd
mod_unadjust$data$arm <- extra_data$arm
mod_unadjust$data$stnd <-extra_data$stnd

# run jags
fit_unadjust <- crossnma.run(model=mod_unadjust,
                             n.adapt = n.adapt,
                             n.iter=n.iter,
                             n.burnin = n.burnin,
                             thin=thin,
                             n.chains=n.chains,
                             monitor=c("LOR","totresdev.ipd","totresdev.ad","theta"))

totresdev.ipd.unadj <- mean(c(fit_unadjust$samples[[1]][,"totresdev.ipd"],fit_unadjust$samples[[2]][,"totresdev.ipd"])) # 3351866 # 1473048
totresdev.ipd.unadj
#** 2. use NRS as a prior
# jags model: code+data
mod_prior <- crossnma.model(trt=treat,study=study,
                               outcome=r,n=n,design=design,
                               prt.data=myprt.data,
                               std.data=mystd.data,
                               reference='Natalizumab',
                               trt.effect='common',
                               method.bias = 'prior')
mod_prior$jagsmodel <- rrms_jags_prior_dic #   DIC part added to JAGS model
extra_data <- find_r_n_mat_ipd(mod_prior)
mod_prior$data$r.ipd <-extra_data$r.ipd
mod_prior$data$n.ipd <-extra_data$n.ipd
mod_prior$data$arm <- extra_data$arm
mod_prior$data$stnd <-extra_data$stnd
# run jags
fit_prior <- crossnma.run(model=mod_prior,
                             n.adapt = n.adapt,
                             n.iter=n.iter,
                             n.burnin = n.burnin,
                             thin=thin,
                             n.chains=n.chains,
                             monitor=c("LOR","totresdev.ipd","totresdev.ad","theta"))
totresdev.ipd.prior <- mean(c(fit_prior$samples[[1]][,"totresdev.ipd"],fit_prior$samples[[2]][,"totresdev.ipd"])) # 3351866 # 1473048
totresdev.ipd.prior
#** 3. bias-adjust1 
myprt.data$design <- "rct" # I want to specify for the RCT at high RoB and the unique NRS the same prior of probability of bias (pi)
# jags model: code+data
mod_adjust1 <- crossnma.model(trt=treat,study=study,
                              outcome=r,n=n,design=design,
                              prt.data=myprt.data,
                              std.data=mystd.data,
                              reference='Placebo',
                              trt.effect='common',
                              #---------- bias adjustment ----------
                              method.bias='adjust1',
                              bias=bias,
                              bias.type='add',
                              bias.effect='common',
                              unfav=unfav,
                              bias.group=bias.group,
                              prior =list(pi.high.rct="dbeta(100,1)",
                                          pi.low.rct="dbeta(1,100)")
                              )
mod_adjust1$jagsmodel <- rrms_jags_adj1_dic #   DIC part added to JAGS model
extra_data <- find_r_n_mat_ipd(mod_adjust1)
mod_adjust1$data$r.ipd <-extra_data$r.ipd
mod_adjust1$data$n.ipd <-extra_data$n.ipd
mod_adjust1$data$arm <- extra_data$arm
mod_adjust1$data$stnd <-extra_data$stnd
# run jags
fit_adjust1 <- crossnma.run(model=mod_adjust1,
                            n.adapt = n.adapt,
                            n.iter=n.iter,
                            n.burnin = n.burnin,
                            thin=thin,
                            n.chains=n.chains,
                            monitor=c("LOR","totresdev.ipd","totresdev.ad","theta"))
totresdev.ipd.adjust1 <- mean(c(fit_adjust1$samples[[1]][,"totresdev.ipd"],fit_adjust1$samples[[2]][,"totresdev.ipd"])) # 3351866 # 1473048

#** 4. bias-adjust2 
mod_adjust2 <- crossnma.model(trt=treat,study=study,
                              outcome=r,n=n,design=design,
                              prt.data=myprt.data,
                              std.data=mystd.data,
                              reference='Placebo',
                              trt.effect='common',
                              #---------- bias adjustment ----------
                              method.bias='adjust2',
                              bias=bias,
                              bias.effect='common',
                              unfav=unfav,
                              bias.group=bias.group,
                              prior =list(pi.high.rct="dbeta(100,1)",
                                          pi.low.rct="dbeta(1,100)")
                              )
mod_adjust2$jagsmodel <- rrms_jags_adj2_dic #   DIC part added to JAGS model
extra_data <- find_r_n_mat_ipd(mod_adjust2)
mod_adjust2$data$r.ipd <-extra_data$r.ipd
mod_adjust2$data$n.ipd <-extra_data$n.ipd
mod_adjust2$data$arm <- extra_data$arm
mod_adjust2$data$stnd <-extra_data$stnd
# run jags
fit_adjust2 <- crossnma.run(model=mod_adjust2,
                            n.adapt = n.adapt,
                            n.iter=n.iter,
                            n.burnin = n.burnin,
                            thin=thin,
                            n.chains=n.chains,
                            monitor=c("LOR","totresdev.ipd","totresdev.ad","theta"))


#== SAVE output
jagsfit_rrms_NMA_dic <- list(fit_unadjust,fit_prior,fit_adjust1,fit_adjust2)
save(jagsfit_rrms_NMA_dic,file = "output/RRMS/JAGS/jagsfit_rrms_NMA_dic.RData")


