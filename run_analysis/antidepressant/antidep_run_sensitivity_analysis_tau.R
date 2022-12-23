#  Cross-Network meta-analysis of antidepressants ** Main analysis **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file use antidepressnats to run these analysis
# 1. Unadjusted                                  "Main analysis"
# 2. bias-adjust1 with Beta(1,20) and Beta(20,1) to model bias probability "Main analysis"
# 3. bias-adjust2 with Beta(1,20) and Beta(20,1) to model bias probability "Main analysis"

#-------- load the libray --------#
#install.packages("crossnma")
library(crossnma)

#-------- data --------#
gris <- read.csv("data/antidepressant/final data/antidepressant_for_main_analysis")

#-------- MCMC settings --------#
n.adapt = 2000
n.iter=100000
n.burnin = 40000
thin=1
n.chains=2

#** 1. unadjusted 
# jags model: code+data
mod_unadjust_tau1 <-
  crossnma.model(std.data=gris,
                 trt=drug,
                 study=study,
                 outcome=r,
                 n=n,
                 design=design,
                 reference="Plac",
                 trt.effect="random",
                 method.bias = "naive",
                 prior=list(tau.trt="dunif(0,5)")
                 )


# run jags
fit_unadjust_tau1 <- crossnma(mod_unadjust_tau1,
                              n.adapt = n.adapt,
                              n.iter=n.iter,
                              n.burnin = n.burnin,
                              thin=thin,
                              n.chains=n.chains
                              )
summary(fit_unadjust_tau1) # 0.209, (0.168, 0.250)
#** 2. bias-adjust1 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust1_tau1 <- crossnma.model(std.data=gris,
                                   trt=drug,
                                   study=study,
                                   outcome=r,
                                   n=n,
                                   design=design,
                                   reference="Plac",
                                   trt.effect="random",
                                   #---------- bias adjustment ----------
                                   method.bias='adjust1',
                                   bias=rob,
                                   unfav=unfav,
                                   bias.group = bias.group,
                                   bias.type='add',
                                   bias.effect='random',
                                   prior =list(pi.high.rct="dbeta(20,1)",
                                               pi.low.rct="dbeta(1,20)",
                                               tau.trt="dunif(0,5)"))


# run jags
fit_adjust1_tau1 <- crossnma(mod_adjust1_tau1,
                             n.adapt = n.adapt,
                             n.iter=n.iter,
                             n.burnin = n.burnin,
                             thin=thin,
                             n.chains=n.chains
                             )

summary(fit_adjust1_tau1) # tau.trt: 0.167, (0.087, 0.228) 
#** 3. bias-adjust2 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust2_tau1 <- crossnma.model(std.data=gris,
                                   trt=drug,
                                   study=study,
                                   outcome=r,
                                   n=n,
                                   design=design,
                                   reference="Plac",
                                   trt.effect="random",
                                   #---------- bias adjustment ----------
                                   method.bias='adjust2',
                                   unfav=unfav,
                                   bias = rob,
                                   bias.group = bias.group,
                                   bias.type='add',
                                   bias.effect='random',
                                   prior =list(pi.high.rct="dbeta(20,1)",
                                               pi.low.rct="dbeta(1,20)",
                                               tau.trt="dunif(0,5)"
                                               ))

# run jags
fit_adjust2_tau1 <- crossnma(mod_adjust2_tau1,
                             n.adapt = n.adapt,
                             n.iter=n.iter,
                             n.burnin = n.burnin,
                             thin=thin,
                             n.chains=n.chains,
                             monitor=c("g.act")
                             )
summary(fit_adjust2_tau1) # 0.212, (0.144,0.289) 
#== SAVE output
jagsfit_antidep_sens_tau1 <- list(fit_unadjust_tau1,fit_adjust1_tau1,fit_adjust2_tau1)
save(jagsfit_antidep_sens_tau1,file = "output/antidepressant/JAGS/jagsfit_antidep_sens_tau1.RData")
