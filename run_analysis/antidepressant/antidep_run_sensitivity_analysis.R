#  Cross-Network meta-analysis of antidepressants ** Sensitivity analysis **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file use antidepressnats to run these analysis (check robust of the estimates for changeing the prior of bias probability)
# 1. bias-adjust1 with Beta(1,10) and Beta(10,1) to model bias probability "Sensitivity analysis"
# 2. bias-adjust2 with Beta(1,10) and Beta(10,1) to model bias probability "Sensitivity analysis"


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

#-------- Main analysis --------#

#** 1. bias-adjust1 with Beta(1,10) and Beta(10,1)
# jags model: code+data
mod_adjust1_sens <- crossnma.model(prt.data=NULL,
                                   std.data=gris,
                                   trt=drug,
                                   study=study,
                                   outcome=r,
                                   n=n,
                                   design=design,
                                   reference='Plac',
                                   trt.effect='random',
                                   #---------- bias adjustment ----------
                                   method.bias='adjust1',
                                   bias=rob,
                                   unfav=unfav,
                                   bias.group = bias.group,
                                   bias.type='add',
                                   bias.effect='random',
                                   prior =list(pi.high.rct="dbeta(10,1)",
                                               pi.low.rct="dbeta(1,10)"))


# run jags
fit_adjust1_sens <- crossnma(mod_adjust1_sens,
                                 n.adapt = n.adapt,
                                 n.iter=n.iter,
                                 n.burnin = n.burnin,
                                 thin=thin,
                                 n.chains=n.chains,
                                 monitor="g.act")

#** 2. bias-adjust2 with Beta(1,10) and Beta(10,1)
# jags model: code+data
mod_adjust2_sens <- crossnma.model(prt.data=NULL,
                                   std.data=gris,
                                   trt=drug,
                                   study=study,
                                   outcome=r,
                                   n=n,
                                   design=design,
                                   reference='Plac',
                                   trt.effect='random',
                                   #---------- bias adjustment ----------
                                   method.bias='adjust2',
                                   bias=rob,
                                   unfav=unfav,
                                   bias.group = bias.group,
                                   # bias.type='add',
                                   bias.effect='random',
                                   prior =list(pi.high.rct="dbeta(10,1)",
                                               pi.low.rct="dbeta(1,10)")
)


# run jags
fit_adjust2_sens <- crossnma(mod_adjust2_sens,
                                 n.adapt = n.adapt,
                                 n.iter=n.iter,
                                 n.burnin = n.burnin,
                                 thin=thin,
                                 n.chains=n.chains,
                                 monitor="g.act")

#== SAVE output
jagsfit_antidep_sens <- list(fit_adjust1_sens,fit_adjust2_sens)
save(jagsfit_antidep_sens,file = "output/antidepressant/JAGS/jagsfit_antidep_sens.RData")





