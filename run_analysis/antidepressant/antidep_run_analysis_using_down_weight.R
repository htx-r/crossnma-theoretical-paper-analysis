#  Cross-Network meta-analysis of antidepressants ** Additional analysis with q **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file use antidepressnats to run these analysis

# 1. bias-adjust1 with q=0.25 "Additional analysis with q"
# 2. bias-adjust1 with q=0.8 "Additional analysis with q"
# 3. bias-adjust2 with q=0.25  "Additional analysis with q"
# 4. bias-adjust2 with q=0.8  "Additional analysis with q"

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

#-------- Analysis using the down_weight factor q --------#

#** 1. bias-adjust1 with q=0.25
# jags model: code+data
mod_adjust1_q1 <- crossnma.model(prt.data=NULL,
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
                                 down.wgt=0.25,
                                 unfav=unfav,
                                 bias.group = bias.group,
                                 bias.type='add',
                                 bias.effect='random')

# run jags
fit_adjust1_q1 <- crossnma(mod_adjust1_q1,
                                 n.adapt = n.adapt,
                                 n.iter=n.iter,
                                 n.burnin = n.burnin,
                                 thin=thin,
                                 n.chains=n.chains,
                                 monitor="g.act")

#** 2. bias-adjust1 with q=0.8
# jags model: code+data
mod_adjust1_q2 <- crossnma.model(prt.data=NULL,
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
                                 down.wgt=0.8,
                                 unfav=unfav,
                                 bias.group = bias.group,
                                 bias.type='add',
                                 bias.effect='random')


# run jags
fit_adjust1_q2 <- crossnma(mod_adjust1_q2,
                                 n.adapt = n.adapt,
                                 n.iter=n.iter,
                                 n.burnin = n.burnin,
                                 thin=thin,
                                 n.chains=n.chains,
                                 monitor="g.act")
#** 3. bias-adjust2 with q=0.25
mod_adjust2_q1 <- crossnma.model(prt.data=NULL,
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
                       down.wgt=0.25,
                       unfav=unfav,
                       bias.group = bias.group,
                       bias.effect='random')



# run jags
fit_adjust2_q1 <- crossnma(mod_adjust2_q1,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor="g.act")
#** 4. bias-adjust2 with q=0.8
mod_adjust2_q2 <- crossnma.model(prt.data=NULL,
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
                                 down.wgt=0.8,
                                 unfav=unfav,
                                 bias.group = bias.group,
                                 bias.effect='random')



# run jags
fit_adjust2_q2 <- crossnma(mod_adjust2_q2,
                               n.adapt = n.adapt,
                               n.iter=n.iter,
                               n.burnin = n.burnin,
                               thin=thin,
                               n.chains=n.chains,
                               monitor="g.act")
#== SAVE output
jagsfit_antidep_down_weight <- list(fit_adjust1_q1,fit_adjust1_q2,fit_adjust2_q1,fit_adjust2_q2)
save(jagsfit_antidep_down_weight,file = "output/antidepressant/JAGS/jagsfit_antidep_down_weight.RData")





