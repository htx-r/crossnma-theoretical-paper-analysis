#  Cross-Network meta-analysis of antidepressants ** Main analysis **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file use antidepressnats to run these analysis
# 1. Unadjusted                                  "Main analysis"
# 2. bias-adjust1 with Beta(1,20) and Beta(20,1) to model bias probability "Main analysis"
# 3. bias-adjust2 with Beta(1,20) and Beta(20,1) to model bias probability "Main analysis"

#-------- load the libray --------#
devtools::install_github("htx-r/crossnma",force = TRUE)
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

#** 1. unadjusted 
# jags model: code+data
mod_unadjust <- crossnma.model(prt.data=NULL,
                       std.data=gris,
                       trt=c('drug'),
                       study=c('study'),
                       outcome=c('r'),
                       n='n',
                       design=c('design'),
                       reference='Plac',
                       trt.effect='random',
                       method.bias = 'naive'
                       )

# run jags
fit_unadjust <- crossnma.run(model=mod_unadjust,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains
                         )

#** 2. bias-adjust1 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust1_main <- crossnma.model(prt.data=NULL,
                       std.data=gris,
                       trt=c('drug'),
                       study=c('study'),
                       outcome=c('r'),
                       n='n',
                       design=c('design'),
                       reference='Plac',
                       trt.effect='random',
                       #---------- bias adjustment ----------
                       method.bias='adjust1',
                       bias=c('rob'),
                       unfav=c("unfav"),
                       bias.group = "bias.group",
                       bias.type='add',
                       bias.effect='random',
                       prior =list(pi.high.rct="dbeta(20,1)",
                                   pi.low.rct="dbeta(1,20)"))


# run jags
fit_adjust1_main <- crossnma.run(model=mod_adjust1_main,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor="g.act")

#** 3. bias-adjust2 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust2_main <- crossnma.model(prt.data=NULL,
                       std.data=gris,
                       trt=c('drug'),
                       study=c('study'),
                       outcome=c('r'),
                       n='n',
                       design=c('design'),
                       reference='Plac',
                       trt.effect='random',
                       #---------- bias adjustment ----------
                       method.bias='adjust2',
                       bias=c('rob'),
                       unfav=c("unfav"),
                       bias.group = "bias.group",
                       bias.effect='random',
                       prior =list(pi.high.rct="dbeta(20,1)",
                                   pi.low.rct="dbeta(1,20)")
                       )


# run jags
fit_adjust2_main <- crossnma.run(model=mod_adjust2_main,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor="g.act")

#== SAVE output
jagsfit_antidep_main <- list(fit_unadjust,fit_adjust1_main,fit_adjust2_main)
save(jagsfit_antidep_main,file = "output/antidepressant/JAGS/jagsfit_antidep_main.RData")





