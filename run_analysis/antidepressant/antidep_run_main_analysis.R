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
source("functions/JAGSmodels/antidep_jags_add_dic.R") # JAGS models with deviance
#** 1. unadjusted 
# jags model: code+data
mod_unadjust <- crossnma.model(prt.data=NULL,
                       std.data=gris,
                       trt=drug,
                       study=study,
                       outcome=r,
                       n=n,
                       design=design,
                       reference="Plac",
                       trt.effect="random",
                       method.bias = "naive"
                       )


mod_unadjust$jagsmodel <- antidep_jags_unadj_dic #  DIC to JAGS model

# run jags
fit_unadjust <- crossnma.run(model=mod_unadjust,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor = c("dev.ipd","resdev.ipd","totresdev.ipd","rhat.ipd",
                                     "dev.ad","resdev.ad","totresdev.ad","rhat.ad")
                         )

#** 2. bias-adjust1 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust1_main <- crossnma.model(prt.data=NULL,
                       std.data=gris,
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
                                   pi.low.rct="dbeta(1,20)"))

mod_adjust1_main$jagsmodel <- antidep_jags_adj1_dic
# run jags
fit_adjust1_main <- crossnma.run(model=mod_adjust1_main,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor=c("g.act","dev.ipd","resdev.ipd","totresdev.ipd","rhat.ipd",
                                   "dev.ad","resdev.ad","totresdev.ad","rhat.ad"))

#** 3. bias-adjust2 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust2_main <- crossnma.model(prt.data=NULL,
                       std.data=gris,
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
                                   pi.low.rct="dbeta(1,20)")
                       )

mod_adjust2_main$jagsmodel <- antidep_jags_adj2_dic
# run jags
fit_adjust2_main <- crossnma.run(model=mod_adjust2_main,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor=c("g.act","dev.ipd","resdev.ipd","totresdev.ipd","rhat.ipd",
                                   "dev.ad","resdev.ad","totresdev.ad","rhat.ad"))

#== SAVE output
jagsfit_antidep_main <- list(fit_unadjust,fit_adjust1_main,fit_adjust2_main)
save(jagsfit_antidep_main,file = "output/antidepressant/JAGS/jagsfit_antidep_main.RData")

#== SAVE output - Revision - with deviance
jagsfit_antidep_main_R <- list(fit_unadjust,fit_adjust1_main,fit_adjust2_main)
save(jagsfit_antidep_main_R,file ="output/antidepressant/JAGS/jagsfit_antidep_main_after_revision.RData")




