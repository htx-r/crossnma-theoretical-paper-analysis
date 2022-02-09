#  Cross-Network meta-analysis of relapsing-remitting multiple sclerosis RRMS ** Main analysis **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file use RRMS to run these analysis

# 1. unadjusted 
# 2. bias-adjust1 
# 3. bias-adjust2 

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

#** 1. unadjusted 
# jags model: code+data
mod_unadjust <- crossnma.model(prt.data=myprt.data,
                                 std.data=mystd.data,
                                 trt='treat',
                                 study='study',
                                 outcome='r',
                                 n='n',
                                 design='design',
                                 reference='Placebo',
                                 trt.effect='common',
                                 method.bias = 'naive'
                                 )

# run jags
fit_unadjust <- crossnma.run(model=mod_unadjust,
                                 n.adapt = n.adapt,
                                 n.iter=n.iter,
                                 n.burnin = n.burnin,
                                 thin=thin,
                                 n.chains=n.chains,
                                 monitor="LOR")

#** 2. bias-adjust1 
myprt.data$design <- "rct" # I want to specify for the RCT at high RoB and the unique NRS the same prior of probability of bias (pi)
# jags model: code+data
mod_adjust1 <- crossnma.model(prt.data=myprt.data,
                                 std.data=mystd.data,
                                 trt='treat',
                                 study='study',
                                 outcome='r',
                                 n='n',
                                 design='design',
                                 reference='Placebo',
                                 trt.effect='common',
                                 #---------- bias adjustment ----------
                                 method.bias='adjust1',
                                 bias='bias',
                                 bias.type='add',
                                 bias.effect='common',
                                 unfav="unfav",
                                 bias.group="bias.group",
                                 prior =list(pi.high.rct="dbeta(100,1)",
                                          pi.low.rct="dbeta(1,100)")
                                 )

# run jags
fit_adjust1 <- crossnma.run(model=mod_adjust1,
                                 n.adapt = n.adapt,
                                 n.iter=n.iter,
                                 n.burnin = n.burnin,
                                 thin=thin,
                                 n.chains=n.chains,
                                 monitor="LOR")
#** 3. bias-adjust2 
mod_adjust2 <- crossnma.model(prt.data=myprt.data,
                                 std.data=mystd.data,
                                 trt='treat',
                                 study='study',
                                 outcome='r',
                                 n='n',
                                 design='design',
                                 reference='Placebo',
                                 trt.effect='common',
                                 #---------- bias adjustment ----------
                                 method.bias='adjust2',
                                 bias='bias',
                                 bias.effect='common',
                                 unfav="unfav",
                                 bias.group="bias.group",
                                 prior =list(pi.high.rct="dbeta(100,1)",
                                          pi.low.rct="dbeta(1,100)")
                                 )

# run jags
fit_adjust2 <- crossnma.run(model=mod_adjust2,
                         n.adapt = n.adapt,
                         n.iter=n.iter,
                         n.burnin = n.burnin,
                         thin=thin,
                         n.chains=n.chains,
                         monitor="LOR")

#== SAVE output
jagsfit_rrms_NMA <- list(fit_unadjust,fit_adjust1,fit_adjust2)
save(jagsfit_rrms_NMA,file = "output/RRMS/JAGS/jagsfit_rrms_NMA.RData")





