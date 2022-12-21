#  Cross-Network meta-regression of relapsing-remitting multiple sclerosis RRMS
#  Authors: Tasnim Hamza and Georgia Salanti

# This file fits NMR model to RRMS using age as a covaraite. we use bias-adjust method to adjust for biases

# -------- load the libray --------#
devtools::install_github("htx-r/crossnma",force = TRUE)
library(crossnma)

#-------- data --------#
myprt.data <- read.csv("data/RRMS/final data/rrms_final_IPD") 
mystd.data <- read.csv("data/RRMS/final data/rrms_final_AD") 

# center age around the mean age of IPDs (38)
mystd.data$age <- mystd.data$age-round(mean(myprt.data$age)) 
myprt.data$age <- myprt.data$age-round(mean(myprt.data$age))

#-------- MCMC settings --------#
n.adapt = 2000
n.iter=100000
n.burnin = 40000
thin=1
n.chains=2

#-------- Network meta-regression using age --------#
myprt.data$design <- "rct" # I want to assign the same prior of probability of bias (pi) for the RCT at high RoB and the unique NRS 
mod_adjust1_age <- crossnma.model(prt.data=myprt.data,
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
                                  unfav="unfav",
                                  bias.type='add',
                                  bias.effect='common',
                                  prior =list(pi.high.rct="dbeta(100,1)",
                                              pi.low.rct="dbeta(1,100)"),
                                  #----------  meta-regression ----------
                                  covariate = 'age',
                                  cov.ref=0,
                                  split.regcoef = F,
                                  regb.effect='common',
                                  regw.effect='common'
                                  )

# run jags
jagsfit_rrms_adjust1_NMR_age <- crossnma.run(model=mod_adjust1_age,
                          n.adapt = n.adapt,
                          n.iter=n.iter,
                          n.burnin = n.burnin,
                          thin=thin,
                          n.chains=n.chains,
                          monitor='LOR')
#== SAVE output
save(jagsfit_rrms_adjust1_NMR_age,file = "output/RRMS/JAGS/jagsfit_rrms_adjust1_NMR_age.RData")



