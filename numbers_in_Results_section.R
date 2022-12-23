# -------- install required libraries --------#
#install.packages("crossnma")
#install.packages("metafor")
#install.packages("dplyr")
#install.packages("ggplot2")

# -------- load libraries --------#
library(crossnma)
library(metafor) # for fig2: forest plot
library(dplyr) # for fig2
library(ggplot2) # fig3 and fig4

#-------- load JAGS output --------#

#** RRMS
load("output/RRMS/JAGS/jagsfit_rrms_NMA.RData") # has jagsfit_rrms_NMA
load("output/RRMS/JAGS/jagsfit_rrms_adjust1_NMR_age.RData") # has jagsfit_rrms_adjust1_NMR_age
load("output/RRMS/JAGS/jagsfit_RRMS_prior.RData") # has jagsfit.prior

#** antidepressant 
load("output/antidepressant/JAGS/jagsfit_antidep_main.RData") # jags output
load("output/antidepressant/JAGS/jagsfit_antidep_sens.RData") # 

#** RRMS 
# mean bias effect g from adjust1 and adjust2
summary(jagsfit_rrms_NMA[[2]])
summary(jagsfit_rrms_NMA[[3]])

# OR of NZ vs PL for studies at high RoB 
0.32*0.705

# age effect in adjust 1
summary(jagsfit_rrms_adjust1_NMR_age)

#** Antidepressant 
# Main analysis: mean bias effect g_p and  g_act from adjust1 and adjust2
summary(jagsfit_antidep_main[[2]]) # adjust1
summary(jagsfit_antidep_main[[3]]) # adjust2
plot(jagsfit_antidep_main[[2]])
plot(jagsfit_antidep_main[[2]]$samples[[2]][,'g'], type = "l",ylim=c(-0.4,0.4))
lines(jagsfit_antidep_main[[2]]$samples[[1]][,'g'],col=2)

FE_NMA_binomial_probitplot(jagsfit_antidep_main[[2]]$samples[[2]][,'g.act'], type = "l",ylim=c(-0.4,0.4))
lines(jagsfit_antidep_main[[2]]$samples[[1]][,'g.act'],col=2)
# Sensitivity analysis: mean bias effect g_p and  g_act from adjust1 and adjust2
summary(jagsfit_antidep_sens[[1]]) # adjust1, q=0.25
summary(jagsfit_antidep_sens[[2]]) # adjust2, q=0.25
