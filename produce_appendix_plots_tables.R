#-------- load the libray --------#
devtools::install_github("htx-r/crossnma",force = TRUE)
library(crossnma)
library(ggplot2) # appx fig1, 2,3
library(dplyr) # appx fig1, 2,3
library(BUGSnet) # for appx tab1,2,3
#-------- data --------# (needed to create Appendix table 1,2,3)
myprt.data <- read.csv("data/RRMS/final data/rrms_final_IPD") 
mystd.data <- read.csv("data/RRMS/final data/rrms_final_AD") 

#-------- load JAGS output --------#
# RRMS
load("output/RRMS/JAGS/jagsfit_rrms_NMA.RData") # for appx fig4 a,b& appx table 4-6
load("output/RRMS/JAGS/jagsfit_RRMS_prior.RData") # for appx fig1

#** antidepressant 
load("output/antidepressant/JAGS/jagsfit_antidep_main.RData") # to get the unadjusted results
load("output/antidepressant/JAGS/jagsfit_antidep_sens.RData") # for appx fig5

#-------- Appendix Figures --------#
# Appendix Figure 1
source("functions/plot and table/rrms_ORvsW_prior.R")
rrms_ORvsW_prior(jagsfit.prior)
# Appendix Figure 2
source("functions/plot and table/rrms_beta_dist.R")
rrms_beta_dist(ab=100)
# Appendix Figure 3
rrms_beta_dist(ab=20)
# Appendix Figure 4 a and b
par(las=1)
coda::traceplot(jagsfit_rrms_NMA[[2]]$samples[,"g"]) # adjust 1
coda::traceplot(jagsfit_rrms_NMA[[3]]$samples[,"g"]) # adjust 2

# Appendix Figure 5
source("functions/plot and table/antidep_forestplot.R")
antidep <- read.csv("data/antidepressant/raw data/antidepressant")# to get drug names from the raw data
jagsfit_antidep_unadjust_sens <- list(jagsfit_antidep_main[[1]],jagsfit_antidep_sens[[1]],jagsfit_antidep_sens[[2]])
antidep_forestplot(jagsfit_antidep_unadjust_sens, antidep)

#-------- Appendix Tables --------#

#-------- Appendix Table 1, 2, 3
source("functions/plot and table/rrms_net_charac_tab.R")
network.char <- rrms_net_charac_tab(myprt.data,mystd.data)
network.char$network
network.char$intervention
network.char$comparison 

# Appendix Table 4, 5, 6, 7
crossnma.league(jagsfit_rrms_NMA[[1]],log.scale = F)$heatplot
crossnma.league(jagsfit.prior[[10]],log.scale = F)$heatplot
crossnma.league(jagsfit_rrms_NMA[[2]],log.scale = F)$heatplot
crossnma.league(jagsfit_rrms_NMA[[3]],log.scale = F)$heatplot


# 

