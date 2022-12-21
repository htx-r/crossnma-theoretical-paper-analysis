# -------- install required libraries --------#
#install.packages("crossnma")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("remotes")
#remotes::install_github("audrey-b/BUGSnet")

# -------- load libraries --------#
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
load("output/antidepressant/JAGS/jagsfit_antidep_main_after_revision.RData") # jags output
load("output/antidepressant/JAGS/jagsfit_antidep_sens.RData") # 

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

# Appendix Figure 5 & Appendix Figure 6
### Traceplot for g and g.act for antidepressants network - adjust 1 and 2

#** bias1
# g (mean bias for active-inactive comparison)
plot(as.vector(jagsfit_antidep_main_R[[2]]$samples[[2]][,'g']), type = "l",
     ylab="", xlab="",las=1, main="Traceplot of g (adjust1)")
lines(as.vector(jagsfit_antidep_main_R[[2]]$samples[[1]][,'g']),col=2)

# g.act (mean bias for active-active comparison, sponsored)
plot(as.vector(jagsfit_antidep_main_R[[2]]$samples[[2]][,'g.act']), type = "l",
     ylab="", xlab="",las=1, main="Traceplot of g.act (adjust1)")
lines(as.vector(jagsfit_antidep_main_R[[2]]$samples[[1]][,'g.act']),col=2)

# g (mean bias for active-inactive comparison)
plot(as.vector(jagsfit_antidep_main_R[[3]]$samples[[2]][,'g']), type = "l",
     ylab="", xlab="",las=1, main="Traceplot of g (adjust2)")
lines(as.vector(jagsfit_antidep_main_R[[3]]$samples[[1]][,'g']),col=2)

# g.act (mean bias for active-active comparison, sponsored)
plot(as.vector(jagsfit_antidep_main_R[[3]]$samples[[2]][,'g.act']), type = "l",
     ylab="", xlab="",las=1, main="Traceplot of g.act (adjust2)")
lines(as.vector(jagsfit_antidep_main[[3]]$samples[[1]][,'g.act']),col=2)

# Appendix Figure 7
source("functions/plot and table/antidep_forestplot.R")
antidep <- read.csv("data/antidepressant/raw data/antidepressant")# to get drug names from the raw data
jagsfit_antidep_unadjust_sens <- list(jagsfit_antidep_main_R[[1]],jagsfit_antidep_sens[[1]],jagsfit_antidep_sens[[2]])
antidep_forestplot(jagsfit_antidep_unadjust_sens, antidep)

#-------- Appendix Tables --------#

#-------- Appendix Table 1, 2, 3
source("functions/plot and table/rrms_net_charac_tab.R")
network.char <- rrms_net_charac_tab(myprt.data,mystd.data)
network.char$network
network.char$intervention
network.char$comparison 

# Appendix Table 4, 5, 6, 7
heatplot(jagsfit_rrms_NMA[[1]],exp = T)
heatplot(jagsfit.prior[[10]],exp = T,
         order = c("Placebo", "Dimethyl fumarate",
                   "Glatiramer acetate", "Natalizumab"))
heatplot(jagsfit_rrms_NMA[[2]],exp = T)
heatplot(jagsfit_rrms_NMA[[3]],exp = T)


# Appendix Table 8
source("functions/plot and table/rrms_conv_stat.R")
conv.tab1 <- rrms_conv_stat_tab(jagsfit_rrms_NMA[[1]], method = "unadjusted")
conv.tab2 <-rrms_conv_stat_tab(jagsfit.prior[[10]], method = "prior")
conv.tab3 <-rrms_conv_stat_tab(jagsfit_rrms_NMA[[2]], method = "bias-adjusted 1")
conv.tab4 <-rrms_conv_stat_tab(jagsfit_rrms_NMA[[3]], method = "bias-adjusted 2")

rbind(conv.tab1, conv.tab2, conv.tab3, conv.tab4)
