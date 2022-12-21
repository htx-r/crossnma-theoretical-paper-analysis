# -------- load the libray --------#
#Install required packages ----
#install.packages("crossnma")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("metafor")
#install.packages("coda")

library(crossnma)
library(metafor) # for fig2: forest plot
library(dplyr) # for fig2
library(ggplot2) # fig3 and fig4
library(coda)

#-------- data --------# (needed to create Fig1a and Appendix table 1,2,3)
myprt.data <- read.csv("data/RRMS/final data/rrms_final_IPD") 
mystd.data <- read.csv("data/RRMS/final data/rrms_final_AD") 

#-------- load JAGS output --------#

#** RRMS 
load("output/RRMS/JAGS/jagsfit_rrms_NMA.RData") # has jagsfit_rrms_NMA
load("output/RRMS/JAGS/jagsfit_rrms_adjust1_NMR_age.RData") # has jagsfit_rrms_adjust1_NMR_age
load("output/RRMS/JAGS/jagsfit_RRMS_prior.RData") # has jagsfit.prior

#** antidepressant 
load("output/antidepressant/JAGS/jagsfit_antidep_main_after_revision.RData") # jags output
load("output/antidepressant/JAGS/jagsfit_antidep_sens.RData") # Table 4
load("output/antidepressant/JAGS/jagsfit_antidep_down_weight.RData") # Table 4
# Figure 1: network plot

# a. RRMS
mod_rct <- crossnma.model(prt.data=myprt.data[myprt.data$design!="nrs",],
                               std.data=mystd.data,
                               trt=treat,
                               study=study,
                               outcome=r,
                               n=n,
                               design=design,
                               reference='Placebo',
                               trt.effect='common',
                               method.bias = 'naive') ## needed to create netplot for only RCTs to be able to add manually dashed lines to represent the NRS evidence 
#par(family = "Arial")
netgraph(mod_rct,
        start = "prcomp",scale=0.8,
        cex=1.95,
        labels=c("dimethyl fumarate",
        "glatiramer acetate",
        "natalizumab",
        "placebo"),
        plastic = FALSE, number = FALSE,
        lwd=10,
        points = TRUE, cex.points = 5, col.points = "darkred",
        multiarm = FALSE,
        col='grey'
        )

# b. antidepressant
#-------- data --------#
gris <- read.csv("data/antidepressant/final data/antidepressant_for_main_analysis")

mod_antidep <- crossnma.model(prt.data=NULL,
                               std.data=gris,
                               trt=drug,
                               study=study,
                               outcome=r,
                               n=n,
                               design=design,
                               reference="Plac",
                               trt.effect="random",
                               method.bias = "naive")
netgraph(mod_antidep, 
        cex=1.25,
        plastic = FALSE, thickness = "se.f", number = FALSE,
        points = TRUE, cex.points = 7, col.points = "darkred",
        multiarm = FALSE,
        col='grey',
        bg.number.of.studies =NULL,
        col.number.of.studies = 'black'
        )

# Figure 2: forest plot of OR (PL vs NZ, GA, DF) using 3 different methods of estimation
load("output/RRMS/JAGS/jagsfit_rrms_NMA_dic.RData") # 
source("functions/plot and table/rrms_forestplot_with_study_effect.R")
rrms_forestplot_with_study_effect(jagsfit_rrms_NMA_dic)
#source("functions/plot and table/rrms_forestplot.R")
#rrms_forestplot(jagsfit.nma=jagsfit_rrms_NMA)

# Figure 3
source("functions/plot and table/rrms_ORvsAge.R")
rrms_ORvsAge(jagsfit_rrms_adjust1_NMR_age, add.cri=T)

# Figure 4
source("functions/plot and table/antidep_forestplot.R")
#antidep <- read.csv("data/antidepressant/raw data/antidepressant")# to get drug names from the raw data
antidep_forestplot(jagsfit_antidep_main_R, gris)

# Table 3

#** For each study, find
# assigned treatments 
table(myprt.data$study, myprt.data$treat)
table(mystd.data$study, mystd.data$treat)
# number of relapsed patients 
myprt.data%>% group_by(study) %>% summarise_at("r", sum)
mystd.data%>% group_by(study) %>% summarise_at("r", sum)

# sample size
myprt.data%>% count(study) 
mystd.data%>% group_by(study) %>% summarise_at("n", sum)

# bias 
myprt.data%>% group_by(study, bias) %>% group_keys()
mystd.data%>% group_by(study, bias) %>% group_keys()

# mean age
myprt.data%>% group_by(study) %>% summarise_at("age", mean)
mystd.data%>% group_by(study) %>% summarise_at("age", mean)


# Table 4
bias1_main <- summary(jagsfit_antidep_main_R[[2]])
bias2_main <- summary(jagsfit_antidep_main_R[[3]])

bias1_sens <- summary(jagsfit_antidep_sens[[1]])
bias2_sens <- summary(jagsfit_antidep_sens[[2]])

bias1_q1 <- summary(jagsfit_antidep_down_weight[[1]])
bias1_q2 <- summary(jagsfit_antidep_down_weight[[2]])
bias2_q1 <- summary(jagsfit_antidep_down_weight[[3]])
bias2_q2 <- summary(jagsfit_antidep_down_weight[[4]])

list(bias1_main=bias1_main[c("g","g.act","tau.gamma"),c("Mean","2.5%","97.5%")],
     bias2_main=bias2_main[c("g","g.act","tau.gamma"),c("Mean","2.5%","97.5%")],
     
     bias1_sens=bias1_sens[c("g","g.act","tau.gamma"),c("Mean","2.5%","97.5%")],
     bias2_sens=bias2_sens[c("g","g.act","tau.gamma"),c("Mean","2.5%","97.5%")],
     
     bias1_q1=bias1_q1[c("g","g.act"),c("Mean","2.5%","97.5%")],
     bias2_q1=bias2_q1[c("g","g.act"),c("Mean","2.5%","97.5%")],
     bias1_q2=bias1_q2[c("g","g.act"),c("Mean","2.5%","97.5%")],
     bias2_q2=bias2_q2[c("g","g.act"),c("Mean","2.5%","97.5%")])

# Table 6 DIC values 
# Finding the DIC for the 3 models presenetd in the Main analysis
source("functions/plot and table/dic_fun.R")

# Calculate DIC 
dic.fun(jagsfit_antidep_main_R[[1]]) # unadjusted analysis   #  DIC= 2666.473, D_res=966.4946 (exactly what we got in BUGSnet)
dic.fun(jagsfit_antidep_main_R[[2]]) # bias-adjsuted 1 model #  DIC= , D_res=
dic.fun(jagsfit_antidep_main_R[[3]]) # bias-adjsuted 2 model #  DIC=  -6275.857

# Calculate total residual deviance
res.dev.unadjust.antidep <- mean(c(jagsfit_antidep_main_R[[1]]$samples[[1]][,"totresdev.ad"],jagsfit_antidep_main_R[[1]]$samples[[2]][,"totresdev.ad"])) # 966.4946
res.dev.adjust1.antidep <- mean(c(jagsfit_antidep_main_R[[2]]$samples[[1]][,"totresdev.ad"],jagsfit_antidep_main_R[[2]]$samples[[2]][,"totresdev.ad"]))  # 956.4426 
res.dev.adjust2.antidep <- mean(c(jagsfit_antidep_main_R[[3]]$samples[[1]][,"totresdev.ad"],jagsfit_antidep_main_R[[3]]$samples[[2]][,"totresdev.ad"]))  # 966.4664

#** RRMS: Calculate DIC for 
# load JAGS output
load("output/RRMS/JAGS/jagsfit_rrms_NMA_dic.RData") # 
# effective number of parameters (as we run common-effect model, pD=number of data points)
rrms_ipd <- read.csv("data/RRMS/final data/rrms_final_IPD") 
rrms_ad <- read.csv("data/RRMS/final data/rrms_final_AD") 

# 1. unadjust
sum_unadj <- summary(jagsfit_rrms_NMA_dic[[1]])

# Total residual deviance: IPD and AD
totresdev.ipd.unadj <- mean(c(jagsfit_rrms_NMA_dic[[1]]$samples[[1]][,"totresdev.ipd"],jagsfit_rrms_NMA_dic[[1]]$samples[[2]][,"totresdev.ipd"])) # 90482.96
# totresdev.ipd.unadj = INF! (sum_unadj["totresdev.ipd","Mean"])
# I compute the totresdev.ipd.unadj by calculating the mean of samples
totresdev.ad.unadj <- sum_unadj["totresdev.ad","Mean"] # 183.7683

# Effective number of parameters: IPD and AD: number of basic parameters(num of treat -1)+number of baselines (studies)+add 1 for each 3-arm trial
pD.ipd.unadj <- length(unique(rrms_ipd$study))+length(unique(rrms_ipd$treat))-1+2
pD.ad.unadj <- length(unique(rrms_ad$id))+length(unique(rrms_ad$treat))-1

# DIC
dic.ipd.unadj <- totresdev.ipd.unadj+pD.ipd.unadj
dic.ad.unadj <- totresdev.ad.unadj+pD.ad.unadj

# 2. prior 
sum_prior <- summary(jagsfit_rrms_NMA_dic[[2]])
# Total residual deviance: IPD and AD
totresdev.ad.prior <- sum_prior["totresdev.ad","Mean"] # 139.0034
totresdev.ipd.prior <- mean(c(jagsfit_rrms_NMA_dic[[2]]$samples[[1]][,"totresdev.ipd"],
                              jagsfit_rrms_NMA_dic[[2]]$samples[[2]][,"totresdev.ipd"])) 
# 87135.31
# Effective number of parameters: IPD and AD: number of basic parameters(num of treat -1)+number of baselines (studies)+add 1 for each 3-arm trial
pD.ipd.prior <- length(unique(rrms_ipd[rrms_ipd$design=="rct",]$study))+length(unique(rrms_ipd$treat))-1+1
pD.ad.prior <- length(unique(rrms_ad[rrms_ad$design=="rct",]$id))+length(unique(rrms_ad$treat))-1

# DIC
dic.ipd.prior <- totresdev.ipd.prior+pD.ipd.prior
dic.ad.prior <- totresdev.ad.prior+pD.ad.prior
# 3. adjust1 
sum_adj1 <- summary(jagsfit_rrms_NMA_dic[[3]])
totresdev.ad.adj1 <- sum_adj1["totresdev.ad","Mean"] # 240.0245
totresdev.ipd.adj1 <- mean(c(jagsfit_rrms_NMA_dic[[3]]$samples[[1]][,"totresdev.ipd"],
                             jagsfit_rrms_NMA_dic[[3]]$samples[[2]][,"totresdev.ipd"])) 
#  90491.8
# Effective number of parameters: IPD and AD: number of basic parameters(num of treat -1)+number of baselines (studies)+add 1 for each 3-arm trial+1 bias parameter (g)+2 bias probabilties (pi) #!+ study-specfic bias indicator (R_j)
pD.ipd.adjust1 <- length(unique(rrms_ipd$study))+length(unique(rrms_ipd$treat))-1+2 +1+2#+length(unique(rrms_ipd$study))
pD.ad.adjust1 <- length(unique(rrms_ad$id))+length(unique(rrms_ad$treat))-1+1+2#+length(unique(rrms_ad$id))

# DIC
dic.ipd.adjust1 <- totresdev.ipd.adj1+pD.ipd.adjust1
dic.ad.adjust1 <- totresdev.ad.adj1+pD.ad.adjust1

# 4. adjust2
sum_adj2 <- summary(jagsfit_rrms_NMA_dic[[4]])
totresdev.ad.adj2 <- sum_adj2["totresdev.ad","Mean"] # 149.8184
totresdev.ipd.adj2 <- mean(c(jagsfit_rrms_NMA_dic[[4]]$samples[[1]][,"totresdev.ipd"],
                             jagsfit_rrms_NMA_dic[[4]]$samples[[2]][,"totresdev.ipd"])) 
#  90353.37
# Effective number of parameters: IPD and AD: number of basic parameters(num of treat -1)+number of baselines (studies)+add 1 for each 3-arm trial+1 bias parameter (g)+2 bias probabilties (pi)
pD.ipd.adjust2 <- length(unique(rrms_ipd$study))+length(unique(rrms_ipd$treat))-1+2+1+2
pD.ad.adjust2 <- length(unique(rrms_ad$id))+length(unique(rrms_ad$treat))-1+1+2+length(unique(rrms_ad$id))

# DIC
dic.ipd.adjust2 <- totresdev.ipd.adj2+pD.ipd.adjust2
dic.ad.adjust2<- totresdev.ad.adj2+pD.ad.adjust2





