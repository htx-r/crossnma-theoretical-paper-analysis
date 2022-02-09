# -------- load the libray --------#
devtools::install_github("htx-r/crossnma",force = TRUE)
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

# Figure 1: network plot

# a. RRMS
crossnma::netplot(jagsfit_rrms_NMA[[1]]$model,
        #seq = "optimal", start = "prcomp",
        cex=1.25,
        plastic = FALSE, thickness = "se.f", number = FALSE,
        points = TRUE, cex.points = 7, col.points = "darkred",
        multiarm = FALSE,
        col='grey',
        bg.number.of.studies =NULL,
        col.number.of.studies = 'black')
# b. antidepressant
netplot(jagsfit_antidep_main[[1]]$model, 
        cex=1.25,
        plastic = FALSE, thickness = "se.f", number = FALSE,
        points = TRUE, cex.points = 7, col.points = "darkred",
        multiarm = FALSE,
        col='grey',
        bg.number.of.studies =NULL,
        col.number.of.studies = 'black'
        )

# Figure 2: forest plot of OR (PL vs NZ, GA, DF) using 3 different methods of estimation
source("functions/plot and table/rrms_forestplot.R")
rrms_forestplot(jagsfit.nma=jagsfit_rrms_NMA)

# Figure 3
source("functions/plot and table/rrms_ORvsAge.R")
rrms_ORvsAge(jagsfit_rrms_adjust1_NMR_age, add.cri=F)

# Figure 4
source("functions/plot and table/antidep_forestplot.R")
antidep <- read.csv("data/antidepressant/raw data/antidepressant")# to get drug names from the raw data
antidep_forestplot(jagsfit_antidep_main, antidep)

# Table 3
#-------- data --------# (needed to create Appendix table 1,2,3)
myprt.data <- read.csv("data/RRMS/final data/rrms_final_IPD") 
mystd.data <- read.csv("data/RRMS/final data/rrms_final_AD") 




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


