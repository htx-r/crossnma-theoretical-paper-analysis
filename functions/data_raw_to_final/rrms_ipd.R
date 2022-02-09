#!!! error when I run the clean function, I need to import dataset from the right panel
        # but without actually execute in the end
# ------------ # prepare RCT-IPD (BIOGEN) # ------------ #
source("functions/data_raw_to_final/rrms_numericalDataRisk.R")
# load libraries
library(devtools)
install_github("htx-r/CleaningData",force=TRUE)
library(CleaningData)
library(dplyr) # for numericalDataRisk.fun: recode
library(vcd) # for numericalDataRisk.fun: assocstats

# your data path
datapath="data/RRMS/raw data"

# cleaned IPD(BIOGEN): add RELAPSE02Year, RELAPSE01Year and re-name studies and drugs
cleanBIOGENtrials<-cleanBIOGENtrials.fun(datapath=datapath)
RCTs0<-cleanBIOGENtrials$adsl01

# exclude studies
  # SENTTINEL - has combination of treatments and
  # ADVANCE - no data in 2-year relapse (only 1-year)
  # MSCRG -
RCTs<-RCTs0[!RCTs0$STUDYID%in% c("SENTINEL","ADVANCE","MSCRG"),]

# manage variables:
 # exclude variable with more than 30% NA's
 # exclude factors with just one category (as SAFFL=Safety population flag)
 # exclude factors that are transformations from already existing variables (e.g. TRT01AN is the same as TRTO1A but numeric)
 # exclude highly correlated variables
 # recode factors as numerical values (e.g. Male=1, Female=0)
 # transformations of continuous variables to approximate normal distribution
rct.ipd <-numericalDataRisk.fun(RCTs)  ##final full dataset
rct.ipd$design <- 'rct'
rct.ipd$RELAPSE2year <- as.numeric(rct.ipd$RELAPSE2year)-1



