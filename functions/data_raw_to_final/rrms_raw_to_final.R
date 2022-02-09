
library(magrittr)
library(plyr)
# raw DATA
source("functions/data_raw_to_final/rrms_ad.R") # rct.ad
source("functions/data_raw_to_final/rrms_ipd.R") # rct.ipd
source("functions/data_raw_to_final/rrms_nrs.R") # nrs.ipd: wait to the new data
# combine the two designs: NRS and RCT

# IPD
prt.data <- rbind.data.frame(rct.ipd[,c("STUDYID",
                                        "RELAPSE2year",
                                        "TRT01A",
                                        "design",
                                        "AGE",
                                        "SEX")], nrs.ipd)


prt.data$TRT01A <- as.character(prt.data$TRT01A)

#prt.data$STUDYID <- as.character(prt.data$study)
prt.data%<>%mutate(bias=mapvalues(STUDYID,
                                  from=unique(STUDYID),
                                  to=c('low','low','low','high')))
prt.data$unfav <- prt.data$bias.group <- NA
study_plac <- unique(prt.data$STUDYID[prt.data$TRT01A=="Placebo"]) # index of studies with placebo
for(i in 1:length(study_plac)){
  prt.data[prt.data$STUDYID==study_plac[i]&prt.data$TRT01A=="Placebo",]["unfav"] <- 0
  prt.data[prt.data$STUDYID==study_plac[i]&prt.data$TRT01A!="Placebo",]["unfav"]  <- 1
  prt.data[prt.data$STUDYID==study_plac[i],]["bias.group"]<- 1
}

# NRS
prt.data[prt.data$STUDYID=="NRS"&prt.data$TRT01A=="Glatiramer acetate",]["unfav"] <- 0
prt.data[prt.data$STUDYID=="NRS"&prt.data$TRT01A!="Glatiramer acetate",]["unfav"]  <- 1
prt.data[prt.data$STUDYID=="NRS",]["bias.group"]<- 1
prt.data <- prt.data[,c("STUDYID","RELAPSE2year","TRT01A","design",  "AGE","bias","unfav","bias.group")]
names(prt.data) <- c("study","r","treat","design","age","bias","unfav","bias.group")

# AD
rct.ad$design <- 'rct'
#rct.ad$cov <-c(34.3,34.3,30,30)
rct.ad$bias <- c('high','high','high','high')
rct.ad$sex <- c(0.2,0.3,0.4,0.5)
rct.ad$age <- c(34.3,34.3,30,30)
rct.ad$edss <- c(2,3,5,6)

rct.ad$unfav <- rct.ad$bias.group <- NA

for(i in 1:length(unique(rct.ad$study))){
  rct.ad[rct.ad$study==unique(rct.ad$study)[i]&rct.ad$treat=="Placebo",]["unfav"] <- 0
  rct.ad[rct.ad$study==unique(rct.ad$study)[i]&rct.ad$treat!="Placebo",]["unfav"]  <- 1
  rct.ad[rct.ad$study==unique(rct.ad$study)[i],]["bias.group"]<- 1
}

# Save files
write.csv(prt.data, "data/RRMS/final data/rrms_final_IPD")
write.csv(rct.ad, "data/RRMS/final data/rrms_final_AD")

