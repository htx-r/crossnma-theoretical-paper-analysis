library(dplyr) # recode
library(readxl)

# read
smsc0 <- read_excel("data/RRMS/raw data/SMSC_new.xlsx")
smsc0 <- as.data.frame(smsc0)

# empty data frame + nr relpases -->  0/1
smsc <- data.frame(relapse01=ifelse(smsc0$nr.relapses.2y.after.treatment>0,1,0))

# change trt names
smsc$trt <- recode(smsc0$ms.specific.treatment,
                    'Copaxone'= 'Glatiramer acetate',
                    'Tecfidera'='Dimethyl fumarate',
                    'Tysabri'='Natalizumab'
                   )

# Women = 0 , Men = 1
smsc$sex <- recode(smsc0$gender,'Men'=1,'Women'=0 )

# age from date of birth
smsc$age <- floor(as.numeric(difftime(Sys.Date(),smsc0$birth.date, units = "weeks"))/52.25)


# change names
nrs.ipd <-with(smsc,data.frame(
  STUDYID='NRS',
  RELAPSE2year=relapse01,
  TRT01A=trt,
  design="nrs",
  AGE=age,
  SEX=sex
))


