# ------------ # prepare RCT-AD # ------------ #
require(readxl)
# read
MSaggregate <- read_excel("data/RRMS/raw data/MSaggregate.xlsx")
MSaggregate <- as.data.frame(MSaggregate)

# exclude IPDs
toexclude.IPDstudies <- c("MSCRG","AFFIRM","DEFINE","CONFIRM")
rct.ad0 <- MSaggregate[!MSaggregate$study%in%toexclude.IPDstudies,]
# only Johnson and Borenstien
rct.ad <- rct.ad0[rct.ad0$study%in%c("Johnson","Bornstein"),]
