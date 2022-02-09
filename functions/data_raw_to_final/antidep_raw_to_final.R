library(readxl)
library(dplyr)
gris0<- read_excel("data/antidepressant/raw data/GRISdataRESPmergeALL.xlsx")
gris0 <- as.data.frame(gris0)

#** select columns and rename them **#
gris <- gris0 %>% select(studyID, drug_name,
                         Responders,Ntotal,
                         Sponsoredbythisdrugcompanyo,
                         OverallRoB)
colnames(gris) <- c("study","drug", "r","n","spon","rob")
gris$design <- "rct"
gris$unfav <- NA
gris$bias.group <- NA

# change names of RoB
gris$rob[gris$rob=="Low risk"] <- "low"
gris$rob[gris$rob=="Moderate risk"] <- "high"
std.na <- unique(gris$study[is.na(gris$rob)])
gris <- gris[!gris$study%in%std.na,]
#** Bias assumption for studies without placebo

# studies without placebo have different sponsorship forms:
 # 1. 1 No and 1 Yes (most cases)
 # 2. 2 No (both are not sponsored)
 # 3. 1 Unclear 1 Yes 
 # 4. 2 or 3 Unclear
 # 5. 1 No 1 Unclear
 # 6. 2 Yes (only one study)
 # 7. 1 No and 2 Yes 
table(gris$spon)

study_plac <- gris$study[gris$drug=="Plac"] # index of studies with placebo
study_no_plac <- unique(gris$study)[!unique(gris$study)%in%study_plac]
tab.unfav <- list()
for (i in study_no_plac) {
  tab.unfav[[which(study_no_plac==i)]] <- table(gris$spon[gris$study==i])
}
tab.unfav

# I will consider the following in my analysis:
 # Unclear = Yes (consider unclear trt as sponsored)
 # studies with only Yes's or with onyl No's --> no need for bias adjustment
gris$spon[gris$spon=="Unclear"] <- "Yes"


#** add two columns to the dataset **#
# 1. unfav: 0 (treatment is unfavoured) & 1 (treatment is favoured)
# 2. bias.group: 1 (study has inactive treatment and need bias adjustment)
             # 0 (study has only active treatments and no need for bias adjustment)     
             # 2 (study has only active treatments and need another bias adjustment) 

# studies with placebo: 0 to placebo and 1 to non-palcebo

study_plac <- gris$study[gris$drug=="Plac"] # index of studies with placebo

for(i in study_plac){
  gris[gris$study==i&gris$drug=="Plac",]["unfav"] <- 0
  gris[gris$study==i&gris$drug!="Plac",]["unfav"] <- 1
  gris[gris$study==i,]["bias.group"]<- 1
  }

# studies without placebo
#**# So these are the cases that I will be analysing
# 1.+ 5. + 7.: "No" will be set as unfavoured, 2=bias.group
# 2.+ 3. + 4. + 6.: the first will be set as unfavoured but no bias adjustment will be applied

#  0 to not sponsored and 1 otherwise

# index of studies without placebo and has a unique not sponsored treatment
study_no_plac <- unique(gris$study)[!unique(gris$study)%in%study_plac]
study_no_plac_unsp0 <- c()
for (i in study_no_plac) {
  study_no_plac_unsp0[which(study_no_plac==i)] <- ifelse(length(table(gris$spon[gris$study==i]))==1, FALSE, TRUE)
}
study_no_plac_unsp <- study_no_plac[study_no_plac_unsp0]
for(i in study_no_plac_unsp){
  gris[gris$study==i&gris$spon=="No",]["unfav"] <- 0
  gris[gris$study==i&gris$spon!="No",]["unfav"] <- 1
  gris[gris$study==i,]["bias.group"]<- 2
}

# studies without placebo and all either sponsored or not sponsored: 1 to all except first treatment 0 and no bias adjustment 
study_no_plac_all_sp_or_unsp <-study_no_plac[!study_no_plac_unsp0]

for(i in study_no_plac_all_sp_or_unsp){
  gris[gris$study==i,]["unfav"] <- 1
  gris[gris$study==i,][1,]["unfav"] <- 0
  gris[gris$study==i,]["bias.group"]<- 0
}

# check which study has 2 zeros (No's) 
ns <- length(unique(gris$study))
std.0unfav0 <- c()

for (i in 1:ns) {
  std.0unfav0[i] <- nrow(gris[(gris$study==unique(gris$study)[i])&gris$unfav==0,])!=1
}

std.0unfav <- unique(gris$study)[std.0unfav0] # 202
gris[gris$study%in%std.0unfav,]

#!! we have one study with 2 not sponsored treatments and one sponsored
# give the first No a zero and the other 1
gris[gris$study%in%std.0unfav & gris$unfav==0,][2,]["unfav"] <- 1

#** End: save the dataset that will be used in the analysis
write.csv(gris, file = 'data/antidepressant/final data/antidepressant_for_main_analysis')


