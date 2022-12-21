# Compute the data needed to calculate DIC and residual deviance
library(dplyr)
library(tidyr)
library(magrittr)
library(plyr)
# we need to run the model to find this matrix :) 
find_r_n_mat_ipd <- function(mod){
prt.data.ad <- myprt.data%>%arrange(study,treat)%>% group_by(study,treat)%>%
  do(r=sum(.$r),
     n=nrow(.),
     design=unique(.$design)
  )%>%unnest(cols=c(r,n,design))%>%as.data.frame()

#** convert study to the study key as used in JAGS (to make sure they have the same order)
# attach ".ipd" to study names
prt.data.ad$study <- paste0(prt.data.ad$study,".ipd")
prt.data.ad %<>% mutate(study=mapvalues(study,
                                       from=mod$study.key$std.id,
                                       to=mod$study.key$study.jags,
                                       warn_missing = FALSE) %>% as.integer)

prt.data.ad %<>% mutate(treat=mapvalues(treat,
                                        from=mod$trt.key$trt.ini,
                                        to=mod$trt.key$trt.jags,
                                        warn_missing = FALSE) %>% as.integer)

# 2. represent r and n as a matrix with dim: study X treatment arm
r_n_data <- prt.data.ad %>% arrange(study,treat) %>% group_by(study) %>%
dplyr::mutate(arm = row_number()) %>%  ungroup()%>%
  dplyr::select(-c(treat,design))  %>%
  gather("variable", "value",-study, -arm) %>% spread(arm, value)
r.ipd_mat <- as.matrix(r_n_data %>% filter(variable == "r") %>% select(-study, -variable))
n.ipd_mat <- as.matrix(r_n_data %>% filter(variable == "n") %>% select(-study, -variable))
stnd_mat <- log(choose(n.ipd_mat,r.ipd_mat))

# add a column "arm" indicate the participant treatment arm
myprt.data.arrange <- myprt.data%>%arrange(study,treat)
#table(myprt.data.arrange$treat,myprt.data.arrange$study)
arm <- c(rep(1,408),rep(2,826),
         rep(1,363),rep(2,703),rep(3,351),
         rep(1,312),rep(2,627),
         rep(1,138),rep(2,24),rep(3,44))

list(r.ipd=r.ipd_mat, n.ipd=n.ipd_mat,stnd=stnd_mat, arm=arm)
}







