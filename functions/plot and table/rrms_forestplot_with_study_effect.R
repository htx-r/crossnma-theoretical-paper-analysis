#jagsfit.nma <- jagsfit_rrms_NMA_dic
source("run_analysis/RRMS/R1_estimate_per_study.R")
rrms_forestplot_with_study_effect <- function(jagsfit.nma){
  x1 <- data.frame(summary(jagsfit.nma[[1]],exp=F),method='unadjusted') 
  x2 <- data.frame(summary(jagsfit.nma[[3]],exp=F),method='bias-adjusted 1') 
  x3 <- data.frame(summary(jagsfit.nma[[4]],exp=F),method='bias-adjusted 2') 
  
  x<- rbind(x1,x2,x3)
  drug.lab <- jagsfit.nma[[1]]$trt.key$trt.ini
  #-- for comparison of Natalizumab against other drugs
  plotdata1 <- x %>% t()%>%data.frame()  %>% select(starts_with(c(# NZ vs each of the 3 drugs
                                                                  # NZ vs Pl
                                                                  #"theta.3.4.1", # AFFIRM
                                                                  "LOR.4.1",
                                                                  # NZ vs DF
                                                                  #"theta.4.4.1", # NRS
                                                                  "LOR.4.2",
                                                                  # NZ vs GA
                                                                  #"theta.4.2.1", # NRS
                                                                  "LOR.4.3",
                                                                  # GA vs PL and vs DF
                                                                  # "theta.2.3.1", # CONFIRM
                                                                  # "theta.5.3.1", # Born
                                                                  # "theta.6.3.1", # Joh
                                                                  "LOR.3.1",
                                                                  #  DF vs GA
                                                                  # "theta.2.2.1", # CONFIRM
                                                                  # "theta.4.2.2", # Swiss cohort
                                                                  "LOR.3.2",
                                                                  # DF vs PL 
                                                                  # "theta.1.2.1", # DEFINE
                                                                  # "theta.2.3.1" , # CONFIRM
                                                                  "LOR.2.1"
                                                                  )))%>%t()%>% data.frame()
  plotdata11 <- type.convert(plotdata1, as.is = TRUE) # choose appropriate class for columns  
  plotdata11 <- plotdata11[,c("Mean","SD","method")]
  plotdata11[c("LOR.3.2.","LOR.3.2.1","LOR.3.2.2"),"Mean"] <- -plotdata11[c("LOR.3.2.","LOR.3.2.1","LOR.3.2.2"),"Mean"]
  plotdata11$comp <- rep(c("NZ vs PL",
                           "NZ vs DF",
                           "NZ vs GA",
                           "GA vs PL",
                           "DF vs GA",
                           "DF vs PL"),each=3)
  rownames(plotdata11) <- NULL
  # bind the crossNMA estimates with study-specific estimates (calculated in different file)
  plotdata <- rbind(plotdata11,study_specific_est)
  plotdata <- plotdata%>% arrange(comp,method)
  plotdata_final <- plotdata[c(3,4,5,1,2,
             8,9,10,6,7,
             13,14,15,16,11,12,
             19,20,17,18,
             23,24,21,22,
             25,28,26,27
             ),]

  

  
  #--- Plot
  metafor::forest(x=as.numeric(as.character(plotdata_final$Mean)),
                  sei=as.numeric(as.character(plotdata_final$SD)), 
                  transf=exp,
                  refline=1,
                  alim=c(-1,2),
                  xlim=c(-1, 2.5), 
                  at=c(0,0.20,.5, 1,1.5,2),
                  slab=plotdata_final$method,
                  #pch=rep(c(22,23,23,23),3),
                  #col=c("green",rep("black",3),"red",rep("black",3), "red",rep("black",3)),
                  ylim=c(1,210),
                  rows=c(196,192,184,180,176,
                         160,156,148,144,140,
                         122,118,114,106,102,98,# GA vs PL
                         82,74,70,66, # NZ vs DF
                         50,42,38,34, # NZ vs GA
                         18,10,6,2 # NZ vs PL
                         ), # the position of each CrI !! should start from up (high num) to down
                  #rows=rev(seq(1,82,by=3)),
                  psize=1, # the size of the circle of point estimate 
                  header=c("Treatment and method of synthesis","OR [95% CrI]"),#,,
                  cex.axis=1.2,
                  cex=0.8,
                  cex.lab=1,
                  xlab="Odds Ratio (OR)")
  
  # add treatment names 
  text(x=-1,#x=rep(c(-1,-0.9,-0.9),6), 
       y=c(204,200,188,
           168,164,152,
           130,126,110,
           90,86,78,
           58,54,46,
           26,22,14), 
       pos=4, 
       font = 2,
       labels=c("Dimethyl fumarate vs Placebo",
                "Study-specific estimates",
                "NMA estimate",
                
                "Dimethyl fumarate vs Glatiramer acetate",
                "Study-specific estimates",
                "NMA estimate",
                
                "Glatiramer acetate vs Placebo",
                "Study-specific estimates",
                "NMA estimate",
                
                "Natalizumab vs Dimethyl fumarate",
                "Study-specific estimates",
                "NMA estimate",
                
                "Natalizumab vs Glatiramer acetate",
                "Study-specific estimates",
                "NMA estimate",
                
                "Natalizumab vs Placebo",
                "Study-specific estimates",
                "NMA estimate"
                ),
       cex=0.85)
}









# NZ vs GA = theta[NZ, DF]-theta[GA,DF] #! (check that with Georgia)
# NZ vs others
# plotdata11["theta.4.2.2",] <- plotdata11["theta.4.2.1",]  
# plotdata11["theta.4.2.1","Mean"] <- plotdata11["theta.4.4.1","Mean"]- plotdata11["theta.4.2.1","Mean"]
# plotdata11["theta.4.2.1","SD"] <- sqrt(plotdata11["theta.4.4.1","SD"]^2+ plotdata11["theta.4.2.1","SD"]^2)
# plotdata11["theta.3.4.1","method"] <- "AFFIRM"
# plotdata11["theta.4.4.1","method"] <- "Swiss Cohort"
# plotdata11["theta.4.2.1","method"] <- "Swiss Cohort"
# 
# # GA vs others
# plotdata11["theta.2.2.1","Mean"] <- plotdata11["theta.2.3.1","Mean"]- plotdata11["theta.2.2.1","Mean"]
# plotdata11["theta.2.2.1","SD"] <- sqrt(plotdata11["theta.2.3.1","SD"]^2+ plotdata11["theta.2.2.1","SD"]^2)
# plotdata11["theta.2.3.1","method"] <- "CONFIRM"
# plotdata11["theta.5.3.1","method"] <- "Bornstein"
# plotdata11["theta.6.3.1","method"] <- "Johnson"
# plotdata11["theta.4.2.2","method"] <- "Swiss Cohort"
# plotdata11["theta.2.2.1","method"] <- "CONFIRM"
# # GA vs others
# plotdata11["theta.1.2.1","method"] <- "DEFINE"
# rrms_forestplot_with_study_effect_GA_DF <- function(jagsfit.nma){
#   x1 <- data.frame(summary(jagsfit.nma[[1]],expo=F),method='unadjusted') 
#   x2 <- data.frame(summary(jagsfit.nma[[3]],expo=F),method='bias-adjusted 1') 
#   x3 <- data.frame(summary(jagsfit.nma[[4]],expo=F),method='bias-adjusted 2') 
#   
#   x<- rbind(x1,x2,x3)
#   drug.lab <- jagsfit.nma[[1]]$trt.key$trt.ini
#     #-- for comparison of Glatiramer acetate against other drugs
#   plotdata1 <- x %>% t()%>%data.frame()  %>% select(starts_with(c("theta.2.3.1",
#                                                                   "theta.5.3.1",
#                                                                   "theta.6.3.1",
#                                                                   "LOR.3.1",
#                                                                   "theta.2.2.1",
#                                                                   "theta.4.2.1",
#                                                                   "LOR.3.2")))%>%t()%>% data.frame()
#   plotdata11 <- type.convert(plotdata1, as.is = TRUE) # choose appropriate class for columns  
#   
#   
#   plotdata11["theta.2.2.1","Mean"] <- plotdata11["theta.2.3.1","Mean"]- plotdata11["theta.2.2.1","Mean"]
#   plotdata11["theta.2.2.1","SD"] <- sqrt(plotdata11["theta.2.3.1","SD"]^2+ plotdata11["theta.2.2.1","SD"]^2)
#   plotdata11["theta.2.3.1","method"] <- "CONFIRM"
#   plotdata11["theta.5.3.1","method"] <- "Bornstein"
#   plotdata11["theta.6.3.1","method"] <- "Johnson"
#   metafor::forest(x=as.numeric(as.character(plotdata11$Mean)),
#                   sei=as.numeric(as.character(plotdata11$SD)), 
#                   transf=exp,
#                   refline=1,
#                   xlim=c(-1, 4), at=c(0.20,.5, 1,1.5,2,2.5,3),
#                   #xlim=c(-0.1, 5), at=c(0.5, 1,2,4),
#                   slab=plotdata11$method,
#                   pch=rep(c(22,23,23,23),3),
#                   cex=1.5,
#                   col=c("green",rep("black",3),"red",rep("black",3), "red",rep("black",3)),
#                   ylim=c(1,25),
#                   rows=c(20,18:16,14,12:10,7,5:3), # the position of each CrI !! should start from up (high num) to down
#                   psize=1, # the size of the circle of point estimate 
#                   header=c("Treatment and method of synthesis","OR [95% CrI]"),#,,
#                   cex.axis=1.5,
#                   xlab="Odds Ratio (OR)")
#   
#   # add treatment names 
#   text(x=-1, 
#        y=c(21,19,15,13,8,6), 
#        pos=4, 
#        font = 2,
#        labels=c("Natalizumab vs Placebo",
#                 "Summary measure",
#                 "Natalizumab vs Dimethyl fumarate",
#                 "Summary measure",
#                 "Natalizumab vs Glatiramer acetate",
#                 "Summary measure"
#        ),
#        cex=1.5)
# }

# trt.ini trt.jags
#             Placebo        1
#   Dimethyl fumarate        2
#  Glatiramer acetate        3
#         Natalizumab        4
# std.id study.jags
# 1   DEFINE.ipd          1
# 2  CONFIRM.ipd          2
# 3   AFFIRM.ipd          3
# 4      NRS.ipd          4
# 5 Bornstein.ad          5
# 6   Johnson.ad          6


