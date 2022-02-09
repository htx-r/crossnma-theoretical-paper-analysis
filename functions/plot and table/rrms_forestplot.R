rrms_forestplot <- function(jagsfit.nma=jagsfit_rrms_NMA){
  x1 <- data.frame(summary(jagsfit.nma[[1]],expo=F),method='unadjusted') 
  x2 <- data.frame(summary(jagsfit.nma[[2]],expo=F),method='adjust 1') 
  x3 <- data.frame(summary(jagsfit.nma[[3]],expo=F),method='adjust 2') 
  
  x<- rbind(x1,x2,x3)
  drug.lab <- jagsfit.nma[[1]]$trt.key$trt.ini
  
  plotdata <- x %>% t() %>% data.frame() %>% select(starts_with('d'))%>%t()%>% data.frame()
  plotdata$trt <- rep(drug.lab, 3) 
  plotdata <- plotdata[order(plotdata$trt,decreasing = TRUE),] # order them by drug not method
  plotdata <- plotdata[plotdata$trt!="Placebo",] # don't display placebo
  
  metafor::forest(x=as.numeric(as.character(plotdata$Mean)),
                  sei=as.numeric(as.character(plotdata$SD)), 
                  transf=exp,
                  refline=1,
                  xlim=c(-0.1, 1.30), at=c(0.20,.5, 1),
                  #xlim=c(-0.1, 5), at=c(0.5, 1,2,4),
                  slab=plotdata$method,
                  cex=1.5,
                  ylim=c(1,20),
                  rows=c(15:13,10:8,5:3), # the position of each CrI !! should start from up (high num) to down
                  psize=1, # the size of the circle of point estimate 
                  header=c("Treatment and method of synthesis","Estimate OR [95% CrI]"),#,,
                  cex.axis=1.5,
                  xlab="Odds Ratio (OR)")
  
  # add treatment names 
  text(x=-0.1, 
       y=c(6,11,16), 
       pos=4, 
       font = 2,
       labels=c("Dimethyl fumarate",
                "Glatiramer acetate",
                "Natalizumab"),
       cex=1.5)
  }


