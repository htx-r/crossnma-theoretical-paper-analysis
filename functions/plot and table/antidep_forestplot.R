library(crossnma)
library(ggplot2)
antidep_forestplot <- function(jagsfit_antidep_main, antidep){
  jagsfit1 <- jagsfit_antidep_main[[1]]
  jagsfit2 <- jagsfit_antidep_main[[2]]
  jagsfit3 <- jagsfit_antidep_main[[3]]

  # Data to plot
  
  # unadjust 
  sum1 <-summary(jagsfit1,expo=T)
  sum1_1 <- as.data.frame(sum1[-c(1,nrow(sum1)),])
  sum1_1$drug <- sort(unique(as.character(antidep$agent)))[-16]
  sum1_1$method <- "unadjusted"
  
  # adjust1
  sum2 <-summary(jagsfit2,expo=T)
  sum2_1 <- as.data.frame(sum2[2:22,])
  sum2_1$drug <- sort(unique(as.character(antidep$agent)))[-16]
  sum2_1$method <- "bias-adjusted 1"
  
  # adjust2
  sum3 <-summary(jagsfit3,expo=T)
  sum3_1 <- as.data.frame(sum3[2:22,])
  sum3_1$drug <- sort(unique(as.character(antidep$agent)))[-16]
  sum3_1$method <- "bias-adjusted 2"
  
  # gather them in one dataframe
  plotdata <- rbind.data.frame(sum1_1,sum2_1,sum3_1)
  plotdata$method <- factor(plotdata$method, levels = c("bias-adjusted 2","bias-adjusted 1","unadjusted"))
  
  # 
  theme_set(
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(face='bold',size=14),
          axis.title.x=element_text(size=16,face = "bold"),axis.title.y=element_text(size=16,face = "bold"),
          strip.text.x = element_text(size = 16),
          panel.spacing.x = unit(4, "mm")
    )
  )
  g_split <- ggplot(plotdata, aes(y=Mean, x=method,shape=method)) +
    geom_point()+ # point estimate
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), 
                  width=0.4,size=0.5) +
    coord_flip()+
    facet_wrap(~drug)
  
  
  g <- g_split + 
    xlab("") +
    ylab("Odds Ratio (OR)")
  
  g 
}

