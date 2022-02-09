rrms_beta_dist <- function(ab){
  # data to go to ggplot
  x <- seq(0,1,length=200)
  beta_dist <- data.frame(cbind(x, dbeta(x,1,ab),dbeta(x,ab,1)))
  colnames(beta_dist) <-c("x",paste0("a=1, b=",ab," (low RoB)"),paste0("a=",ab,", b=1 (high RoB)"))
  beta_dist <- reshape2::melt(beta_dist,x)
  colnames(beta_dist)[2] <- "parameters"
  
  # beta distribution plot
  #** hese are color-blind-friendly palettes, one with gray, and one with black.   
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  theme_set(theme_minimal())
 ggplot(beta_dist, aes(x,value, color=parameters))+
   geom_line(size=1.4) + labs(title="Beta Distribution") + 
   labs(x="Bias Probability", y="")+
            theme(
              legend.position = "right",
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
              axis.title.x=element_text(size=16,face = "bold"),
              axis.title.y=element_text(size=16,face = "bold"),
              strip.text.x = element_text(size = 14),
              plot.title = element_text(hjust = 0.5),
              legend.key.size = unit(1, 'cm'), #change legend key size
              legend.key.height = unit(1, 'cm'), #change legend key height
              legend.key.width = unit(1, 'cm'), #change legend key width
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=10), #change legend text font size
              panel.border = element_rect(colour = "grey40", fill=NA, size=1.5))+
            scale_colour_manual(values=cbbPalette)
}
