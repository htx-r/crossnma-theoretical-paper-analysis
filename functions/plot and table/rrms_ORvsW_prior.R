rrms_ORvsW_prior <- function(jagsfit.prior){
  #** data to plot: OR at each value of w
  w <- seq(0,1,l=11)[-1]
  df.plot0 <- list()
  for (i in 1:length(w)) {
    plotdata <- summary(jagsfit.prior[[i]], expo = FALSE) %>% t() %>% data.frame() %>% select(starts_with('LOR'))%>%t()%>% data.frame()
    df.plot0[[i]] <- data.frame(w=w[i],
                                OR=exp(plotdata[,1]),
                                # X2.5.=exp(plotdata[,'X2.5.']),
                                # X97.5.=exp(plotdata[,'X97.5.']),
                                Contrast=c("Dimethyl fumarate vs Natalizumab",
                                           "Glatiramer acetate vs Natalizumab",
                                           "Placebo vs Natalizumab",
                                           "Glatiramer acetate vs Dimethyl fumarate",
                                           "Placebo vs Dimethyl fumarate",
                                           "Placebo vs Glatiramer acetate")
    )
  }
  
  df.plot <- do.call(rbind,df.plot0)
  
  df.plot <- df.plot[df.plot$Contrast%in%c("Dimethyl fumarate vs Natalizumab",
                                           "Glatiramer acetate vs Dimethyl fumarate",
                                           "Glatiramer acetate vs Natalizumab"),]
  df.plot$Contrast <- factor(df.plot$Contrast, levels =c("Glatiramer acetate vs Natalizumab",
                                                         "Dimethyl fumarate vs Natalizumab",
                                                         "Glatiramer acetate vs Dimethyl fumarate"
                                                         ) 
                             )
  #** 2. These are color-blind-friendly palettes, one with gray, and one with black.   
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #** 3. plot OR vs w
  theme_set(theme_minimal())
  ggplot(data=df.plot, aes(x=w, y=OR, col=Contrast,shape=Contrast)) +
    geom_line(size=1.2)+
    geom_point(size=3)+
    ylab("")+
    xlab("")+
    theme(
      legend.position = "right",
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16),
      axis.title.x=element_text(size=16,face = "bold"),
      axis.title.y=element_text(size=16,face = "bold"),
      strip.text.x = element_text(size = 14),
      legend.key.size = unit(1, 'cm'), #change legend key size
      legend.key.height = unit(1, 'cm'), #change legend key height
      legend.key.width = unit(1, 'cm'), #change legend key width
      legend.title = element_text(size=14), #change legend title font size
      legend.text = element_text(size=10), #change legend text font size
      panel.border = element_rect(colour = "grey40", fill=NA, size=1.5))+
    scale_colour_manual(values=cbbPalette)
  }
