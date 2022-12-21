# use d instead of LOR - no need  to all comparisons to be displayed
rrms_ORvsAge <- function(jagsfit_rrms_adjust1_NMR_age, add.cri=T){
  #** 1. create the data to plot
  age <- seq(18,56, l=1000) # the range of age in IPDs
  plotdata <- summary(jagsfit_rrms_adjust1_NMR_age, expo = FALSE) %>% t() %>% data.frame() %>% select(starts_with('d')|starts_with('b'))%>%t()%>% data.frame()
  
  df.plot0 <- list()
  for (i in 1:3) {
    df.plot0[[i]] <- data.frame(age=age,
                                OR=exp(plotdata[i+1,1]+(plotdata[5,1]*(age-38))),
                                low.cri=exp(plotdata[i+1,'X2.5.']+(plotdata[5,'X2.5.']*(age-38))),
                                up.cri=exp(plotdata[i+1,'X97.5.']+(plotdata[5,'X97.5.']*(age-38))),
                                Treatment=c("Dimethyl fumarate",
                                             "Glatiramer acetate",
                                             "Natalizumab")[i]
    )       
  }
  
  df.plot <- do.call(rbind,df.plot0)
  df.plot$Treatment <- factor(df.plot$Treatment, 
                              levels = c("Glatiramer acetate",
                                         "Dimethyl fumarate",
                                         "Natalizumab")
                              ) # to control the order of trt in legend
  
  #** 2. These are color-blind-friendly palettes, one with gray, and one with black.   
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  theme_set(theme_minimal())
  
  
  if(add.cri){ # add credible bands
    #** 3. Appx fig: OR vs age with confidence bands around them 
  ggplot(aes(x = age,y = OR, linetype=Treatment), data = df.plot) +
  geom_smooth(aes(ymin=low.cri, 
                  ymax=up.cri,
                  colour = Treatment, fill=Treatment),
                stat = "identity")+ # this part add the CrI bands
  xlab("Age (in years)")+
    ylab("Odds Ratio (OR)")+
    scale_colour_manual(values=cbbPalette)+ # control the line color
    scale_fill_manual(values=cbbPalette)+  # control the bands color
    theme( # modify text in axis
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.title.x=element_text(size=16,face = "bold"),
          axis.title.y=element_text(size=16,face = "bold"),
          strip.text.x = element_text(size = 14),
          # control legned 
          legend.position = "right",
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14), #change legend text font size
          # add borders to the panel
          panel.border = element_rect(colour = "grey40", fill=NA, size=1.5)
          )+
    coord_cartesian(ylim = c(0.1, 1.8)) # cut the panel, we keep just 0.1 to 1.8 in y-axis
    
  } else { # without credible bands
    #** 4. Fig3: OR vs age without credible bands 
  ggplot(aes(x = age,y = OR,col = Treatment, linetype=Treatment), data = df.plot) +  
    geom_line(size=1.2)+
    xlab("Age (in years)")+
    ylab("Odds Ratio (OR)")+
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
}


