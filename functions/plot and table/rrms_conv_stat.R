rrms_conv_stat_tab <- function(mod, method){
  x0 <- data.frame(summary(mod,expo=F, digits = 4),method=method) 
  x <- x0 %>% t() %>% data.frame() %>% select(!starts_with('LOR'))%>%t()%>% data.frame()
  x[,c("Rhat","n.eff","method")]
}
