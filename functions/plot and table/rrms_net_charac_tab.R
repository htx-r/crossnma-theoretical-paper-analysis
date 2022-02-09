rrms_net_charac_tab <- function(myprt.data,mystd.data){
  # convert IPD to arm-level data and combine it with AD 
  prt.data.ad0 <- sapply(1:length(unique(myprt.data$study)), 
                         function(i){
                           dd.stdy <- subset(myprt.data, myprt.data$study==unique(as.character(myprt.data$study))[i])
                           dd.mat <- as.matrix(table(as.character(dd.stdy$treat),dd.stdy$r))
                           # subset of study i
                           data.frame(
                             study=2+i, 
                             treat=rownames(dd.mat),
                             r=dd.mat[,"1"],
                             n =rowSums(dd.mat),
                             design=unique(as.character(dd.stdy$design))
                           )
                         }, simplify = F)
  prt.data.ad <- do.call(rbind,prt.data.ad0)
  mystd.data$study <- c(1,1,2,2)
  all.data <- rbind(mystd.data[c('study','treat','r','n','design')],prt.data.ad)
  
  
  dich.slr <- data.prep(arm.data = all.data,
                        varname.t = "treat",
                        varname.s = "study")
  
  network.char <- net.tab(data = dich.slr,
                          outcome = "r",
                          N = "n", 
                          type.outcome = "binomial",
                          time = NULL)
  network.char
}