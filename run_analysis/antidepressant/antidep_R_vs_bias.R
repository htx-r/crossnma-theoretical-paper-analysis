#  Cross-Network meta-analysis of antidepressants 
# ** Main analysis - monitor R and compare estimated R to **
#  Authors: Tasnim Hamza and Georgia Salanti

# This file use antidepressnats to run these analysis
# 1. Unadjusted                                  "Main analysis"
# 2. bias-adjust1 with Beta(1,20) and Beta(20,1) to model bias probability "Main analysis"
# 3. bias-adjust2 with Beta(1,20) and Beta(20,1) to model bias probability "Main analysis"

#-------- load the libray --------#
#install.packages("crossnma")
library(crossnma)

#-------- data --------#
gris <- read.csv("data/antidepressant/final data/antidepressant_for_main_analysis")

#-------- MCMC settings --------#
n.adapt = 2000
n.iter=100000
n.burnin = 40000
thin=1
n.chains=2

#** 2. bias-adjust1 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust1_main <- crossnma.model(std.data=gris,
                                   trt=drug,
                                   study=study,
                                   outcome=r,
                                   n=n,
                                   design=design,
                                   reference="Plac",
                                   trt.effect="random",
                                   #---------- bias adjustment ----------
                                   method.bias='adjust1',
                                   bias=rob,
                                   unfav=unfav,
                                   bias.group = bias.group,
                                   bias.type='add',
                                   bias.effect='random',
                                   prior =list(pi.high.rct="dbeta(20,1)",
                                               pi.low.rct="dbeta(1,20)"))

#mod_adjust1_main$jagsmodel <- antidep_jags_adj1_dic
# run jags
fit_adjust1_main_R <- crossnma(mod_adjust1_main,
                               n.adapt = n.adapt,
                               n.iter=n.iter,
                               n.burnin = n.burnin,
                               thin=thin,
                               n.chains=n.chains,
                               monitor=c("g.act","R"))

#** 3. bias-adjust2 with Beta(1,20) and Beta(20,1)
# jags model: code+data
mod_adjust2_main <- crossnma.model(std.data=gris,
                                   trt=drug,
                                   study=study,
                                   outcome=r,
                                   n=n,
                                   design=design,
                                   reference="Plac",
                                   trt.effect="random",
                                   #---------- bias adjustment ----------
                                   method.bias='adjust2',
                                   unfav=unfav,
                                   bias = rob,
                                   bias.group = bias.group,
                                   bias.type='add',
                                   bias.effect='random',
                                   prior =list(pi.high.rct="dbeta(20,1)",
                                               pi.low.rct="dbeta(1,20)")
)

#mod_adjust2_main$jagsmodel <- antidep_jags_adj2_dic
# run jags
fit_adjust2_main_pi <- crossnma(mod_adjust2_main,
                                n.adapt = n.adapt,
                                n.iter=n.iter,
                                n.burnin = n.burnin,
                                thin=thin,
                                n.chains=n.chains,
                                monitor=c("g.act","pi"))


# Reviewer 2 suggested to check the agreement between ( Will the study that is deemed at high/low be assumed as so in the analysis)
  # our judgment about the RoB (high/low) and 
  # the estimated bias indicator (R) (adjust1) or the bias probablity (pi) (adjust2)

# adjust1
sum_adjust1 <- summary(fit_adjust1_main_R,exp=FALSE)
R_study <- sum_adjust1[startsWith(rownames(sum_adjust1),"R"),"Mean"]
R_bias_mat <- cbind(round(R_study),mod_adjust1_main$data$bias_index )

std_high <- R_bias_mat[R_bias_mat[,2]==1,]
nrow(std_high) # Num of high RoB
nrow(std_high[std_high[,1]==1,]) # all high RoB studies are given R=1 
std_low <- R_bias_mat[R_bias_mat[,2]==2,]
nrow(std_low) # Num of low RoB
nrow(std_low[std_low[,1]==0,]) # all low RoB studies are given R=0 

# adjust2
sum_adjust2 <- summary(fit_adjust2_main_pi,exp=FALSE)
pi_designRoB <- sum_adjust2[startsWith(rownames(sum_adjust2),"pi"),"Mean"]
# pi[1]= 0.94306678      
# pi[2]=0.09472230

#== SAVE output
jagsfit_antidep_main_R_pi <- list(fit_adjust1_main_R,fit_adjust2_main_pi)
save(jagsfit_antidep_main_R_pi,file = "output/antidepressant/JAGS/jagsfit_antidep_main_R_pi.RData")
