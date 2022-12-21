#-- Use glm to compute the estimates from IPDs
rrms_ipd <- read.csv("data/RRMS/final data/rrms_final_IPD") 
rrms_ad <- read.csv("data/RRMS/final data/rrms_final_AD") 
rrms_ipd0 <- within(rrms_ipd, treat <- relevel(treat, ref = "Placebo"))

# DEFINE: DF vs PL
mod_std1 <- glm(r~treat,rrms_ipd0[rrms_ipd0$study=="DEFINE",],family = binomial(link='logit'))
exp(coef(mod_std1))
exp(confint(mod_std1))
study_specific_est <- data.frame(Mean=summary(mod_std1)$coefficients[2,"Estimate"], 
                  SD=summary(mod_std1)$coefficients[2,"Std. Error"],
                  method="DEFINE",
                  comp="DF vs PL")
# CONFIRM: GA vs PL, DF vs PL and GA vs DF (or vice versa)
mod_std2 <- glm(r~treat,rrms_ipd0[rrms_ipd0$study=="CONFIRM",],family = binomial(link='logit'))
exp(coef(mod_std2))
exp(confint(mod_std2))
study_specific_est <- rbind(study_specific_est,data.frame(Mean=summary(mod_std2)$coefficients[-1,"Estimate"], 
                                                          SD=summary(mod_std2)$coefficients[-1,"Std. Error"],
                                                          method="CONFIRM",
                                                          comp=c("DF vs PL","GA vs PL")))
# AFFIRM: NZ vs Pl
mod_std3 <- glm(r~treat,rrms_ipd0[rrms_ipd0$study=="AFFIRM",],family = binomial(link='logit'))
exp(coef(mod_std3))
exp(confint(mod_std3))
study_specific_est <- rbind(study_specific_est,data.frame(Mean=summary(mod_std3)$coefficients[-1,"Estimate"], 
                                                          SD=summary(mod_std3)$coefficients[-1,"Std. Error"],
                                                          method="AFFIRM",
                                                          comp="NZ vs PL"))

# NRS: NZ vs DF, NZ vs GA, GA vs DF (or vice versa)
rrms_ipd0 <- within(rrms_ipd, treat <- relevel(treat, ref = "Dimethyl fumarate"))
mod_std4 <- glm(r~treat,rrms_ipd0[rrms_ipd0$study=="NRS",],family = binomial(link='logit'))
exp(coef(mod_std4))
exp(confint(mod_std4))
study_specific_est <- rbind(study_specific_est,data.frame(Mean=summary(mod_std4)$coefficients[3,"Estimate"], 
                                                          SD=summary(mod_std4)$coefficients[3,"Std. Error"],
                                                          method="Swiss Cohort",
                                                          comp=c("NZ vs DF")))
# NRS: NZ vs GA, DF vs GA
rrms_ipd0 <- within(rrms_ipd, treat <- relevel(treat, ref = "Glatiramer acetate"))
mod_std44 <- glm(r~treat,rrms_ipd0[rrms_ipd0$study=="NRS",],family = binomial(link='logit'))
exp(coef(mod_std44))
exp(confint(mod_std44))
study_specific_est <- rbind(study_specific_est,data.frame(Mean=summary(mod_std44)$coefficients[-1,"Estimate"], 
                                                          SD=summary(mod_std44)$coefficients[-1,"Std. Error"],
                                                          method="Swiss Cohort",
                                                          comp=c("NZ vs GA","DF vs GA")))

# CONFIRM: DF vs GA (or vice versa)
rrms_ipd0 <- within(rrms_ipd, treat <- relevel(treat, ref = "Glatiramer acetate"))
mod_std5 <- glm(r~treat,rrms_ipd0[rrms_ipd0$study=="CONFIRM",],family = binomial(link='logit'))
exp(coef(mod_std5))
exp(confint(mod_std5))
study_specific_est <- rbind(study_specific_est,data.frame(Mean=summary(mod_std5)$coefficients[2,"Estimate"], 
                                                          SD=summary(mod_std5)$coefficients[2,"Std. Error"],
                                                          method="CONFIRM",
                                                          comp="DF vs GA"))


# Use metabin to compute OR and the CI 
library(meta)
mod_std6 <- metabin(rrms_ad$r[2],rrms_ad$n[2],rrms_ad$r[1],rrms_ad$n[1],sm="OR")
study_specific_est <- rbind(study_specific_est,data.frame(Mean=mod_std6$TE, 
                                                          SD=mod_std6$seTE,
                                                          method="Bornstein",
                                                          comp="GA vs PL"))


mod_std7 <- metabin(rrms_ad$r[4],rrms_ad$n[4],rrms_ad$r[3],rrms_ad$n[3],sm="OR")

study_specific_est <- rbind(study_specific_est,data.frame(Mean=mod_std7$TE, 
                                                          SD=mod_std7$seTE,
                                                          method="Johnson",
                                                          comp="GA vs PL"))
rownames(study_specific_est) <- NULL

study_specific_est

