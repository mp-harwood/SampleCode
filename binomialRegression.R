#Manuscript Figure 3A

#####################
# 	Load libraries 	#
#####################
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggeffects)
library(lme4)
library(lmerTest)

#####################
# 	Load data 	#
#####################
ASEwasp<- read.table("/Users/mharwood/Documents/ASE/Aging/ASE_WASPcorrected_ageIRS_December2022", header=TRUE, sep="\t")
eset<- load("/Users/mharwood/Documents/ASE/eset.star.997.17224.Nov2018.rda")

#####################
# 	File setup 	#
#####################
#Fix mismatching IDs
new_id <- paste('C_',sprintf("%04g", as.numeric(sub('C_','',eset.star.997.17224$validated_RNA_ID))), sep='')

#####################
#  Figure 3A 	#
#####################
#Generate table of proportions
sig<- ASEwaso[ASEwasp$BH<0.05,]
sig_ID<- data.frame(table(sig$ID))
colnames(sig_ID)<- c("ID", "sig")
all_ID<- data.frame(table(all$ID))
colnames(all_ID)<- c("ID", "total")
freq<- merge(sig_ID, all_ID, by=c("ID"))
freq$prop_sig<- freq$sig/freq$total
freq$nosig<- freq$total-freq$sig

#merge metadata
age<- ASEwasp %>% dplyr::select(ID, P_AGE, sex, IRSnoAge)
age<- unique(age)
freq_age<- merge(age, freq, by=c("ID"))

#binomial linear regression with age and risk score 
age.model <- glm(cbind(sig,nosig)~P_AGE*risk, freq_age,family="binomial")
age.predict<- ggpredict(age.model,c("P_AGE [all]", "risk"), ci.lvl=0.95)

ggplot(age.predict, aes(x = x, y = predicted, fill=group), label=sprintf("%0.2f", round(a, digits = 2))) +
  stat_smooth(method = "glm", aes(colour=group))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  theme_pubr()+
  xlab("Age")+
  ylab("Proportion of sites with significant ASE")+
  scale_fill_manual(values=c("#D55E00", "#0072B2"), labels=c("High Risk", "Low Risk"))+
  scale_colour_manual(values=c("#D55E00", "#0072B2"), labels=c("High Risk", "Low Risk"))

