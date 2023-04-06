

### Relaimpo...script using root shoot 


#### USING PMVD ########### 

# RELAIMPO _Groemping  ONLY ROOT TRAITS (abajo con leaf, root and soil properties)

# Utilizo las transformadas porq estas fueron las q use en el path analyses. Si no, salen cosas diferentes.  

dataF <- read.csv ("traits_fungal.csv")

dataF$log_rootshoot_ratio= log(dataF$rootshoot_ratio)

#Root_mutual ab (root mass)
Root_Mutualists2.RI <-lm(Mutualists2 ~ log_rootN+ log_rootC+ sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)

metrics_rootMab1 <- calc.relimp(Root_Mutualists2.RI, type = c("pmvd"))
metrics_rootMab <- calc.relimp(Root_Mutualists2.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMab@pmvd = metrics_rootMab@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMab=data.frame(variables=names(metrics_rootMab@pmvd),varexpl.percent=metrics_rootMab@pmvd) # nueva tabla

Tab_metrics_rootMab$variables <- factor(Tab_metrics_rootMab$variables,
                                        levels=c("sqrt_root","log_rootC","sqrt_diameter","log_SRL","log_rootN","sqrt_SRSA","log_RTD"))  # order panels

REL_MUT_rootAB<- ggplot(data=Tab_metrics_rootMab, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Mutualists abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootAB

#Root_mutual ab (root:shoot )
Root_Mutualists2RS.RI <-lm(Mutualists2 ~ log_rootN+ log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)

metrics_rootMab1RS <- calc.relimp(Root_Mutualists2RS.RI, type = c("pmvd"))
metrics_rootMabRS <- calc.relimp(Root_Mutualists2RS.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMabRS@pmvd = metrics_rootMabRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMabRS=data.frame(variables=names(metrics_rootMabRS@pmvd),varexpl.percent=metrics_rootMabRS@pmvd) # nueva tabla

Tab_metrics_rootMabRS$variables <- factor(Tab_metrics_rootMabRS$variables,
                                          levels=c("log_SRL","rootshoot_ratio","sqrt_diameter","log_rootN","log_rootC","sqrt_SRSA", "log_RTD"))  # order panels

REL_MUT_rootABRS<- ggplot(data=Tab_metrics_rootMabRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Mutualists abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootABRS

#########

## Root_S_Mutual  (root mass)

Root_S_Mutual.RI <-lm(log_S_Mutual ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)

metrics_rootMS <- calc.relimp(Root_S_Mutual.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMS@pmvd = metrics_rootMS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMS=data.frame(variables=names(metrics_rootMS@pmvd),varexpl.percent=metrics_rootMS@pmvd) # nueva tabla

Tab_metrics_rootMS$variables <- factor(Tab_metrics_rootMS$variables,
                                       levels=c("log_rootN", "sqrt_diameter","log_RTD","sqrt_SRSA","log_SRL","log_rootC","log_rootshoot_ratio"))  # order panels

REL_MUT_rootMS<-ggplot(data=Tab_metrics_rootMS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Mutualists richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMS


#### # Root_S_Mutual  (INcluding root:shoot  en lugar de root mass)

Root_S_Mutual.RIRS <-lm(log_S_Mutual ~ log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA+log_rootN, data = dataF)

metrics_rootMSRS <- calc.relimp(Root_S_Mutual.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMSRS@pmvd = metrics_rootMSRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMSRS=data.frame(variables=names(metrics_rootMSRS@pmvd),varexpl.percent=metrics_rootMSRS@pmvd) # nueva tabla

Tab_metrics_rootMSRS$variables <- factor(Tab_metrics_rootMSRS$variables,
                                         levels=c("log_RTD", "rootshoot_ratio","log_rootC","log_SRL","sqrt_SRSA", "sqrt_diameter","log_rootN"))  # order panels

REL_MUT_rootMSRS<-ggplot(data=Tab_metrics_rootMSRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Mutualists richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMSRS

#########

### Root_ PCO1_Mutual (root masss)
Root_PCo1_Mutual.RI <-lm(PCo1_Mutual ~ log_rootN + log_rootC+ sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootMC <- calc.relimp(Root_PCo1_Mutual.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMC@pmvd = metrics_rootMC@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMC=data.frame(variables=names(metrics_rootMC@pmvd),varexpl.percent=metrics_rootMC@pmvd) # nueva tabla

Tab_metrics_rootMC$variables <- factor(Tab_metrics_rootMC$variables,
                                       levels=c("log_RTD","sqrt_root","log_rootC","log_rootN","sqrt_diameter","sqrt_SRSA","log_SRL"))  # order panels

REL_MUT_rootMC<-ggplot(data=Tab_metrics_rootMC, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Mutualists composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMC

## Root_ PCO1_Mutual (root:shoot)

Root_PCo1_Mutual.RIRS <-lm(PCo1_Mutual ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootMCRS <- calc.relimp(Root_PCo1_Mutual.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMCRS@pmvd = metrics_rootMCRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMCRS=data.frame(variables=names(metrics_rootMCRS@pmvd),varexpl.percent=metrics_rootMCRS@pmvd) # nueva tabla

Tab_metrics_rootMCRS$variables <- factor(Tab_metrics_rootMCRS$variables,
                                         levels=c("sqrt_diameter","rootshoot_ratio","log_rootC","log_rootN","sqrt_SRSA","log_RTD","log_SRL"))  # order panels

REL_MUT_rootMCRS<-ggplot(data=Tab_metrics_rootMCRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Mutualists composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMCRS



############
#Root_ SAPRO_AB   (root mass)

Root_Saprotrophs.RI <-lm(Saprotrophs2 ~ log_rootN+log_rootC+sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootSab <- calc.relimp(Root_Saprotrophs.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSab@pmvd = metrics_rootSab@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSab=data.frame(variables=names(metrics_rootSab@pmvd),varexpl.percent=metrics_rootSab@pmvd) # nueva tabla

Tab_metrics_rootSab$variables <- factor(Tab_metrics_rootSab$variables,
                                        levels=c("sqrt_root","log_rootN","sqrt_diameter","log_SRL","log_rootC","sqrt_SRSA","log_RTD"))  # order panels

REL_SAP_rootAB<- ggplot(data=Tab_metrics_rootSab, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Saprotrophs abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootAB

#Root_ SAPRO_AB   (root:shoot)

Root_Saprotrophs.RIRS <-lm(Saprotrophs2 ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootSabRS <- calc.relimp(Root_Saprotrophs.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSabRS@pmvd = metrics_rootSabRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSabRS=data.frame(variables=names(metrics_rootSabRS@pmvd),varexpl.percent=metrics_rootSabRS@pmvd) # nueva tabla

Tab_metrics_rootSabRS$variables <- factor(Tab_metrics_rootSabRS$variables,
                                          levels=c("log_rootN","log_SRL","rootshoot_ratio","log_rootC","sqrt_diameter","sqrt_SRSA","log_RTD"))  # order panels

REL_SAP_rootABRS<- ggplot(data=Tab_metrics_rootSabRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Saprotrophs abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootABRS

################
#Root_SAPRO_S  (root mass)
Root_S_Saprotrophs.RI <-lm(S_Sapro ~log_rootN+log_rootC+sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootSS <- calc.relimp(Root_S_Saprotrophs.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSS@pmvd = metrics_rootSS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSS=data.frame(variables=names(metrics_rootSS@pmvd),varexpl.percent=metrics_rootSS@pmvd) # nueva tabla

Tab_metrics_rootSS$variables <- factor(Tab_metrics_rootSS$variables,
                                       levels=c("sqrt_root","log_rootN","log_SRL","log_rootC","sqrt_SRSA","sqrt_diameter","log_RTD"))  # order panels

REL_SAP_rootSS<- ggplot(data=Tab_metrics_rootSS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Saprotrophs richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootSS

##
#Root_SAPRO_S  (root:shoot)
Root_S_Saprotrophs.RIRS <-lm(S_Sapro ~log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootSSRS <- calc.relimp(Root_S_Saprotrophs.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSSRS@pmvd = metrics_rootSSRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSSRS=data.frame(variables=names(metrics_rootSSRS@pmvd),varexpl.percent=metrics_rootSSRS@pmvd) # nueva tabla

Tab_metrics_rootSSRS$variables <- factor(Tab_metrics_rootSSRS$variables,
                                         levels=c("log_rootN","sqrt_diameter","log_rootC","sqrt_SRSA","log_SRL","log_RTD","rootshoot_ratio"))  # order panels

REL_SAP_rootSSRS<- ggplot(data=Tab_metrics_rootSSRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Saprotrophs richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootSSRS

#Root SAPRO_PCOA1 (root mass)
names(dataF)
Root_PCo1_Saprotrophs.RI <-lm(PCo1_Sapro ~ log_rootN+log_rootC+sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootSC <- calc.relimp(Root_PCo1_Saprotrophs.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSC@pmvd = metrics_rootSC@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSC=data.frame(variables=names(metrics_rootSC@pmvd),varexpl.percent=metrics_rootSC@pmvd) # nueva tabla

Tab_metrics_rootSC$variables <- factor(Tab_metrics_rootSC$variables,
                                       levels=c("sqrt_SRSA","sqrt_root","log_rootC","log_rootN","log_RTD","log_SRL","sqrt_diameter")) # order panels

REL_SAP_rootSC<- ggplot(data=Tab_metrics_rootSC, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Saprotrophs composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootSC


##Root SAPRO_PCOA1   ### INcluding Root: shoot 
names(dataF)
Root_PCo1_Saprotrophs.RIRS <-lm(PCo1_Sapro ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootSCRS <- calc.relimp(Root_PCo1_Saprotrophs.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSCRS@pmvd = metrics_rootSCRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSCRS=data.frame(variables=names(metrics_rootSCRS@pmvd),varexpl.percent=metrics_rootSCRS@pmvd) # nueva tabla

Tab_metrics_rootSCRS$variables <- factor(Tab_metrics_rootSCRS$variables,
                                         levels=c("sqrt_diameter","log_RTD","log_rootC","log_rootN","log_SRL","sqrt_SRSA","rootshoot_ratio")) # order panels

REL_SAP_rootSCRS<- ggplot(data=Tab_metrics_rootSCRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Saprotrophs composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootSCRS

#Root PATHO_AB (root mass)

Root_Pathogens.RI <-lm(Pathogens2 ~ log_rootN+log_rootC+sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPab <- calc.relimp(Root_Pathogens.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPab@pmvd = metrics_rootPab@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPab=data.frame(variables=names(metrics_rootPab@pmvd),varexpl.percent=metrics_rootPab@pmvd) # nueva tabla

Tab_metrics_rootPab$variables <- factor(Tab_metrics_rootPab$variables,
                                        levels=c("sqrt_diameter","log_SRL","log_RTD","log_rootC","sqrt_SRSA","sqrt_root","log_rootN"))  # order panels

REL_rootPAT_AB<- ggplot(data=Tab_metrics_rootPab, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Pathogens abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_AB

#Root PATHO_AB (root; shoot )

Root_Pathogens.RIRS <-lm(Pathogens2 ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPabRS <- calc.relimp(Root_Pathogens.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPabRS@pmvd = metrics_rootPabRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPabRS=data.frame(variables=names(metrics_rootPabRS@pmvd),varexpl.percent=metrics_rootPabRS@pmvd) # nueva tabla

Tab_metrics_rootPabRS$variables <- factor(Tab_metrics_rootPabRS$variables,
                                          levels=c("log_rootC","log_RTD","rootshoot_ratio","log_rootN","sqrt_SRSA","log_SRL","sqrt_diameter"))  # order panels

REL_rootPAT_ABRS<- ggplot(data=Tab_metrics_rootPabRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Pathogens abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_ABRS

#RELAIMPO PATHO_S (root mass)

Root_S_Pathogens.RI <-lm(S_Patho ~ log_rootN+log_rootC+sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPS <- calc.relimp(Root_S_Pathogens.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPS@pmvd <- metrics_rootPS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPS <- data.frame(variables=names(metrics_rootPS@pmvd),varexpl.percent=metrics_rootPS@pmvd) # nueva tabla

Tab_metrics_rootPS$variables <- factor(Tab_metrics_rootPS$variables,
                                       levels=c("log_rootN","sqrt_root","log_SRL","sqrt_SRSA","log_rootC","sqrt_diameter","log_RTD"))  # order panels

REL_rootPAT_PS<- ggplot(data=Tab_metrics_rootPS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Pathogens richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_PS

#RELAIMPO PATHO_S (root;shoot)

Root_S_Pathogens.RIRS <-lm(S_Patho ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPSRS <- calc.relimp(Root_S_Pathogens.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPSRS@pmvd <- metrics_rootPSRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPSRS <- data.frame(variables=names(metrics_rootPSRS@pmvd),varexpl.percent=metrics_rootPSRS@pmvd) # nueva tabla

Tab_metrics_rootPSRS$variables <- factor(Tab_metrics_rootPSRS$variables,
                                         levels=c("log_rootN","log_rootC","sqrt_SRSA","log_RTD","sqrt_diameter","log_SRL","rootshoot_ratio"))  # order panels

REL_rootPAT_PSRS<- ggplot(data=Tab_metrics_rootPSRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Pathogens richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_PSRS

#RELAIMPO PATHO_PCOA2  (root mass)
Root_PCo2_Pathogens.RI <-lm(PCo2_Patho ~ log_rootN+log_rootC+sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPC <- calc.relimp(Root_PCo2_Pathogens.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPC@pmvd = metrics_rootPC@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPC=data.frame(variables=names(metrics_rootPC@pmvd),varexpl.percent=metrics_rootPC@pmvd) # nueva tabla

Tab_metrics_rootPC$variables <- factor(Tab_metrics_rootPC$variables,
                                       levels=c("sqrt_diameter","log_SRL","sqrt_root","log_RTD","log_rootC","log_rootN","sqrt_SRSA"))  # order panels

REL_rootPAT_PC<- ggplot(data=Tab_metrics_rootPC, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Pathogens composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_PC

#RELAIMPO PATHO_PCOA2 #### Including Root:shoot 
names(dataF)
Root_PCo2_Pathogens.RIRS <-lm(PCo2_Patho ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPCRS <- calc.relimp(Root_PCo2_Pathogens.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPCRS@pmvd = metrics_rootPCRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPCRS=data.frame(variables=names(metrics_rootPCRS@pmvd),varexpl.percent=metrics_rootPCRS@pmvd) # nueva tabla

Tab_metrics_rootPCRS$variables <- factor(Tab_metrics_rootPCRS$variables,
                                         levels=c("log_SRL","log_rootN","sqrt_diameter","log_rootC","log_RTD","sqrt_SRSA","rootshoot_ratio"))  # order panels

REL_rootPAT_PCRS<- ggplot(data=Tab_metrics_rootPCRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Pathogens composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_PCRS

#Using root mass rather than root_shoot
ggarrange (REL_rootPAT_AB, REL_rootPAT_PS, REL_rootPAT_PC,REL_SAP_rootAB,
           REL_SAP_rootSS,REL_SAP_rootSC, REL_MUT_rootAB, REL_MUT_rootMS, REL_MUT_rootMC,
           labels = c("A", "B","C","D","E","F","G", "H", "I"),
           common.legend = FALSE, align = c("hv"))

## USing Root:Shoot rather than root mass
ggarrange (REL_rootPAT_ABRS, REL_rootPAT_PSRS, REL_rootPAT_PCRS,REL_SAP_rootABRS,
           REL_SAP_rootSSRS,REL_SAP_rootSCRS, REL_MUT_rootABRS, REL_MUT_rootMSRS, REL_MUT_rootMCRS,
           labels = c("A", "B","C","D","E","F","G", "H", "I"),
           common.legend = FALSE, align = c("hv"))
###################################################