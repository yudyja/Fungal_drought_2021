install.packages("car")
install.packages("carData")
install.packages("corrplot")
install.packages ("vegan")
install.packages("biclust")
install.packages("devtools")
install.packages("biomod2")
install.packages("ade4")
install.packages('VennDiagram')
install.packages("ggpubr")
install.packages(c("FactoMineR", "factoextra"))
install.packages("MVN")
install.packages ("RVAideMemoire")
install.packages ("mvnormtest")
install.packages("Rmisc")
install.packages ("doBy")
install.packages("Rtools")
install.packages("relaimpo")
install.packages("C:/Users/aurel/Desktop/CURSO IGAC/relaimpo_2.2-3.zip", repos = NULL, type = "source")# type = "win.binary")
# installed Relaimpo from Linux console   https://anaconda.org/conda-forge/r-relaimpo

install.packages("C:/Users/aurel/Desktop/CURSO IGAC/relaimpo_2.2-3.zip", repos=NULL, type="win.binary")

library(ggplot2)
library(car)
require(corrplot)
library (vegan)
library(ade4)
library(lavaan)
library(ggpubr)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(nlme)
library(lme4)
library(MASS)
library(biclust)
library(MVN)
library(mvnormtest)
library(Rmisc)
library(BiodiversityR)
library(relaimpo)
library(tidyr)
library(multcomp)
library(sandwich)
library(relaimpo)

setwd("/home/yudyja/Escritorio/FUNGAL_drought/Fungaltraits")
getwd()
#data = read.csv("mydatafungalt.csv")  # DATA OLD. Sol ousando Funguilg but no  FUNGUILD y DATAB Carlos
dataF <- read.csv ("traits_fungal.csv")


dataF$log_Pathogen <- log(dataF$Pathogens2+5)
dataF$log_Sapro <- log(dataF$Saprotrophs2+5)
dataF$log_Mutual <- log(dataF$Mutualists2+5)

dataF$log_shootmass <- log(dataF$shootmass)

dataF$log_S_Patho <- log(dataF$S_Patho+5)
dataF$log_S_Sapro <- log(dataF$S_Sapro+5)
dataF$log_S_Mutual <- log(dataF$S_Mutual+5)

dataF$log_PCo1_Patho <- log (dataF$PCo1_Patho+5)  # los PCoAs salieron de Infostat
dataF$log_PCo1_Sapro <- log(dataF$PCo1_Sapro+5)
dataF$log_PCo1_Mutual <- log(dataF$PCo1_Mutual+5)

dataF$log_PCo2_Patho <- log (dataF$PCo2_Patho+5) 
dataF$log_PCo2_Sapro <- log(dataF$PCo2_Sapro+5)
dataF$log_PCo2_Mutual <- log(dataF$PCo2_Mutual+5)

dataF$log_rootN <- log(dataF$rootN)
dataF$log_rootC <- log(dataF$rootC)
dataF$sqrt_root <- sqrt(dataF$rootmass)
dataF$sqrt_diameter <- sqrt(dataF$diameter)
dataF$log_RTD <- log(dataF$RTD)
dataF$log_SRL <- log(dataF$SRL)
dataF$sqrt_SRSA <- sqrt(dataF$SRSA)

dataF$rootshoot_ratio <- dataF$rootmass/dataF$shootmass

dataF<-unite(dataF, sptto,c("Specie", "Tto"),  sep = "_", remove = FALSE) # create a new column uniting these two. 
dataF$sptto<-as.factor(dataF$sptto)

##### Root shoot ratio

#Este orden de species salio de total fungal y lo use para todas las graficas 
dataF$Specie <- factor(dataF$Specie, levels=c("H_Rumex","H_Ranunculus","G_Lolium","L_Trifolium", "G_Dactylis","H_Armeria","H_Silene","G_F.rubra",
                                              "G_Poa","H_Galium","L_Medicago",  "G_F.brevipila","H_Berteroa","H_Achillea","H_Artemisia","G_Holcus","G_Anthoxanthum", 
                                              "H_Hypericum","L_Vicia","H_Potentilla","H_Daucus","H_Hieracium","G_Arrhenatherum", "H_Plantago"))  # order panels


rootshootratio <- summarySE(dataF, measurevar="rootshoot_ratio", groupvars=c("Specie","Tto"))
rootshootratio2 <- summarySE(dataF, measurevar="rootshoot_ratio", groupvars=c("Tto"))

rootshootratioplotOK<-ggplot(rootshootratio, aes(x=Specie, y=rootshoot_ratio, colour=Tto, group =Tto))+# ylim (-0.5,5)+
  scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
                             "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
                            "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=rootshoot_ratio-se, ymax=rootshoot_ratio+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+ #jitter using summarySE
  scale_color_manual(values=c("lightblue","darkgoldenrod3"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 10))+
  labs(x='', y="Root:Shoot") #+
rootshootratioplotOK

ggarrange (rootshootratioplotOK, rootshootratioplotOK,
          labels = c("A", "B"), legend = c("top"), common.legend = TRUE, ncol = 1, nrow = 2)

#NORMALITY

# the most normal data (infostat R value are: 
#abundancia :log_pathogen5, sqrt_sapro, sqrt_AMF (0.79),  
#richness: spathogen, sqrt_S_sapro, sqrt_S_AMF (0.84), 
# #shoot mass. it is better log shoot mass


##################    PER-MANOVA  ##################

#SPECIES AS FACTOR

### PERMANOVA Total
NP_total <- read.csv ("ESV_totalFungComposition_mantel.csv")
esv_total <-NP_total %>% dplyr::select(starts_with("esv"))
View(esv_total)

#esv_total.adj <-subset(esv_total, select = -c(esv1462,esv2447,esv2818,esv3183,esv2531,esv2604,esv1811,esv2013,
#                                              esv2603,esv2954,esv3286,esv3395,esv3396,esv3440)) # quitar unos ESV que son EMF (-c quitar column)
# NO las quite porq tb son FUNGI 
adonis (esv_total ~ Species*Tto, data = NP_total, permutations = 999,  method = "bray")

# Permanova
fungal_var <- cbind(dataF$Pathogens2, dataF$Saprotrophs2, dataF$log_mutual, dataF$S_Patho, dataF$S_Sapro, dataF$S_Mutual,
                   dataF$PCo1_Patho, dataF$PCo1_Sapro, dataF$PCo1_Mutual)
adonis (fungal_var ~ Specie*Tto, data = dataF, permutations = 99)

###
fungal_var2 <- cbind(dataF$log_Pathogen, dataF$log_Sapro, dataF$log_Mutual, dataF$log_S_Patho, dataF$log_S_Sapro, dataF$log_S_Mutual,
                    dataF$log_PCo1_Patho, dataF$log_PCo1_Sapro, dataF$log_PCo1_Mutual)
adonis (fungal_var2 ~ Specie*Tto, data = dataF, permutations = 999) ## ESCOGIDO! 

### PERMANOVA Pathogens
NP_patho<- read.csv("PAST_Pathogen_composition_NEW.csv")

NP_patho2<-NP_patho[-which(NP_patho$ID== "34" | NP_patho$ID== "61" | NP_patho$ID== "97" ),] # quitar estas plantas q tiene solo zeros 
esv_patho <-NP_patho2 %>% dplyr::select(starts_with("esv"))
adonis (esv_patho ~ species*tto, data = NP_patho2, permutations = 999,  method = "bray")

### PERMANOVA Saprotrophs
NP_sapro<- read.csv("PAST_Saprotro_Composition_NEW.csv")
esv_sapro <-NP_sapro %>% dplyr::select(starts_with("esv"))
adonis (esv_sapro ~ Species*Tto, data = NP_sapro, permutations = 999,  method = "bray")

### PERMANOVA Mutualists
NP_mutual<- read.csv("PAST_Mutualists_Composition_NEW.csv")
NP_mutual2<-NP_mutual[-which(NP_mutual$ID== "34" | NP_mutual$ID== "141"),]  # those plants that only have zeros 
esv_mutual <-NP_mutual2 %>% dplyr::select(starts_with("esv"))

esv_mutual.adj <-subset(esv_mutual, select = -c(esv1462,esv2447,esv2818,esv3183,esv2531,esv2604,esv1811,esv2013,
                                              esv2603,esv2954,esv3286,esv3395,esv3396,esv3440)) # quitar unos ESV que son EMF (-c quitar column)
adonis (esv_mutual.adj ~ Specie*Tto, data = NP_mutual2, permutations = 999,  method = "bray")


#PLANT GROUP AS FACTOR

#abundance, richness and composition.. ###Results similar to the previous one
manPG <- manova(cbind(Pathog_comp_1, Sapro_Comp_1, AMF_comp_1, log_pathogen5, sqrt_sapro, sqrt_AMF, S_pathog, sqrt_Ssapro, sqrt_Samf) ~ PlantGroup + Tto + species:PlantGroup, data = data)
summary(manPG)  

#########################################
#####################
###################  NMDS y las barras (el centroide de drought menos el de control por plant specie)

dataESV <- read.csv("ESV_totalFungComposition_mantel.csv")  # la tabla donde esta la "total" fungal composition. The ESV selected
totESV<-dataESV[, grep(pattern="^esv", colnames(dataESV))]  # seleccionar las columnas que empiezan con "esv" 
dis <- vegdist(totESV, method = "bray") 

groups <- factor(c(rep(1,5),rep(2,5),rep(3,5), rep(4,5),rep(5,5),rep(6,5),rep(7,4),rep(8,5),rep(9,5),rep(10,5),
                    rep(11,5),rep(12,5),rep(13,5),rep(14,5),rep(15,5),rep(16,5),rep(17,5),rep(18,5),rep(19,5),rep(20,5),
                    rep(21,5),rep(22,5),rep(23,5),rep(24,5),rep(25,5),rep(26,5),rep(27,5),rep(28,5),rep(29,5),rep(30,5),
                    rep(31,5),rep(32,5),rep(33,5),rep(34,4),rep(35,5),rep(36,5),rep(37,5),rep(38,5),rep(39,5),rep(40,5),
                    rep(41,5),rep(42,5),rep(43,4),rep(44,5),rep(45,5),rep(46,4),rep(47,5),rep(48,5)), #(#grupo, how many rows)
                  labels = c("G_AnthoxantControl","G_AnthoxantDrought","BerteroaControl","BerteroaDrought","DaucusControl","DaucusDrought",
                             "FbreviControl","FbreviDrought","AchilleaControl","AchilleaDrought","PotentillaControl","PotentillaDrought",
                             "PlantagoControl","PlantagoDrought","HieraciumControl","HieraciumDrought","ArtemisiaControl","ArtemisiaDrought",
                             "HolcusControl","HolcusDrought","ArrhenatherumControl","ArrhenatherumDrought","ViciaControl","ViciaDrought",
                             "HypericumControl","HypericumDrought","FrubraControl","FrubraDrought","GaliumControl","GaliumDrought","PoaControl",
                             "PoaDrought","SileneControl","SileneDrought","TrifoliumControl","TrifoliumDrought","DactylisControl","DactylisDrought",
                             "RumexControl","RumexDrought","MedicagoControl","MedicagoDrought","LoliumControl","LoliumDrought","RanunculusControl",
                             "RanunculusDrought","ArmeriaControl","ArmeriaDrought"))

betadisper(dis,groups) # Rsto control -drought y se si la diferencia entre las fungal comun entre drought y control por cada plant sp.

#options(max.print=10000) # lo pongo para q muestre toooodas las filas en la consola. Pero al convertirlo a data.frame y exportarlo salen todas anyway
TukeyHSD(betadisper(dis,groups)) # para ver todas las posibles combinaciones
as.data.frame(Difference$group)
write.csv(Difference$group,"Difference_centroids_NMDS.csv")  #group es como se llama en el result los datos que necesito extraer 


#### PCOA (Figure 2 paper). Usamos este en lugar del NMDS ( porq siempre habiamos usado PCoA)

spec.pco <- wcmdscale(vegdist(OTU_abundance6[,c(1:2113)],"bray"), eig=TRUE) ## NO lo he corrido .. revisarlo..algun dia!!

# en el paper hay dos PCoA( taxonomic) que incluye todas las ESV, y (functional) que incluye aquellas asignadas como 
#patho, sapro, mutual, mix. Aqui no estan incluidas 38 ESVs q solo fueron clasificadas como phyla. Ademas, en el "functional" 
#solo se tienen 4 columnas de datos ( se sumaron todas las abundancias brutas de los pathog, o sea el numero de reads after bioinform)
#####################################################
##ABUNDANCE

#Este orden es con base en Total fungal diveristy (descendent).Usare este orden para todos los graphs de AB y S
dataF$Specie <- factor(dataF$Specie, levels=c("H_Rumex","H_Ranunculus","G_Lolium","L_Trifolium", "G_Dactylis","H_Armeria","H_Silene","G_F.rubra",
                 "G_Poa","H_Galium","L_Medicago",  "G_F.brevipila","H_Berteroa","H_Achillea","H_Artemisia","G_Holcus","G_Anthoxanthum", 
                 "H_Hypericum","L_Vicia","H_Potentilla","H_Daucus","H_Hieracium","G_Arrhenatherum", "H_Plantago"))  # order panels

#ABUNDANCE LM
str(dataF)
#Pathogens
dataF$log_Pathogens2 <- log(dataF$Pathogens2+5)
pathoA <-lm(log_Pathogens2~Specie*Tto, data=dataF) 
qqnorm(residuals(pathoA))
plot((pathoA), add.smooth = FALSE, which = 1)
res.pathoA<-residuals(pathoA)
bartlett.test(res.pathoA, dataF$Specie) #no
bartlett.test(res.pathoA, dataF$Tto) #ok

vf1 <- varIdent(form= ~1|Specie)
Mgls.pathoA <- gls(log_Pathogens2~Specie*Tto, data =dataF,weights = vf1)
plot(Mgls.pathoA, add.smooth=FALSE, which=1)              
E.pathoA<-resid(Mgls.pathoA, type="normalized") 
bartlett.test(E.pathoA, dataF$Specie)  #ok  
bartlett.test(E.pathoA, dataF$Tto)    #Ok      
shapiro.test(E.pathoA) #ok
qqnorm(E.pathoA)

anova(Mgls.pathoA) # new model
#ggplot pathoAB..... run Lines before to get the right order
PathogensAB <- summarySE(dataF, measurevar="Pathogens2", groupvars=c("Specie","Tto"))
PathogensAB2 <- summarySE(dataF, measurevar="Pathogens2", groupvars=c("Specie"))

PathogenABplotOK<-ggplot(PathogensAB, aes(x=Specie, y=Pathogens2, colour=Tto, group =Tto))+ ylim (0,15)+
  #scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
  #                            "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
  #                            "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=Pathogens2-se, ymax=Pathogens2+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+ #jitter using summarySE
  scale_color_manual(values=c("blue","red"))+
 theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5,size = 10))+
  labs(x='', y="Pathogen abundance (%)") #+
PathogenABplotOK

#Saprotros
dataF$log_Saprotrophs2 <- log(dataF$Saprotrophs2+5)
saproA <-lm(log_Saprotrophs2~Specie*Tto, data=dataF) 
qqnorm(residuals(saproA))
plot((saproA), add.smooth = FALSE, which = 1)
res.saproA<-residuals(saproA)
bartlett.test(res.saproA, dataF$Specie) #no
bartlett.test(res.saproA, dataF$Tto) #yes
ggplot(dataF, aes(x=Specie, y=Saprotrophs, fill=Tto)) + geom_boxplot()

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.saproA <- gls(log_Saprotrophs2~Specie*Tto, data =dataF,weights = vf1)
plot(Mgls.saproA, add.smooth=FALSE, which=1)              
E.saproA<-resid(Mgls.saproA, type="normalized") 
bartlett.test(E.saproA, dataF$Specie)  #ok  
bartlett.test(E.saproA, dataF$Tto)    #no      
shapiro.test(E.saproA) #ok
qqnorm(E.saproA)

anova(Mgls.saproA) # new model
#ggplot SaproAB  .....run Lines before to get the right order
SaprotrophsAB <- summarySE(dataF, measurevar="Saprotrophs2", groupvars=c("Specie","Tto"))
saproABplotOK<-ggplot(SaprotrophsAB, aes(x=Specie, y=Saprotrophs2, colour=Tto, group =Tto)) +ylim (0,45)+
  #scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
  #                            "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
  #                            "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=Saprotrophs2-se, ymax=Saprotrophs2+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5,size = 10))+
  labs(x='', y="Saprotroph abundance (%)")#+
saproABplotOK

#Mutualists (log)
dataF$log_mutual2 <- log(dataF$Mutualists2+5)  
mutualA <-lm(log_mutual2~Specie*Tto, data=dataF) 
qqnorm(residuals(mutualA))
plot((mutualA), add.smooth = FALSE, which = 1)
res.mutualA<-residuals(mutualA)
bartlett.test(res.mutualA, dataF$Specie) #no
bartlett.test(res.mutualA, dataF$Tto) #yes
ggplot(dataF, aes(x=Specie, y=Mutualists, fill=Tto)) + geom_boxplot()

vf1 <- varIdent(form= ~1|Specie)
Mgls.mutualA <- gls(log_mutual2~Specie*Tto, data =dataF,weights = vf1)  # mejor con LOG que con sqrt q con nada 
plot(Mgls.mutualA, add.smooth=FALSE, which=1)              
E.mutualA<-resid(Mgls.mutualA, type="normalized") 
bartlett.test(E.mutualA, dataF$Specie)  #ok  
bartlett.test(E.mutualA, dataF$Tto)    #Ok      
shapiro.test(E.mutualA) #ok
qqnorm(E.mutualA)

anova(Mgls.mutualA) # new model
#ggplot MutualAB ## run Lines before to get the right order
MutualistsAB <- summarySE(dataF, measurevar="Mutualists2", groupvars=c("Specie","Tto"))
mutualABplotOK<-ggplot(MutualistsAB, aes(x=Specie, y=Mutualists2, colour=Tto, group =Tto)) + ylim (0,15)+
 # scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
 #                           "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
  #                          "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=Mutualists2-se, ymax=Mutualists2+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, size = 10))+
  labs(x='', y="Mutualist abundnace (%)")#+
mutualABplotOK

ggarrange (PathogenABplotOK,saproABplotOK, mutualABplotOK ,
           labels = c("A", "B", "C"), legend = c("none"), common.legend = TRUE, ncol = 1, nrow = 3)

##########################################
##########ABUNDANCE  analysis  GLM (species and TTo as FACTOR)

# variables have to be "as.factor"  DA MAS SIGNIFCNACIAS CON EL LINEAL MODEL ATRAS 

#pathogen
#dataF$Pathogens2_01 <- dataF$Pathogens2/100
patho1 <- glm(Pathogens2_01 ~ Specie*Tto, family=quasibinomial(link='logit'),data=dataF)
aovpatho <-anova(patho1, test = "F")
aovpatho
aovpatho$Deviance/aovpatho$Df
summary (aovpatho)
#saprotro
#dataF$Saprotrophs2 <- dataF$Saprotrophs2/100
sapro1 <- glm(Saprotrophs2 ~ Specie*Tto, family=quasibinomial(link='logit'),data=dataF)
aovsapro <-anova(sapro1, test = "F")
aovsapro
aovsapro$Deviance/aovsapro$Df
summary (aovsapro)
#Mutualist
#dataF$Mutualists2 <- dataF$Mutualists2/100
mutual1 <- glm(Mutualists2 ~ Specie*Tto, family=quasibinomial(link='logit'),data=dataF)
aovmutual <-anova(mutual1, test = "F")
aovmutual
aovmutual$Deviance/aovmutual$Df
summary (aovmutual)

####################################################
#PLANT SPECIES AS RANDOM  (Abundance and richness)

names(dataF)
dataF$fspecies = factor(dataF$Specie)

d_2<-dataF$Tto:dataF$PlantGroup:dataF$Specie

### ABUNDANCIA
#pathogens (sp as random)
summary(anova(lmer(Pathogens2 ~ Tto+(1|d_2),data = dataF))) #no p values
summary(aov(Pathogens2 ~ Tto+Error(d_2),data = dataF))
Mlme1 <- lme(Pathogens2 ~ Tto, random = ~1 | fspecies, data = dataF)
summary(Mlme1)

plot(Mlme1, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.p<-resid(Mlme1, type="normalized") 
qqnorm(E.p)
anova(Mlme1)  # no signif
 
#sapro (sp as random)
summary(aov(Saprotrophs2 ~ Tto+Error(d_2),data = dataF))
Mlme2 <- lme(Saprotrophs2 ~ Tto, random = ~1 | fspecies, data = dataF)
summary(Mlme2)
plot(Mlme2, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.s<-resid(Mlme2, type="normalized") 
qqnorm(E.s)
anova(Mlme2) # dio lo mismo usando sqrt_sapro y log_sapro5  (no signif)

#Mutualists (sp as random)
summary(aov(log_mutual2 ~ Tto+Error(d_2),data = dataF))
Mlme3 <- lme(log_mutual2 ~ Tto, random = ~1 | fspecies, data = dataF)
summary(Mlme3)
plot(Mlme3, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.a<-resid(Mlme3, type="normalized") 
qqnorm(E.a)
anova(Mlme3) # con sqrt y log y nada no da signif ( mas cercano sin nada p=0.22)

#pathogen richness (sp as random)
names(data)

Mlme4 <- lme(S_Patho ~ Tto, random = ~1 | fspecies, data = dataF)
summary(Mlme4)
plot(Mlme4, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.Sp<-resid(Mlme4, type="normalized") 
qqnorm(E.Sp)
anova(Mlme4)

#sapro richness (sp as random)

Mlme5 <- lme(S_Sapro ~ Tto, random = ~1 | fspecies, data = dataF)
summary(Mlme5)
plot(Mlme5, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.Ss<-resid(Mlme5, type="normalized") 
qqnorm(E.Ss)
M.gls <- gls(S_Sapro ~ Tto, method = "REML",
               correlation = corCompSymm(form =~ 1 | fspecies), data = dataF)
anova(Mlme5)
anova(M.gls)

#Mutualists richness. (sp as random) Este gls y lme (como estan escritos) son lo mismo
M.gls2 <- gls(S_Mutual ~ Tto, method = "REML",
             correlation = corCompSymm(form =~ 1 | fspecies), data = dataF)
anova(M.gls2)
plot(M.gls2, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.Sa<-resid(M.gls2, type="normalized") 
qqnorm(E.Sa)

Mlme6 <- lme( S_Mutual~ Tto, random = ~1 | fspecies, data = dataF)
anova(Mlme6)

#Total richness. (sp as random)  # All were assigned as FUNGI (even if they were not with some funtional group)
M.gls8 <- gls(S_ALL ~ Tto, method = "REML",
              correlation = corCompSymm(form =~ 1 | fspecies), data = dataF)
anova(M.gls8)
plot(M.gls8, add.smooth=FALSE, which=1)    # resplot(E.p, add.smooth=FALSE, which=1)
E.ALL<-resid(M.gls8, type="normalized") 
qqnorm(E.ALL)

Mlme8 <- lme( S_ALL~ Tto, random = ~1 | fspecies, data = dataF)
anova(Mlme8)

#

##################################################################

# RICHNESS USANDO GLM  (solo "species" es signficativo)

dataF$S_Patho01 <- dataF$S_Patho/100
R_patho <- glm(S_Patho01 ~ Specie*Tto, family=quasipoisson(link='logit'),data=dataF)
R_patho_aov <-anova(R_patho, test = "F")
R_patho_aov

dataF$S_sapro01 <- dataF$S_Sapro/100
R_sapro <- glm(S_sapro01 ~ Specie*Tto, family=quasipoisson(link='logit'),data=dataF)
R_sapro_aov <-anova(R_sapro, test = "F")
R_sapro_aov

dataF$S_mutual01 <- dataF$S_Mutual/100
R_mutual <- glm(S_mutual01 ~ Specie*Tto, family=quasipoisson(link='logit'),data=dataF)
R_mutual_aov <-anova(R_mutual, test = "F")
R_mutual_aov

#################### RICHNESS USANDO LM ### CHOOOOOSEEENNNNNNNN 

#Multiple comparisons
### Run the following lines. These introduce methods for 'gls' objects
#So that we can use multiple comparisons. whitout this the mc no funciona!!!! 

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}

# S_patho
R_patholm <-lm(S_Patho~Specie*Tto, data=dataF) 
qqnorm(residuals(R_patholm))
plot((R_patholm), add.smooth = FALSE, which = 1)
res.R_patholm<-residuals(R_patholm)
bartlett.test(res.R_patholm, dataF$Specie) #no
bartlett.test(res.R_patholm, dataF$Tto) #yes

vf1 <- varIdent(form= ~1|Specie)
Mgls.R_patholm <- gls(S_Patho~Specie*Tto, data =dataF,weights = vf1)  
plot(Mgls.R_patholm, add.smooth=FALSE, which=1)              
E.R_patholm<-resid(Mgls.R_patholm, type="normalized") 
bartlett.test(E.R_patholm, dataF$Specie)  #ok  
bartlett.test(E.R_patholm, dataF$Tto)    #Ok      
shapiro.test(E.R_patholm) #ok
qqnorm(E.R_patholm)

anova(Mgls.R_patholm) 
#Ggplot Patho richness   ( RUN THE ORDER OF PLANTS ABOVE.... )
PathoRIC <- summarySE(dataF, measurevar="S_Patho", groupvars=c("Specie","Tto"))
PathoRICplotOK<-ggplot(PathoRIC, aes(x=Specie, y=S_Patho, colour=Tto, group =Tto)) +ylim (0,21)+
  #scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
  #                           "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
   #                         "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=S_Patho-se, ymax=S_Patho+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 0, hjust=1,size = 10))+
  labs(x='', y="Pathogen richness (#)")#+
PathoRICplotOK

# S_sapro

R_saprolm <-lm(S_Sapro~Specie*Tto, data=dataF) 
qqnorm(residuals(R_saprolm))
plot((R_saprolm), add.smooth = FALSE, which = 1)
res.R_saprolm<-residuals(R_saprolm)
bartlett.test(res.R_saprolm, dataF$Specie) #no
bartlett.test(res.R_saprolm, dataF$Tto) #yes

vf1 <- varIdent(form= ~1|Specie)
Mgls.R_saprolm <- gls(S_Sapro~Specie*Tto, data =dataF,weights = vf1)  
plot(Mgls.R_saprolm, add.smooth=FALSE, which=1)              
E.R_saprolm<-resid(Mgls.R_saprolm, type="normalized") 
bartlett.test(E.R_saprolm, dataF$Specie)  #ok  
bartlett.test(E.R_saprolm, dataF$Tto)    #Ok      
shapiro.test(E.R_saprolm) #ok
qqnorm(E.R_saprolm)

anova(Mgls.R_saprolm) 

#ggplot _Sapro_Richness

str(dataF)
SaproRIC <- summarySE(dataF, measurevar="S_Sapro", groupvars=c("Specie","Tto"))
SaproRICplotOK<-ggplot(SaproRIC, aes(x=Specie, y=S_Sapro, colour=Tto, group =Tto)) +ylim (0,76)+
  scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
                              "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
                              "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=S_Sapro-se, ymax=S_Sapro+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+
  scale_color_manual(values=c("lightblue","darkgoldenrod3"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 10))+
  labs(x='', y="Saprotroph richness (#)")#+
SaproRICplotOK


# S_mutual

#dataF$sqrt_S_Mutual <- sqrt(dataF$S_Mutual)  # con log y sqrt no mejora 
R_mutualm <-lm(S_Mutual~Specie*Tto, data=dataF) 
qqnorm(residuals(R_mutualm))
plot((R_mutualm), add.smooth = FALSE, which = 1)
res.R_mutualm<-residuals(R_mutualm)
bartlett.test(res.R_mutualm, dataF$Specie) #no
bartlett.test(res.R_mutualm, dataF$Tto) #no

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.R_mutualm <- gls(S_Mutual~Specie*Tto, data =dataF,weights = vf1,vf2)  
plot(Mgls.R_mutualm, add.smooth=FALSE, which=1)              
E.R_mutualm<-resid(Mgls.R_mutualm, type="normalized") 
bartlett.test(E.R_mutualm, dataF$Specie)  #ok  
bartlett.test(E.R_mutualm, dataF$Tto)    #Ok      
shapiro.test(E.R_mutualm) #ok
qqnorm(E.R_mutualm)

anova(Mgls.R_mutualm) 

#ggplot _Mutual_Richness
MutualRIC <- summarySE(dataF, measurevar="S_Mutual", groupvars=c("Specie","Tto"))
MutualRICplotOK<-ggplot(MutualRIC, aes(x=Specie, y=S_Mutual, colour=Tto, group =Tto)) + ylim (0,30)+
  scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
                              "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
                              "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=S_Mutual-se, ymax=S_Mutual+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 1, data = dataF)+
  scale_color_manual(values=c("lightblue","darkgoldenrod3"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 10))+
  labs(x='', y="Mutualist richness (#)")#+
MutualRICplotOK

# S_ALL 
R_ALLlm <-lm(S_ALL~Specie*Tto, data=dataF) 
qqnorm(residuals(R_ALLlm))
plot((R_ALLlm), add.smooth = FALSE, which = 1)
res.R_ALLlm<-residuals(R_ALLlm)
bartlett.test(res.R_ALLlm, dataF$Specie) #no
bartlett.test(res.R_ALLlm, dataF$Tto) #Yes

vf1 <- varIdent(form= ~1|Specie)
Mgls.R_ALLlm <- gls(S_ALL~Specie*Tto, data =dataF,weights = vf1)  
plot(Mgls.R_ALLlm, add.smooth=FALSE, which=1)              
E.R_ALLlm<-resid(Mgls.R_ALLlm, type="normalized") 
bartlett.test(E.R_ALLlm, dataF$Specie)  #ok  
bartlett.test(E.R_ALLlm, dataF$Tto)    #Ok      
shapiro.test(E.R_ALLlm) #ok
qqnorm(E.R_ALLlm)

anova(Mgls.R_ALLlm) 

#ggplot richness
dataF$Specie <- factor(dataF$Specie, levels=c("H_Rumex","H_Ranunculus","G_Lolium","L_Trifolium",
                "G_Dactylis","H_Armeria","H_Silene","G_F.rubra","G_Poa","H_Galium","L_Medicago",
                "G_F.brevipila","H_Berteroa","H_Achillea","H_Artemisia","G_Holcus","G_Anthoxanthum",
                "H_Hypericum","L_Vicia","H_Potentilla","H_Daucus","H_Hieracium","G_Arrhenatherum", "H_Plantago"))  # order panels

ALLRIC <- summarySE(dataF, measurevar="S_ALL", groupvars=c("Specie","Tto"))
ALLRICPP <- summarySE(dataF, measurevar="S_ALL", groupvars=c("Specie"))

ALLRICplotOK<-ggplot(ALLRIC, aes(x=Specie, y=S_ALL, colour=Tto, group =Tto)) +ylim (0,400)+
 # scale_x_discrete(labels = c("Rumex","Ranunculus","Lolium","Trifolium", "Dactylis","Armeria","Silene",
#    "F.rubra","Poa","Galium","Medicago","F.brevipila","Berteroa","Achillea","Artemisia","Holcus",
#    "Anthoxanthum","Hypericum","Vicia","Potentilla","Daucus","Hieracium","Arrhenatherum","Plantago"))+
  geom_errorbar(aes(ymin=S_ALL-se, ymax=S_ALL+se), width=0,position = position_dodge(0.8))+
  geom_point(aes(color = Tto), size=3, position = position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex = 0.6, data = dataF)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=8, angle=0))+
  theme(axis.text.x = element_text(angle = 0, hjust=1,size = 10))+
  labs(x='', y="Total fungal richness (#)")+
  theme(legend.position='top')
ALLRICplotOK

ggarrange (ALLRICplotOK,ALLRICplotOK,ALLRICplotOK,
           labels = c("", "B"), legend = c("top"), common.legend = TRUE, ncol = 1, nrow = 3)

ggarrange (PathoRICplotOK,PathogenABplotOK,
           labels = c("A", "B"), legend = c("top"), common.legend = TRUE, ncol = 1, nrow = 2)
  
ggarrange (SaproRICplotOK,saproABplotOK,
           labels = c("A", "B"), legend = c("top"), common.legend = TRUE, ncol = 1, nrow = 2)

ggarrange (MutualRICplotOK,mutualABplotOK,
           labels = c("A", "B"), legend = c("top"), common.legend = TRUE, ncol = 1, nrow = 2)



################### COMPOSITION using PCoA 1

#Patho composit
PCo1_Patho <-lm(PCo1_Patho~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo1_Patho))
plot((PCo1_Patho), add.smooth = FALSE, which = 1)
res.PCo1_Patho<-residuals(PCo1_Patho)
bartlett.test(res.PCo1_Patho, dataF$Specie) #no
bartlett.test(res.PCo1_Patho, dataF$Tto) #Yes

vf1 <- varIdent(form= ~1|Specie)
Mgls.PCo1_Patho <- gls(PCo1_Patho~Specie*Tto, data =dataF,weights = vf1)  
plot(Mgls.PCo1_Patho, add.smooth=FALSE, which=1)              
E.PCo1_Patho<-resid(Mgls.PCo1_Patho, type="normalized") 
bartlett.test(E.PCo1_Patho, dataF$Specie)  #ok  
bartlett.test(E.PCo1_Patho, dataF$Tto)    #Ok      
shapiro.test(E.PCo1_Patho) #ok
qqnorm(E.PCo1_Patho)

anova(Mgls.PCo1_Patho) 
ggplot(dataF, aes(x=Specie, y=S_Total, fill=Tto)) + geom_boxplot()

#PCo2 patho
PCo2_Patho <-lm(PCo2_Patho~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo2_Patho))
plot((PCo2_Patho), add.smooth = FALSE, which = 1)
res.PCo2_Patho<-residuals(PCo2_Patho)
bartlett.test(res.PCo2_Patho, dataF$Specie) #ok
bartlett.test(res.PCo2_Patho, dataF$Tto) #Yes
shapiro.test(res.PCo2_Patho) #ok
qqnorm(res.PCo2_Patho)

anova(PCo2_Patho) 

vf1 <- varIdent(form= ~1|Specie)
Mgls.PCo2_Patho <- gls(PCo2_Patho~Specie*Tto, data =dataF,weights = vf1)  
plot(Mgls.PCo2_Patho, add.smooth=FALSE, which=1)              
E.PCo2_Patho<-resid(Mgls.PCo2_Patho, type="normalized") 
bartlett.test(E.PCo2_Patho, dataF$Specie)  #ok  
bartlett.test(E.PCo2_Patho, dataF$Tto)    #Ok      
shapiro.test(E.PCo2_Patho) #ok
qqnorm(E.PCo2_Patho)

anova(Mgls.PCo2_Patho) 

ggplot(dataF, aes(x=Specie, y=S_Patho, fill=Tto)) + geom_boxplot()

#Sapro composit
PCo1_Sapro <-lm(PCo1_Sapro~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo1_Sapro))
plot((PCo1_Sapro), add.smooth = FALSE, which = 1)
res.PCo1_Sapro<-residuals(PCo1_Sapro)
bartlett.test(res.PCo1_Sapro, dataF$Specie) #no
bartlett.test(res.PCo1_Sapro, dataF$Tto) #no

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.PCo1_Sapro <- gls(PCo1_Sapro~Specie*Tto, data =dataF,weights = vf1,vf2)  
plot(Mgls.PCo1_Sapro, add.smooth=FALSE, which=1)              
E.PCo1_Sapro<-resid(Mgls.PCo1_Sapro, type="normalized") 
bartlett.test(E.PCo1_Sapro, dataF$Specie)  #ok  
bartlett.test(E.PCo1_Sapro, dataF$Tto)    #Ok      
shapiro.test(E.PCo1_Sapro) #ok
qqnorm(E.PCo1_Sapro)

anova(Mgls.PCo1_Sapro) 
#PCo2
PCo2_Sapro <-lm(PCo2_Sapro~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo2_Sapro))
plot((PCo2_Sapro), add.smooth = FALSE, which = 1)
res.PCo2_Sapro<-residuals(PCo2_Sapro)
bartlett.test(res.PCo2_Sapro, dataF$Specie) #no
bartlett.test(res.PCo2_Sapro, dataF$Tto) #no

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.PCo2_Sapro <- gls(PCo2_Sapro~Specie*Tto, data =dataF,weights = vf1,vf2)  
plot(Mgls.PCo2_Sapro, add.smooth=FALSE, which=1)              
E.PCo2_Sapro<-resid(Mgls.PCo2_Sapro, type="normalized") 
bartlett.test(E.PCo2_Sapro, dataF$Specie)  #ok  
bartlett.test(E.PCo2_Sapro, dataF$Tto)    #Ok      
shapiro.test(E.PCo2_Sapro) #ok
qqnorm(E.PCo2_Sapro)

anova(Mgls.PCo2_Sapro) 
ggplot(dataF, aes(x=Specie, y=S_Total, fill=Tto)) + geom_boxplot()

#Mutualists composit
PCo1_Mutual <-lm(PCo1_Mutual~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo1_Mutual))
plot((PCo1_Mutual), add.smooth = FALSE, which = 1)
res.PCo1_Mutual<-residuals(PCo1_Mutual)
bartlett.test(res.PCo1_Mutual, dataF$Specie) #no
bartlett.test(res.PCo1_Mutual, dataF$Tto) #no

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.PCo1_Mutual <- gls(PCo1_Mutual~Specie*Tto, data =dataF,weights = vf1,vf2)  
plot(Mgls.PCo1_Mutual, add.smooth=FALSE, which=1)              
E.PCo1_Mutual<-resid(Mgls.PCo1_Mutual, type="normalized") 
bartlett.test(E.PCo1_Mutual, dataF$Specie)  #ok  
bartlett.test(E.PCo1_Mutual, dataF$Tto)    #Ok      
shapiro.test(E.PCo1_Mutual) #ok
qqnorm(E.PCo1_Mutual)

anova(Mgls.PCo1_Mutual) 
#PCo2
PCo2_Mutual <-lm(PCo2_Mutual~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo2_Mutual))
plot((PCo2_Mutual), add.smooth = FALSE, which = 1)
res.PCo2_Mutual<-residuals(PCo2_Mutual)
bartlett.test(res.PCo2_Mutual, dataF$Specie) #no
bartlett.test(res.PCo2_Mutual, dataF$Tto) #ok

vf1 <- varIdent(form= ~1|Specie)
Mgls.PCo2_Mutual <- gls(PCo2_Mutual~Specie*Tto, data =dataF,weights = vf1)  
plot(Mgls.PCo2_Mutual, add.smooth=FALSE, which=1)              
E.PCo2_Mutual<-resid(Mgls.PCo2_Mutual, type="normalized") 
bartlett.test(E.PCo2_Mutual, dataF$Specie)  #ok  
bartlett.test(E.PCo2_Mutual, dataF$Tto)    #Ok      
shapiro.test(E.PCo2_Mutual) #ok
qqnorm(E.PCo2_Mutual)

anova(Mgls.PCo2_Mutual) 

ggplot(dataF, aes(x=Specie, y=S_Total, fill=Tto)) + geom_boxplot()


# ALL composition (Include also the Non asigned )
PCo1_ALL <-lm(PCo1_ALL~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo1_ALL))
plot((PCo1_ALL), add.smooth = FALSE, which = 1)
res.PCo1_ALL<-residuals(PCo1_ALL)
bartlett.test(res.PCo1_ALL, dataF$Specie) #no
bartlett.test(res.PCo1_ALL, dataF$Tto) #no

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.PCo1_ALL <- gls(PCo1_ALL~Specie*Tto, data =dataF,weights = vf1,vf2)  
plot(Mgls.PCo1_ALL, add.smooth=FALSE, which=1)              
E.PCo1_ALL<-resid(Mgls.PCo1_ALL, type="normalized") 
bartlett.test(E.PCo1_ALL, dataF$Specie)  #ok  
bartlett.test(E.PCo1_ALL, dataF$Tto)    #Ok      
shapiro.test(E.PCo1_ALL) #ok
qqnorm(E.PCo1_ALL)

anova(Mgls.PCo1_ALL) 

#PCo2
PCo2_ALL <-lm(PCo2_ALL~Specie*Tto, data=dataF) 
qqnorm(residuals(PCo2_ALL))
plot((PCo2_ALL), add.smooth = FALSE, which = 1)
res.PCo2_ALL<-residuals(PCo2_ALL)
bartlett.test(res.PCo2_ALL, dataF$Specie) #no
bartlett.test(res.PCo2_ALL, dataF$Tto) #no

vf1 <- varIdent(form= ~1|Specie)
vf2 <- varIdent(form= ~1|Tto)
Mgls.PCo2_ALL <- gls(PCo2_ALL~Specie*Tto, data =dataF,weights = vf1,vf2)  
plot(Mgls.PCo2_ALL, add.smooth=FALSE, which=1)              
E.PCo2_ALL<-resid(Mgls.PCo2_ALL, type="normalized") 
bartlett.test(E.PCo2_ALL, dataF$Specie)  #ok  
bartlett.test(E.PCo2_ALL, dataF$Tto)    #Ok      
shapiro.test(E.PCo2_ALL) #ok
qqnorm(E.PCo2_ALL)

anova(Mgls.PCo2_ALL) 

ggplot(dataF, aes(x=Specie, y=S_ALL, fill=Tto)) + geom_boxplot()

################################################## 

#################  CONTRASIES 

dataF$log_Pathogen <- log(dataF$Pathogens2+5)
dataF$log_Sapro <- log(dataF$Saprotrophs2+5)
dataF$log_Mutual <- log(dataF$Mutualists2+5)

#seleccionar cada species 
anth<- dataF[which(dataF$Specie=='G_Anthoxanthum'),]
berte<- dataF[which(dataF$Specie=='H_Berteroa'),]
dauc<- dataF[which(dataF$Specie=='H_Daucus'),]
fbrev<- dataF[which(dataF$Specie=='G_F.brevipila'),]
achi<- dataF[which(dataF$Specie=='H_Achillea'),]
pote<- dataF[which(dataF$Specie=='H_Potentilla'),]
plan<- dataF[which(dataF$Specie=='H_Plantago'),]
hiera<- dataF[which(dataF$Specie=='H_Hieracium'),]
artem<- dataF[which(dataF$Specie=='H_Artemisia'),]
holc<- dataF[which(dataF$Specie=='G_Holcus'),]
arrhe<- dataF[which(dataF$Specie=='G_Arrhenatherum'),]
vici<- dataF[which(dataF$Specie=='L_Vicia'),]
hyp<- dataF[which(dataF$Specie=='H_Hypericum'),]
frub<- dataF[which(dataF$Specie=='G_F.rubra'),]
gal<- dataF[which(dataF$Specie=='H_Galium'),]
poac<- dataF[which(dataF$Specie=='G_Poa'),]
sil<- dataF[which(dataF$Specie=='H_Silene'),]
trif<- dataF[which(dataF$Specie=='L_Trifolium'),]
dact<- dataF[which(dataF$Specie=='G_Dactylis'),]
rum<- dataF[which(dataF$Specie=='H_Rumex'),]
med<- dataF[which(dataF$Specie=='L_Medicago'),]
lol<- dataF[which(dataF$Specie=='G_Lolium'),]
ran<- dataF[which(dataF$Specie=='H_Ranunculus'),]
arm<- dataF[which(dataF$Specie=='H_Armeria'),]

# Root : shoot (contrastes)
anth_RS_aov <- aov(rootshoot_ratio~Tto, data= anth)
anth_RS_mc <- glht(anth_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) #POrque atras con el gls vi que necesito correcion heteroc
summary(anth_RS_mc)

berte_RS_aov <- aov(rootshoot_ratio~Tto, data= berte)
berte_RS_mc <- glht(berte_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_RS_mc)

dauc_RS_aov <- aov(rootshoot_ratio~Tto, data= dauc)
dauc_RS_mc <- glht(dauc_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_RS_mc)

fbrev_RS_aov <- aov(rootshoot_ratio~Tto, data= fbrev)
fbrev_RS_mc <- glht(fbrev_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_RS_mc)

achi_RS_aov <- aov(rootshoot_ratio~Tto, data= achi)
achi_RS_mc <- glht(achi_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_RS_mc)

pote_RS_aov <- aov(rootshoot_ratio~Tto, data= pote)
pote_RS_mc <- glht(pote_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_RS_mc)

plan_RS_aov <- aov(rootshoot_ratio~Tto, data= plan)
plan_RS_mc <- glht(plan_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_RS_mc)

hiera_RS_aov <- aov(rootshoot_ratio~Tto, data= hiera)
hiera_RS_mc <- glht(hiera_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_RS_mc)

artem_RS_aov <- aov(rootshoot_ratio~Tto, data= artem)
artem_RS_mc <- glht(artem_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_RS_mc)

holc_RS_aov <- aov(rootshoot_ratio~Tto, data= holc)
holc_RS_mc <- glht(holc_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_RS_mc)

arrhe_RS_aov <- aov(rootshoot_ratio~Tto, data= arrhe)
arrhe_RS_mc <- glht(arrhe_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_RS_mc)

vici_RS_aov <- aov(rootshoot_ratio~Tto, data= vici)
vici_RS_mc <- glht(vici_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_RS_mc)

hyp_RS_aov <- aov(rootshoot_ratio~Tto, data= hyp)
hyp_RS_mc <- glht(hyp_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_RS_mc)

frub_RS_aov <- aov(rootshoot_ratio~Tto, data= frub)
frub_RS_mc <- glht(frub_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_RS_mc)

gal_RS_aov <- aov(rootshoot_ratio~Tto, data= gal)
gal_RS_mc <- glht(gal_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_RS_mc)

poac_RS_aov <- aov(rootshoot_ratio~Tto, data= poac)
poac_RS_mc <- glht(poac_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_RS_mc)

sil_RS_aov <- aov(rootshoot_ratio~Tto, data= sil)
sil_RS_mc <- glht(sil_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_RS_mc)

trif_RS_aov <- aov(rootshoot_ratio~Tto, data= trif)
trif_RS_mc <- glht(trif_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_RS_mc)

dact_RS_aov <- aov(rootshoot_ratio~Tto, data= dact)
dact_RS_mc <- glht(dact_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_RS_mc)

rum_RS_aov <- aov(rootshoot_ratio~Tto, data= rum)
rum_RS_mc <- glht(rum_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_RS_mc)

med_RS_aov <- aov(rootshoot_ratio~Tto, data= med)
med_RS_mc <- glht(med_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_RS_mc)

lol_RS_aov <- aov(rootshoot_ratio~Tto, data= lol)
lol_RS_mc <- glht(lol_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_RS_mc)

ran_RS_aov <- aov(rootshoot_ratio~Tto, data= ran)
ran_RS_mc <- glht(ran_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_RS_mc)

arm_RS_aov <- aov(rootshoot_ratio~Tto, data= arm)
arm_RS_mc <- glht(arm_RS_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_RS_mc)



# Richnnes ALL contrastes
anth_RIC_aov <- aov(S_ALL~Tto, data= anth)
anth_RIC_mc <- glht(anth_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) #POrque atras con el gls vi que necesito correcion heteroc
summary(anth_RIC_mc)

berte_RIC_aov <- aov(S_ALL~Tto, data= berte)
berte_RIC_mc <- glht(berte_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_RIC_mc)

dauc_RIC_aov <- aov(S_ALL~Tto, data= dauc)
dauc_RIC_mc <- glht(dauc_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_RIC_mc)

fbrev_RIC_aov <- aov(S_ALL~Tto, data= fbrev)
fbrev_RIC_mc <- glht(fbrev_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_RIC_mc)

achi_RIC_aov <- aov(S_ALL~Tto, data= achi)
achi_RIC_mc <- glht(achi_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_RIC_mc)

pote_RIC_aov <- aov(S_ALL~Tto, data= pote)
pote_RIC_mc <- glht(pote_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_RIC_mc)

plan_RIC_aov <- aov(S_ALL~Tto, data= plan)
plan_RIC_mc <- glht(plan_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_RIC_mc)

hiera_RIC_aov <- aov(S_ALL~Tto, data= hiera)
hiera_RIC_mc <- glht(hiera_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_RIC_mc)

artem_RIC_aov <- aov(S_ALL~Tto, data= artem)
artem_RIC_mc <- glht(artem_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_RIC_mc)

holc_RIC_aov <- aov(S_ALL~Tto, data= holc)
holc_RIC_mc <- glht(holc_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_RIC_mc)

arrhe_RIC_aov <- aov(S_ALL~Tto, data= arrhe)
arrhe_RIC_mc <- glht(arrhe_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_RIC_mc)

vici_RIC_aov <- aov(S_ALL~Tto, data= vici)
vici_RIC_mc <- glht(vici_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_RIC_mc)

hyp_RIC_aov <- aov(S_ALL~Tto, data= hyp)
hyp_RIC_mc <- glht(hyp_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_RIC_mc)

frub_RIC_aov <- aov(S_ALL~Tto, data= frub)
frub_RIC_mc <- glht(frub_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_RIC_mc)

gal_RIC_aov <- aov(S_ALL~Tto, data= gal)
gal_RIC_mc <- glht(gal_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_RIC_mc)

poac_RIC_aov <- aov(S_ALL~Tto, data= poac)
poac_RIC_mc <- glht(poac_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_RIC_mc)

sil_RIC_aov <- aov(S_ALL~Tto, data= sil)
sil_RIC_mc <- glht(sil_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_RIC_mc)

trif_RIC_aov <- aov(S_ALL~Tto, data= trif)
trif_RIC_mc <- glht(trif_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_RIC_mc)

dact_RIC_aov <- aov(S_ALL~Tto, data= dact)
dact_RIC_mc <- glht(dact_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_RIC_mc)

rum_RIC_aov <- aov(S_ALL~Tto, data= rum)
rum_RIC_mc <- glht(rum_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_RIC_mc)

med_RIC_aov <- aov(S_ALL~Tto, data= med)
med_RIC_mc <- glht(med_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_RIC_mc)

lol_RIC_aov <- aov(S_ALL~Tto, data= lol)
lol_RIC_mc <- glht(lol_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_RIC_mc)

ran_RIC_aov <- aov(S_ALL~Tto, data= ran)
ran_RIC_mc <- glht(ran_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_RIC_mc)

arm_RIC_aov <- aov(S_ALL~Tto, data= arm)
arm_RIC_mc <- glht(arm_RIC_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_RIC_mc)

#Richness pathogen Contrastes
anth_RIC_PATHO_aov <- aov(S_Patho~Tto, data= anth)
anth_RIC_PATHO_mc <- glht(anth_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) #POrque atras con el gls vi que necesito correcion heteroc
summary(anth_RIC_PATHO_mc)

berte_RIC_PATHO_aov <- aov(S_Patho~Tto, data= berte)
berte_RIC_PATHO_mc <- glht(berte_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_RIC_PATHO_mc)

dauc_RIC_PATHO_aov <- aov(S_Patho~Tto, data= dauc)
dauc_RIC_PATHO_mc <- glht(dauc_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_RIC_PATHO_mc)

fbrev_RIC_PATHO_aov <- aov(S_Patho~Tto, data= fbrev)
fbrev_RIC_PATHO_mc <- glht(fbrev_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_RIC_PATHO_mc)

achi_RIC_PATHO_aov <- aov(S_Patho~Tto, data= achi)
achi_RIC_PATHO_mc <- glht(achi_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_RIC_PATHO_mc)

pote_RIC_PATHO_aov <- aov(S_Patho~Tto, data= pote)
pote_RIC_PATHO_mc <- glht(pote_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_RIC_PATHO_mc)

plan_RIC_PATHO_aov <- aov(S_Patho~Tto, data= plan)
plan_RIC_PATHO_mc <- glht(plan_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_RIC_PATHO_mc)

hiera_RIC_PATHO_aov <- aov(S_Patho~Tto, data= hiera)
hiera_RIC_PATHO_mc <- glht(hiera_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_RIC_PATHO_mc)

artem_RIC_PATHO_aov <- aov(S_Patho~Tto, data= artem)
artem_RIC_PATHO_mc <- glht(artem_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_RIC_PATHO_mc)

holc_RIC_PATHO_aov <- aov(S_Patho~Tto, data= holc)
holc_RIC_PATHO_mc <- glht(holc_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_RIC_PATHO_mc)

arrhe_RIC_PATHO_aov <- aov(S_Patho~Tto, data= arrhe)
arrhe_RIC_PATHO_mc <- glht(arrhe_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_RIC_PATHO_mc)

vici_RIC_PATHO_aov <- aov(S_Patho~Tto, data= vici)
vici_RIC_PATHO_mc <- glht(vici_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_RIC_PATHO_mc)

hyp_RIC_PATHO_aov <- aov(S_Patho~Tto, data= hyp)
hyp_RIC_PATHO_mc <- glht(hyp_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_RIC_PATHO_mc)

frub_RIC_PATHO_aov <- aov(S_Patho~Tto, data= frub)
frub_RIC_PATHO_mc <- glht(frub_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_RIC_PATHO_mc)

gal_RIC_PATHO_aov <- aov(S_Patho~Tto, data= gal)
gal_RIC_PATHO_mc <- glht(gal_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_RIC_PATHO_mc)

poac_RIC_PATHO_aov <- aov(S_Patho~Tto, data= poac)
poac_RIC_PATHO_mc <- glht(poac_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_RIC_PATHO_mc)

sil_RIC_PATHO_aov <- aov(S_Patho~Tto, data= sil)
sil_RIC_PATHO_mc <- glht(sil_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_RIC_PATHO_mc)

trif_RIC_PATHO_aov <- aov(S_Patho~Tto, data= trif)
trif_RIC_PATHO_mc <- glht(trif_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_RIC_PATHO_mc)

dact_RIC_PATHO_aov <- aov(S_Patho~Tto, data= dact)
dact_RIC_PATHO_mc <- glht(dact_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_RIC_PATHO_mc)

rum_RIC_PATHO_aov <- aov(S_Patho~Tto, data= rum)
rum_RIC_PATHO_mc <- glht(rum_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_RIC_PATHO_mc)

med_RIC_PATHO_aov <- aov(S_Patho~Tto, data= med)
med_RIC_PATHO_mc <- glht(med_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_RIC_PATHO_mc)

lol_RIC_PATHO_aov <- aov(S_Patho~Tto, data= lol)
lol_RIC_PATHO_mc <- glht(lol_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_RIC_PATHO_mc)

ran_RIC_PATHO_aov <- aov(S_Patho~Tto, data= ran)
ran_RIC_PATHO_mc <- glht(ran_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_RIC_PATHO_mc)

arm_RIC_PATHO_aov <- aov(S_Patho~Tto, data= arm)
arm_RIC_PATHO_mc <- glht(arm_RIC_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_RIC_PATHO_mc)

#Richness saprotrophs contrastes
anth_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= anth)
anth_RIC_SAPRO_mc <- glht(anth_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) #POrque atras con el gls vi que necesito correcion heteroc
summary(anth_RIC_SAPRO_mc)

berte_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= berte)
berte_RIC_SAPRO_mc <- glht(berte_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_RIC_SAPRO_mc)

dauc_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= dauc)
dauc_RIC_SAPRO_mc <- glht(dauc_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_RIC_SAPRO_mc)

fbrev_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= fbrev)
fbrev_RIC_SAPRO_mc <- glht(fbrev_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_RIC_SAPRO_mc)

achi_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= achi)
achi_RIC_SAPRO_mc <- glht(achi_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_RIC_SAPRO_mc)

pote_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= pote)
pote_RIC_SAPRO_mc <- glht(pote_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_RIC_SAPRO_mc)

plan_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= plan)
plan_RIC_SAPRO_mc <- glht(plan_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_RIC_SAPRO_mc)

hiera_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= hiera)
hiera_RIC_SAPRO_mc <- glht(hiera_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_RIC_SAPRO_mc)

artem_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= artem)
artem_RIC_SAPRO_mc <- glht(artem_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_RIC_SAPRO_mc)

holc_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= holc)
holc_RIC_SAPRO_mc <- glht(holc_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_RIC_SAPRO_mc)

arrhe_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= arrhe)
arrhe_RIC_SAPRO_mc <- glht(arrhe_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_RIC_SAPRO_mc)

vici_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= vici)
vici_RIC_SAPRO_mc <- glht(vici_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_RIC_SAPRO_mc)

hyp_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= hyp)
hyp_RIC_SAPRO_mc <- glht(hyp_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_RIC_SAPRO_mc)

frub_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= frub)
frub_RIC_SAPRO_mc <- glht(frub_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_RIC_SAPRO_mc)

gal_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= gal)
gal_RIC_SAPRO_mc <- glht(gal_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_RIC_SAPRO_mc)

poac_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= poac)
poac_RIC_SAPRO_mc <- glht(poac_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_RIC_SAPRO_mc)

sil_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= sil)
sil_RIC_SAPRO_mc <- glht(sil_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_RIC_SAPRO_mc)

trif_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= trif)
trif_RIC_SAPRO_mc <- glht(trif_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_RIC_SAPRO_mc)

dact_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= dact)
dact_RIC_SAPRO_mc <- glht(dact_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_RIC_SAPRO_mc)

rum_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= rum)
rum_RIC_SAPRO_mc <- glht(rum_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_RIC_SAPRO_mc)

med_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= med)
med_RIC_SAPRO_mc <- glht(med_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_RIC_SAPRO_mc)

lol_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= lol)
lol_RIC_SAPRO_mc <- glht(lol_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_RIC_SAPRO_mc)

ran_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= ran)
ran_RIC_SAPRO_mc <- glht(ran_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_RIC_SAPRO_mc)

arm_RIC_SAPRO_aov <- aov(S_Sapro~Tto, data= arm)
arm_RIC_SAPRO_mc <- glht(arm_RIC_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_RIC_SAPRO_mc)

# Richness mutualists contrastes

anth_RIC_MUT_aov <- aov(S_Mutual~Tto, data= anth)
anth_RIC_MUT_mc <- glht(anth_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(anth_RIC_MUT_mc)

berte_RIC_MUT_aov <- aov(S_Mutual~Tto, data= berte)
berte_RIC_MUT_mc <- glht(berte_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_RIC_MUT_mc)

dauc_RIC_MUT_aov <- aov(S_Mutual~Tto, data= dauc)
dauc_RIC_MUT_mc <- glht(dauc_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_RIC_MUT_mc)

fbrev_RIC_MUT_aov <- aov(S_Mutual~Tto, data= fbrev)
fbrev_RIC_MUT_mc <- glht(fbrev_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_RIC_MUT_mc)

achi_RIC_MUT_aov <- aov(S_Mutual~Tto, data= achi)
achi_RIC_MUT_mc <- glht(achi_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_RIC_MUT_mc)

pote_RIC_MUT_aov <- aov(S_Mutual~Tto, data= pote)
pote_RIC_MUT_mc <- glht(pote_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_RIC_MUT_mc)

plan_RIC_MUT_aov <- aov(S_Mutual~Tto, data= plan)
plan_RIC_MUT_mc <- glht(plan_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_RIC_MUT_mc)

hiera_RIC_MUT_aov <- aov(S_Mutual~Tto, data= hiera)
hiera_RIC_MUT_mc <- glht(hiera_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_RIC_MUT_mc)

artem_RIC_MUT_aov <- aov(S_Mutual~Tto, data= artem)
artem_RIC_MUT_mc <- glht(artem_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_RIC_MUT_mc)

holc_RIC_MUT_aov <- aov(S_Mutual~Tto, data= holc)
holc_RIC_MUT_mc <- glht(holc_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_RIC_MUT_mc)

arrhe_RIC_MUT_aov <- aov(S_Mutual~Tto, data= arrhe)
arrhe_RIC_MUT_mc <- glht(arrhe_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_RIC_MUT_mc)

vici_RIC_MUT_aov <- aov(S_Mutual~Tto, data= vici)
vici_RIC_MUT_mc <- glht(vici_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_RIC_MUT_mc)

hyp_RIC_MUT_aov <- aov(S_Mutual~Tto, data= hyp)
hyp_RIC_MUT_mc <- glht(hyp_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_RIC_MUT_mc)

frub_RIC_MUT_aov <- aov(S_Mutual~Tto, data= frub)
frub_RIC_MUT_mc <- glht(frub_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_RIC_MUT_mc)

gal_RIC_MUT_aov <- aov(S_Mutual~Tto, data= gal)
gal_RIC_MUT_mc <- glht(gal_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_RIC_MUT_mc)

poac_RIC_MUT_aov <- aov(S_Mutual~Tto, data= poac)
poac_RIC_MUT_mc <- glht(poac_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_RIC_MUT_mc)

sil_RIC_MUT_aov <- aov(S_Mutual~Tto, data= sil)
sil_RIC_MUT_mc <- glht(sil_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_RIC_MUT_mc)

trif_RIC_MUT_aov <- aov(S_Mutual~Tto, data= trif)
trif_RIC_MUT_mc <- glht(trif_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_RIC_MUT_mc)

dact_RIC_MUT_aov <- aov(S_Mutual~Tto, data= dact)
dact_RIC_MUT_mc <- glht(dact_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_RIC_MUT_mc)

rum_RIC_MUT_aov <- aov(S_Mutual~Tto, data= rum)
rum_RIC_MUT_mc <- glht(rum_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_RIC_MUT_mc)

med_RIC_MUT_aov <- aov(S_Mutual~Tto, data= med)
med_RIC_MUT_mc <- glht(med_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_RIC_MUT_mc)

lol_RIC_MUT_aov <- aov(S_Mutual~Tto, data= lol)
lol_RIC_MUT_mc <- glht(lol_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_RIC_MUT_mc)

ran_RIC_MUT_aov <- aov(S_Mutual~Tto, data= ran)
ran_RIC_MUT_mc <- glht(ran_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_RIC_MUT_mc)

arm_RIC_MUT_aov <- aov(S_Mutual~Tto, data= arm)
arm_RIC_MUT_mc <- glht(arm_RIC_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_RIC_MUT_mc)

# Abundance Pathogens Contrastes
anth_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= anth)
anth_AB_PATHO_mc <- glht(anth_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) #POrque atras con el gls vi que necesito correcion heteroc
summary(anth_AB_PATHO_mc)

berte_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= berte)
berte_AB_PATHO_mc <- glht(berte_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_AB_PATHO_mc)

dauc_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= dauc)
dauc_AB_PATHO_mc <- glht(dauc_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_AB_PATHO_mc)

fbrev_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= fbrev)
fbrev_AB_PATHO_mc <- glht(fbrev_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_AB_PATHO_mc)

achi_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= achi)
achi_AB_PATHO_mc <- glht(achi_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_AB_PATHO_mc)

pote_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= pote)
pote_AB_PATHO_mc <- glht(pote_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_AB_PATHO_mc)

plan_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= plan)
plan_AB_PATHO_mc <- glht(plan_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_AB_PATHO_mc)

hiera_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= hiera)
hiera_AB_PATHO_mc <- glht(hiera_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_AB_PATHO_mc)

artem_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= artem)
artem_AB_PATHO_mc <- glht(artem_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_AB_PATHO_mc)

holc_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= holc)
holc_AB_PATHO_mc <- glht(holc_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_AB_PATHO_mc)

arrhe_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= arrhe)
arrhe_AB_PATHO_mc <- glht(arrhe_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_AB_PATHO_mc)

vici_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= vici)
vici_AB_PATHO_mc <- glht(vici_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_AB_PATHO_mc)

hyp_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= hyp)
hyp_AB_PATHO_mc <- glht(hyp_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_AB_PATHO_mc)

frub_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= frub)
frub_AB_PATHO_mc <- glht(frub_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_AB_PATHO_mc)

gal_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= gal)
gal_AB_PATHO_mc <- glht(gal_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_AB_PATHO_mc)

poac_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= poac)
poac_AB_PATHO_mc <- glht(poac_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_AB_PATHO_mc)

sil_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= sil)
sil_AB_PATHO_mc <- glht(sil_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_AB_PATHO_mc)

trif_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= trif)
trif_AB_PATHO_mc <- glht(trif_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_AB_PATHO_mc)

dact_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= dact)
dact_AB_PATHO_mc <- glht(dact_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_AB_PATHO_mc)

rum_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= rum)
rum_AB_PATHO_mc <- glht(rum_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_AB_PATHO_mc)

med_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= med)
med_AB_PATHO_mc <- glht(med_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_AB_PATHO_mc)

lol_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= lol)
lol_AB_PATHO_mc <- glht(lol_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_AB_PATHO_mc)

ran_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= ran)
ran_AB_PATHO_mc <- glht(ran_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_AB_PATHO_mc)

arm_AB_PATHO_aov <- aov(log_Pathogen~Tto, data= arm)
arm_AB_PATHO_mc <- glht(arm_AB_PATHO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_AB_PATHO_mc)

# Abundance Saprotrophs Contrastes
anth_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= anth)
anth_AB_SAPRO_mc <- glht(anth_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) #POrque atras con el gls vi que necesito correcion heteroc
summary(anth_AB_SAPRO_mc)

berte_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= berte)
berte_AB_SAPRO_mc <- glht(berte_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_AB_SAPRO_mc)

dauc_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= dauc)
dauc_AB_SAPRO_mc <- glht(dauc_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_AB_SAPRO_mc)

fbrev_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= fbrev)
fbrev_AB_SAPRO_mc <- glht(fbrev_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_AB_SAPRO_mc)

achi_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= achi)
achi_AB_SAPRO_mc <- glht(achi_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_AB_SAPRO_mc)

pote_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= pote)
pote_AB_SAPRO_mc <- glht(pote_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_AB_SAPRO_mc)

plan_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= plan)
plan_AB_SAPRO_mc <- glht(plan_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_AB_SAPRO_mc)

hiera_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= hiera)
hiera_AB_SAPRO_mc <- glht(hiera_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_AB_SAPRO_mc)

artem_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= artem)
artem_AB_SAPRO_mc <- glht(artem_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_AB_SAPRO_mc)

holc_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= holc)
holc_AB_SAPRO_mc <- glht(holc_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_AB_SAPRO_mc)

arrhe_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= arrhe)
arrhe_AB_SAPRO_mc <- glht(arrhe_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_AB_SAPRO_mc)

vici_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= vici)
vici_AB_SAPRO_mc <- glht(vici_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_AB_SAPRO_mc)

hyp_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= hyp)
hyp_AB_SAPRO_mc <- glht(hyp_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_AB_SAPRO_mc)

frub_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= frub)
frub_AB_SAPRO_mc <- glht(frub_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_AB_SAPRO_mc)

gal_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= gal)
gal_AB_SAPRO_mc <- glht(gal_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_AB_SAPRO_mc)

poac_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= poac)
poac_AB_SAPRO_mc <- glht(poac_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_AB_SAPRO_mc)

sil_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= sil)
sil_AB_SAPRO_mc <- glht(sil_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_AB_SAPRO_mc)

trif_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= trif)
trif_AB_SAPRO_mc <- glht(trif_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_AB_SAPRO_mc)

dact_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= dact)
dact_AB_SAPRO_mc <- glht(dact_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_AB_SAPRO_mc)

rum_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= rum)
rum_AB_SAPRO_mc <- glht(rum_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_AB_SAPRO_mc)

med_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= med)
med_AB_SAPRO_mc <- glht(med_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_AB_SAPRO_mc)

lol_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= lol)
lol_AB_SAPRO_mc <- glht(lol_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_AB_SAPRO_mc)

ran_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= ran)
ran_AB_SAPRO_mc <- glht(ran_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_AB_SAPRO_mc)

arm_AB_SAPRO_aov <- aov(log_Sapro~Tto, data= arm)
arm_AB_SAPRO_mc <- glht(arm_AB_SAPRO_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_AB_SAPRO_mc)

# Abundance Mutualists Contrastes
anth_AB_MUT_aov <- aov(log_Mutual~Tto, data= anth)
anth_AB_MUT_mc <- glht(anth_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(anth_AB_MUT_mc)

berte_AB_MUT_aov <- aov(log_Mutual~Tto, data= berte)
berte_AB_MUT_mc <- glht(berte_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(berte_AB_MUT_mc)

dauc_AB_MUT_aov <- aov(log_Mutual~Tto, data= dauc)
dauc_AB_MUT_mc <- glht(dauc_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dauc_AB_MUT_mc)

fbrev_AB_MUT_aov <- aov(log_Mutual~Tto, data= fbrev)
fbrev_AB_MUT_mc <- glht(fbrev_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(fbrev_AB_MUT_mc)

achi_AB_MUT_aov <- aov(log_Mutual~Tto, data= achi)
achi_AB_MUT_mc <- glht(achi_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(achi_AB_MUT_mc)

pote_AB_MUT_aov <- aov(log_Mutual~Tto, data= pote)
pote_AB_MUT_mc <- glht(pote_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(pote_AB_MUT_mc)

plan_AB_MUT_aov <- aov(log_Mutual~Tto, data= plan)
plan_AB_MUT_mc <- glht(plan_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(plan_AB_MUT_mc)

hiera_AB_MUT_aov <- aov(log_Mutual~Tto, data= hiera)
hiera_AB_MUT_mc <- glht(hiera_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hiera_AB_MUT_mc)

artem_AB_MUT_aov <- aov(log_Mutual~Tto, data= artem)
artem_AB_MUT_mc <- glht(artem_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(artem_AB_MUT_mc)

holc_AB_MUT_aov <- aov(log_Mutual~Tto, data= holc)
holc_AB_MUT_mc <- glht(holc_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(holc_AB_MUT_mc)

arrhe_AB_MUT_aov <- aov(log_Mutual~Tto, data= arrhe)
arrhe_AB_MUT_mc <- glht(arrhe_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arrhe_AB_MUT_mc)

vici_AB_MUT_aov <- aov(log_Mutual~Tto, data= vici)
vici_AB_MUT_mc <- glht(vici_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(vici_AB_MUT_mc)

hyp_AB_MUT_aov <- aov(log_Mutual~Tto, data= hyp)
hyp_AB_MUT_mc <- glht(hyp_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(hyp_AB_MUT_mc)

frub_AB_MUT_aov <- aov(log_Mutual~Tto, data= frub)
frub_AB_MUT_mc <- glht(frub_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(frub_AB_MUT_mc)

gal_AB_MUT_aov <- aov(log_Mutual~Tto, data= gal)
gal_AB_MUT_mc <- glht(gal_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(gal_AB_MUT_mc)

poac_AB_MUT_aov <- aov(log_Mutual~Tto, data= poac)
poac_AB_MUT_mc <- glht(poac_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(poac_AB_MUT_mc)

sil_AB_MUT_aov <- aov(log_Mutual~Tto, data= sil)
sil_AB_MUT_mc <- glht(sil_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(sil_AB_MUT_mc)

trif_AB_MUT_aov <- aov(log_Mutual~Tto, data= trif)
trif_AB_MUT_mc <- glht(trif_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(trif_AB_MUT_mc)

dact_AB_MUT_aov <- aov(log_Mutual~Tto, data= dact)
dact_AB_MUT_mc <- glht(dact_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(dact_AB_MUT_mc)

rum_AB_MUT_aov <- aov(log_Mutual~Tto, data= rum)
rum_AB_MUT_mc <- glht(rum_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(rum_AB_MUT_mc)

med_AB_MUT_aov <- aov(log_Mutual~Tto, data= med)
med_AB_MUT_mc <- glht(med_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(med_AB_MUT_mc)

lol_AB_MUT_aov <- aov(log_Mutual~Tto, data= lol)
lol_AB_MUT_mc <- glht(lol_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(lol_AB_MUT_mc)

ran_AB_MUT_aov <- aov(log_Mutual~Tto, data= ran)
ran_AB_MUT_mc <- glht(ran_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(ran_AB_MUT_mc)

arm_AB_MUT_aov <- aov(log_Mutual~Tto, data= arm)
arm_AB_MUT_mc <- glht(arm_AB_MUT_aov, linfct = mcp(Tto = "Tukey"), vcov=sandwich) 
summary(arm_AB_MUT_mc)

################

###########   Variation paritioning

dataF <- read.csv ("traits_fungal.csv")

#0. https://rdrr.io/rforge/vegan/man/varpart.html
## Negative values of Ra2 are interpreted as zeros; they correspond to cases where the explanatory variables 
#explain less variation than random normal variables would."
#ojo "b" que es la interseccion NO LA PUEDO TESTAR (no tiene degreess) See numerical ecology p230
#using RDA step by step is the same tha using varpart (ver ejemplo en Varpart_root)...pero ojo...igual q los R2 sin adjust
# If Y are dissimilarities, the decomposition is based on distance-based redundancy analysis (db-RDA, see capscale) from http://finzi.psych.upenn.edu/R/library/vegan/html/varpart.html
# Varpart Con diferent matrices (todas matrices) el orden IMPORTAAAAAAAA... 

#var_trial ### El codigo funciona PERFECTTTTT!!!!!   
trialvarp <- read.csv ("trialvarpartborrar.csv")
variab<-trialvarp[, grep(pattern="^ESV", colnames(trialvarp))] 
View(trialvarp)
#varpart usando "rda" lo hace paso por paso
rda (variab ~ trialvarp$specie + trialvarp$tto) # constrained = 0.94
rda(variab ~trialvarp$specie)                # constrained = 0.02  (proportion)
rda(variab ~trialvarp$tto)                       #constrained = 0.92
#varpart usando la funcion "varpart" sale todo de una vez 
varpart(variab, ~trialvarp$specie, ~trialvarp$tto) #      Rsq= 0.02, 0.92, 0.94 

#Varpart_many matrices
trialvarp2 <- read.csv ("trialvarpartborrar2.csv")
trialvarp3 <- read.csv ("trialvarpartborrar3.csv")

varpart(trialvarp2, ~ aa+bb+cc, variab, data=trialvarp3) # EL ORDEN IMPORTA!!!!!.. si pongo primero la segunda var (variab) no funciona 
varpart(trialvarp2, ~ aa+bb+cc, variab, trialvarp2, data=trialvarp3)  # FUNCIONA ok
varpart(variab, ~ aa+bb+cc, trialvarp2, data=trialvarp3)        ## FUNCIONA OK
varpart(variab, ~ trialvarp3, trialvarp2) ### No funciona..raro!!! 

#############################

#1. Pathogens, sapro, mutual (effect of drought and plant species on them)
#2. Total fungal, root traits, leaf traits (effect of drought and plant species on them)
getwd()

#Varpart_Pathogens
data_patESV <- read.csv("PAST_Pathogen_composition_NEW.csv")
patESV <- data_patESV [, grep(pattern="^esv", colnames(data_patESV))]  #seleccionar las ESV
patESV_hell<- decostand(patESV, method="hellinger")
varpart (patESV_hell, ~ data_patESV$species, ~data_patESV$tto)

anova.cca(rda(patESV_hell, dataF$Specie), step=1000)  #Test of fraction [a] #
anova.cca(rda(patESV_hell, dataF$Tto), step=1000)
#varpart _Pathogens  (PCo1 explains more than PCo2)
varpart(dataF$PCo2_Patho, ~dataF$Specie, ~dataF$Tto)
anova.cca(rda(dataF$PCo1_Sapro, dataF$Specie), step=1000)
anova.cca(rda(dataF$PCo1_Sapro, dataF$Tto), step=1000)

#Varpart_Saprotrophs
data_saprESV <- read.csv("PAST_Saprotro_Composition_NEW.csv")
saprESV <- data_saprESV [, grep(pattern="^esv", colnames(data_saprESV))]  #seleccionar las ESV
saprESV_hell<- decostand(saprESV, method="hellinger")
varpart (saprESV_hell, ~ data_saprESV$Species, ~data_saprESV$Tto)

anova.cca(rda(saprESV_hell, dataF$Specie), step=1000)  #Test of fraction [a] #
anova.cca(rda(saprESV_hell, dataF$Tto), step=1000)
#varpart _Saprotrophs  (PCo1 explains more than PCo2)
varpart(dataF$PCo1_Sapro, ~dataF$Specie, ~dataF$Tto)
anova.cca(rda(dataF$PCo1_Sapro, dataF$Specie), step=1000)
anova.cca(rda(dataF$PCo1_Sapro, dataF$Tto), step=1000)


#Varpart_Mutualists
data_mutESV <- read.csv("PAST_Mutualists_Composition_NEW.csv")  ### Le he quitado las EMF 
mutualESV <- data_mutESV [, grep(pattern="^esv", colnames(data_mutESV))] 
mutualESV_hell<- decostand(mutualESV, method="hellinger")
varpart (mutualESV_hell, ~ data_mutESV$Specie, ~data_mutESV$Tto)

anova.cca(rda(mutualESV_hell, dataF$Specie), step=1000)  #Test of fraction [a] #
anova.cca(rda(mutualESV_hell, dataF$Tto), step=1000)
#varpart _Mutual (PCo2) Porq da mejor que el PCo1 
varpart(dataF$PCo2_Mutual, ~dataF$Specie, ~dataF$Tto)
anova.cca(rda(dataF$PCo2_Mutual, dataF$Specie), step=1000)
anova.cca(rda(dataF$PCo2_Mutual, dataF$Tto), step=1000)

# Varpart_Total fungal
dataESV <- read.csv("ESV_totalFungComposition_mantel.csv")  # la tabla donde esta la "total" fungal composition. The ESV selected
str(dataESV)
#qutiar ESVs que eran EMF   #### NOOO porq tb son Fungi
#dataESV <-subset(dataESV, select = -c(esv1462,esv2447,esv2818,esv3183,esv2531,esv2604,esv1811,esv2013,
#                                              esv2603,esv2954,esv3286,esv3395,esv3396,esv3440)) # quitar unos ESV que son EMF (-c quitar,column)

totESV<-dataESV[, grep(pattern="^esv", colnames(dataESV))]  # seleccionar las columnas que empiezan con "esv" 
varpart(totESV, ~dataF$Specie, ~dataF$Tto)

totESV_hell<- decostand(totESV, method="hellinger") # transformando con Hellinger.. igual no es necesari porq todas son ESV ( no hay q standariz)
varpart(totESV_hell, ~dataF$Specie, ~dataF$Tto)

anova.cca(rda(totESV_hell, dataF$Specie), step=1000)  #Test of fraction [a] #shape#
anova.cca(rda(totESV_hell, dataF$Tto), step=1000)
#varpart _Fungal (PCo2) Da mejor que el PCo1..y es el q PCoA q abajo da mejor en fungal vs root vs leaf
varpart(dataF$PCo2_ALL, ~dataF$Specie, ~dataF$Tto)
anova.cca(rda(dataF$PCo2_ALL, dataF$Specie), step=1000)
anova.cca(rda(dataF$PCo2_ALL, dataF$Tto), step=1000)

#############################

#Varpart_Root traits

# he transformado las variables como en el de root traits-drought paper y da peor el varpart. Explica menos.

Root_var<- dataF[,c("sqrt_diameter","log_SRL","sqrt_SRSA","log_RTD","sqrt_root","log_rootN","log_rootC")]
vif <- diag(solve(cor(Root_var))) # colinealidad? VIFs > 10 should be avoided.  Numerical ecology
Root_var<- dataF[,c("sqrt_diameter","sqrt_SRSA","log_RTD","sqrt_root","log_rootN","log_rootC")] #quite SRL
vif <- diag(solve(cor(Root_var)))

Root_var_hell<- decostand(Root_var, method="hellinger", na.rm = TRUE)  #transformations may be applied PRIOR rda or pca (see GUSTAME)
#                                                               mejoro un poco q sin transform . Paso de 0.62 a 0.72 
#variation partitioning step by step using RDA (es lo mismo que using varpart)
rda (Root_var_hell ~ dataF$Specie + dataF$Tto) #constrained proportion (0.72)
rda(Root_var_hell ~dataF$Specie)               #                       (0.72)
rda(Root_var_hell ~dataF$Tto)                                         #(1.1 E-03)

#Using varpart 
varpart (Root_var_hell, ~dataF$Specie, ~dataF$Tto)  # los resultados son identicos q con RDA step by step

anova.cca(rda(Root_var_hell, dataF$Specie), step=1000)  #Test of fraction [a] #shape#
anova.cca(rda(Root_var_hell, dataF$Tto), step=1000) ##


#Varpart_Roottraits_PHYLOGENIA

rda(Root_var_hell ~ Specie*Tto+Condition(Phylog_1,Phylog_2,Phylog_3,Phylog_4), data=dataF, scale = TRUE)
rda(Root_var_hell ~dataF$Specie+Condition(Phylog_1,Phylog_2,Phylog_3,Phylog_4), data=dataF, scale = TRUE)   
rda(Root_var_hell ~dataF$Tto+Condition(Phylog_1,Phylog_2,Phylog_3,Phylog_4), data=dataF, scale = TRUE) 

varpart (Root_var_hell ~dataF$Specie, ~dataF$Tto | Condition(Phylog_1,Phylog_2,Phylog_3,Phylog_4) ) ## NO FUNCIONA..aun

##########################

#Varpart_Leaf traits
Leaf_var<- dataF[,c("shootmass","leafC","leafN")] #"LDMC","SLA",
Leaf_var_hell<- decostand(Leaf_var, method="hellinger", na.rm = TRUE) 
varpart (Leaf_var_hell, ~dataF$Specie, ~dataF$Tto)

anova.cca(rda(Leaf_var_hell, dataF$Specie), step=1000)  #Test of fraction [a] #
anova.cca(rda(Leaf_var_hell, dataF$Tto), step=1000)

#varpart_Soil traits
Soil_var <-dataF[,c("pH","soilC","soilN", "WSA")]
Soil_var_hell<- decostand(Soil_var, method="hellinger") 
varpart (Soil_var_hell, ~dataF$Specie, ~dataF$Tto)

anova.cca(rda(Soil_var_hell, dataF$Specie), step=1000)  #Test of fraction [a] #
anova.cca(rda(Soil_var_hell, dataF$Tto), step=1000)

##########

#Varpart _fungal matrix vs soil, root, leaf 

varpart(totESV_hell, Leaf_var_hell, Root_var_hell, Soil_var_hell)

anova.cca(rda(totESV_hell, Leaf_var_hell), step=1000)  #Test of fraction [a] #
anova.cca(rda(totESV_hell, Root_var_hell), step=1000)
anova.cca(rda(totESV_hell, Soil_var_hell), step=1000)

#varpart_fungal matrix vs root and leaf

varpart(totESV_hell, Leaf_var_hell, Root_var_hell)     # 0.013 Explicado por leaf, 0.012 por root traits
varpart(dataF$PCo1_ALL, Leaf_var_hell, Root_var_hell)  # 0.00311                  , 0.00019 

#varpart fungal using PCo2 
varpart(dataF$PCo2_ALL, Leaf_var_hell, Root_var_hell)  # 0.00114                   , 0.0204     #### ESCOGIDO!!! 
anova.cca(rda(dataF$PCo2_ALL, Leaf_var_hell), step=1000)
anova.cca(rda(dataF$PCo2_ALL, Root_var_hell), step=1000)
#
rda(Root_var_hell ~dataF$Tto+Condition(Phylog_1,Phylog_2,Phylog_3,Phylog_4), data=dataF, scale = TRUE) # no funciona..no se si tiene sentido...


########################
######################################################################
######################### PATH ANALYSES#############

str(dataF)

#Option 1. HAcer un stepAIC para esocger las variables y luego el path analyses. (como en soil blocks)

# 1. Select the variables que van a ir en el modelo 

# a)Pathogens..He puesto todas las de pathogenos para escoger las mejores predictores.

dataF <- as.data.frame(dataF)
M1path <- lm(log_shootmass ~ Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho, data=dataF)
stepAIC(M1path, direction = 'backward') # selected: PCo2 pathog

#M2path <- lm(CP1_leaft ~ Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho, data=dataF)
#stepAIC(M2path, direction = 'backward') # selected: PCo2 pathog 

#M3path <- lm(diameter ~ Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho, data=dataF)
#stepAIC(M3path, direction = 'backward') # selected: S_patho, log_S_patho

#M4path <- lm(SRSA ~ Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho, data=dataF)
#stepAIC(M4path, direction = 'backward') # selected: PCo2 patho, PCo1 patho, 

#M5path <- lm(CP1_root ~ Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho, data=dataF)
#stepAIC(M5path, direction = 'backward') # selected: log_Pathogen, Pathogen2, log_S_Patho

#Pathogens as response varible (PCo2 PAtho fue la escogida atras (log shootmass), y ahora a esta ver los root traits )

P1 <- lm(PCo2_Patho ~ log_rootN + log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(P1, direction = 'backward') #  sqrt_SRSA, root:shoot ratio , 

dataF <-as.data.frame (dataF)
rcorr(dataF, type="pearson")

str(dataF)
a<- dataF[,c("log_rootC", "log_rootN","sqrt_root","sqrt_diameter","log_RTD","log_SRL","sqrt_SRSA")]
a<- data.matrix(a)
b<- dataF[,c("Pathogens2","log_Pathogen", "S_Patho","log_S_Patho","PCo1_Patho","PCo2_Patho","log_S_Mutual")]
b<- data.matrix(b)
ab <-cor(a,b)
corrplot(ab)  # PCoA2 Patho es la q mejor se relaciona con most of the root traits ( CONFIRMADO... Aunq ya lo habia escogido atras en P1)

#b) saprotrophs
names(dataF)
M1sapro <- lm(log_shootmass ~ Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro, data=dataF)
stepAIC(M1sapro, direction = 'backward') # selected:  Saprotrophs2, S_sapro, PCo2. 

#M2sapro <- lm(CP1_leaft ~ Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro, data=dataF)
#stepAIC(M2sapro, direction = 'backward') # selected: S_sapro

#M3sapro <- lm(diameter ~ Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro, data=dataF)
#stepAIC(M3sapro, direction = 'backward') # selected: S_sapro, log_sapro

#M4sapro <- lm(SRSA ~ Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro, data=dataF)
#stepAIC(M4sapro, direction = 'backward') # selected: Saprotrophs2

#M5sapro <- lm(CP1_root ~ Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro, data=dataF)
#stepAIC(M5sapro, direction = 'backward') # selected: Saprotrophs2

S0 <- lm(Saprotrophs2 ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S0, direction = 'backward') #RTD, SRSA 

S4 <- lm(S_Sapro ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S4, direction = 'backward') # sqrt_SRSA, logSRL, logRTD

S1 <- lm(PCo1_Sapro ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S1, direction = 'backward') #log N y SRSA

S1a <- lm(PCo2_Sapro ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S1a, direction = 'backward') #log SRL

S2 <- lm(log_S_Sapro ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S2, direction = 'backward') # log_RTD, sqrtSRSA, SRL

S3 <- lm(PCo2_Sapro ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S3, direction = 'backward') # log_SRL

c <- dataF[,c("Saprotrophs2","log_Sapro","S_Sapro","log_S_Sapro", "PCo1_Sapro", "PCo2_Sapro")]
c<- data.matrix(c)
ac <- cor(a,c)
corrplot(ac)   # no es tan claro como con patho...pero si se CONFIRMA q escojo Sapro 2. q aparte de estar relacionado con Shoot, lo esta con Root traits

S2 <- lm(log_S_Sapro ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(S2, direction = 'backward') # srsa, srl, rtd,diameter


#c) Mutualists
M1Mutual <- lm(log_shootmass ~ Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual, data=dataF)
stepAIC(M1Mutual, direction = 'backward')  #selected:   PCo2 Mutual, log_S_Mutual

#M2Mutual <- lm(CP1_leaft ~ Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual, data=dataF)
#stepAIC(M2Mutual, direction = 'backward')  #selected: log_mutual, S_mutual, log_S_mutual

#M3Mutual <- lm(diameter ~ Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual, data=dataF)
#stepAIC(M3Mutual, direction = 'backward')  #selected: S_Mutual, PCo1 Mutual

#M4Mutual <- lm(SRSA ~ Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual, data=dataF)
#stepAIC(M4Mutual, direction = 'backward')  #selected: S_mutual, PCo1 Mutual

#M5Mutual <- lm(CP1_root ~ Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual, data=dataF)
#stepAIC(M5Mutual, direction = 'backward')  #selected: S_mutual

MUT1 <- lm(PCo2_Mutual ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF )
stepAIC(MUT1, direction = 'backward') # root:shoot 

MUT2 <- lm(log_S_Mutual ~ log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA+ log_rootN, data = dataF )
stepAIC(MUT2, direction = 'backward') # sqrt_diameter, log_rootN 



d <- dataF[,c("Mutualists2","log_Mutual","S_Mutual","log_S_Mutual","PCo1_Mutual","PCo2_Mutual")]
d<- data.matrix(d)
ad <- cor(a,d)
corrplot(ad)   ### Log_S_Mutual tambien es importante aqui... y ojo q para el stepAIC con shoot tb. 


#d) soil properties
names(dataF)
M1soil <- lm(log_shootmass ~ CP1_soilprop+CP2_soilprop,data=dataF)
stepAIC(M1soil, direction = 'backward')   # CP2 

a<- dataF[,c("log_rootC", "log_rootN","sqrt_root","sqrt_diameter","log_RTD","log_SRL","sqrt_SRSA")]
a<- data.matrix(a)
soil<- dataF[,c("CP1_soilprop","CP2_soilprop")]
soil<- data.matrix(soil)
asoil <-cor(a,soil)
corrplot(asoil) ## aqui tb CP2 ( esta bien relacionado con diameter)

#e) Root traits
M1root <- lm(log_shootmass ~ log_rootN + log_rootC+ sqrt_root+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data=dataF)
stepAIC(M1root, direction = 'backward')   # diameter, root mass, rootN, root C  #### HE ESCOGIDO ESTA 

######################################################
########### PATH ANALYSES

dataF$log_SRL100 <- dataF$log_SRL/100
dataF$CP2_soilprop100 <- dataF$CP2_soilprop/100
dataF$sqrt_root100 <- dataF$sqrt_root/100
dataF$rootshoot_ratio100 <- dataF$rootshoot_ratio/100
dataF$log_rootshoot_ratio <- log(dataF$rootshoot_ratio)+5
dataF$Saprotrophs2_100 <- dataF$Saprotrophs2/100
dataF$S_Sapro100 <- dataF$S_Sapro/100

#Sapro ab with shoot da negativo porq sapro increase con drought y shoot decreased with drought

######
way2<-'
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  PCo2_Patho+Saprotrophs2_100+PCo2_Mutual             #

PCo2_Patho ~  sqrt_SRSA +rootshoot_ratio100 
Saprotrophs2_100 ~ log_RTD
PCo2_Mutual ~ rootshoot_ratio100

PCo2_Patho ~~ Saprotrophs2_100 + PCo2_Mutual
Saprotrophs2_100 ~~ PCo2_Mutual
'
fit.way2 <- sem(way2, data=dataF)
varTable(fit.way2) # ajustarlo dividiendo entre 100 la variable de alta variance

summary(fit.way2, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.way2, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # 0.08

#############################################

###### PEARSON  correlations Root traits and Fungal attributes


# correr atras en el path las transformaciones de estas variables

root<- dataF[,c("log_rootC", "log_rootN","sqrt_root","sqrt_diameter","log_RTD","log_SRL","sqrt_SRSA")]
root<- data.matrix(root)
fungal<- dataF[,c("Pathogens2","S_Patho","PCo1_Patho","PCo2_Patho",
                  "Saprotrophs2", "S_Sapro","PCo1_Sapro", "PCo2_Sapro",
                  "Mutualists2","S_Mutual","PCo1_Mutual","PCo2_Mutual")]
fungal<- data.matrix(fungal)
rootfungal <-cor(root,fungal)
corrplot(rootfungal) 

str(dataF)


log_fungal<- dataF[,c("log_Pathogen","log_S_Patho","log_PCo1_Patho","log_PCo2_Patho",
                  "log_Sapro", "log_S_Sapro","log_PCo1_Sapro", "log_PCo2_Sapro",
                  "log_Mutual","log_S_Mutual","log_PCo1_Mutual","log_PCo2_Mutual")]
log_fungal<- data.matrix(log_fungal)
rootlog_fungal <-cor(root,log_fungal)
corrplot(rootlog_fungal) 

##########################################################################################

#### RELAIMPO USING PMVD ########### 

# RELAIMPO _Groemping  ONLY ROOT TRAITS (abajo con leaf, root and soil properties)

# Utilizo las transformadas porq estas fueron las q use en el path analyses. Si no, salen cosas diferentes.  

#log shoot mass' mutual (root:shoot)

shootmass_mutual <-lm(log_shootmass ~ Mutualists2 +S_Mutual + PCo1_Mutual  + PCo2_Mutual, data = dataF)
metrics_shootmass_mutual <- calc.relimp(shootmass_mutual, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_shootmass_mutual@pmvd = metrics_shootmass_mutual@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_shootmass_mutual=data.frame(variables=names(metrics_shootmass_mutual@pmvd),varexpl.percent=metrics_shootmass_mutual@pmvd) # nueva tabla

Tab_metrics_shootmass_mutual$variables <- factor(Tab_metrics_shootmass_mutual$variables,
                                                 levels=c("Mutualists2","PCo1_Mutual","S_Mutual","PCo2_Mutual"))  # order panels

REL_shootmass_mutual<- ggplot(data=Tab_metrics_shootmass_mutual, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Shoot mass _mutual', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_shootmass_mutual

#########

#### # Root_Mutual_AB 

Root_MutualAB <-lm(Mutualists2 ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)

metrics_rootMutAB <- calc.relimp(Root_MutualAB, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMutAB@pmvd = metrics_rootMutAB@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMutAB=data.frame(variables=names(metrics_rootMutAB@pmvd),varexpl.percent=metrics_rootMutAB@pmvd) # nueva tabla

Tab_metrics_rootMutAB$variables <- factor(Tab_metrics_rootMutAB$variables,
                                          levels=c("log_SRL","sqrt_diameter","rootshoot_ratio","log_rootN","log_rootC","sqrt_SRSA","log_RTD"))  # order panels
REL_MUT_rootMutAB<-ggplot(data=Tab_metrics_rootMutAB, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Mutualists abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMutAB


#### # Root_S_Mutual 

Root_S_Mutual.RIRS <-lm(log_S_Mutual ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)

metrics_rootMSRS <- calc.relimp(Root_S_Mutual.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMSRS@pmvd = metrics_rootMSRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMSRS=data.frame(variables=names(metrics_rootMSRS@pmvd),varexpl.percent=metrics_rootMSRS@pmvd) # nueva tabla

Tab_metrics_rootMSRS$variables <- factor(Tab_metrics_rootMSRS$variables,
                                         levels=c("log_RTD","log_SRL","log_rootN","sqrt_SRSA","log_rootC","sqrt_diameter","rootshoot_ratio"))  # order panels
REL_MUT_rootMS<-ggplot(data=Tab_metrics_rootMSRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Mutualists richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMS

########

## Root_ PCO2_Mutual 
Root_PCo2_Mutual.RIRS2 <-lm(PCo2_Mutual ~ log_rootN + log_rootC+ rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootMCRS2 <- calc.relimp(Root_PCo2_Mutual.RIRS2, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootMCRS2@pmvd = metrics_rootMCRS2@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootMCRS2=data.frame(variables=names(metrics_rootMCRS2@pmvd),varexpl.percent=metrics_rootMCRS2@pmvd) # nueva tabla

Tab_metrics_rootMCRS2$variables <- factor(Tab_metrics_rootMCRS2$variables,
                                          levels=c("log_rootC","log_SRL","sqrt_SRSA","sqrt_diameter","log_rootN","log_RTD","rootshoot_ratio"))  # order panels

REL_MUT_rootMC2<-ggplot(data=Tab_metrics_rootMCRS2, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Mutualists composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_rootMC2


############## saprotrophs

#log shoot mass_sapro 

shootmass_sapro <-lm(log_shootmass ~ Saprotrophs2+S_Sapro + PCo1_Sapro  + PCo2_Sapro, data = dataF)
metrics_shootmass_sapro <- calc.relimp(shootmass_sapro, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_shootmass_sapro@pmvd = metrics_shootmass_sapro@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_shootmass_sapro=data.frame(variables=names(metrics_shootmass_sapro@pmvd),varexpl.percent=metrics_shootmass_sapro@pmvd) # nueva tabla

#Tab_metrics_shootmass_sapro$variables <- factor(Tab_metrics_shootmass_sapro$variables,
#                                        levels=c("log_rootN","log_SRL","rootshoot_ratio","log_rootC","sqrt_diameter","sqrt_SRSA","log_RTD"))  # order panels

REL_shootmass_sapro<- ggplot(data=Tab_metrics_shootmass_sapro, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Shoot mass _sapro', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_shootmass_sapro


############

#Root_ SAPRO_AB   

Root_Saprotrophs.RIRS <-lm(Saprotrophs2 ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootSabRS <- calc.relimp(Root_Saprotrophs.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSabRS@pmvd = metrics_rootSabRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSabRS=data.frame(variables=names(metrics_rootSabRS@pmvd),varexpl.percent=metrics_rootSabRS@pmvd) # nueva tabla

Tab_metrics_rootSabRS$variables <- factor(Tab_metrics_rootSabRS$variables,
                                          levels=c("log_SRL","log_rootN","rootshoot_ratio","log_rootC","sqrt_diameter","sqrt_SRSA","log_RTD"))  # order panels

REL_SAP_rootAB<- ggplot(data=Tab_metrics_rootSabRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Saprotrophs abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootAB

################
#Root_SAPRO_S  
Root_S_Saprotrophs.RIRS <-lm(S_Sapro ~log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA, data = dataF)
metrics_rootSSRS <- calc.relimp(Root_S_Saprotrophs.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSSRS@pmvd = metrics_rootSSRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSSRS=data.frame(variables=names(metrics_rootSSRS@pmvd),varexpl.percent=metrics_rootSSRS@pmvd) # nueva tabla

Tab_metrics_rootSSRS$variables <- factor(Tab_metrics_rootSSRS$variables,
                                         levels=c("rootshoot_ratio","log_rootN","sqrt_diameter","log_rootC","sqrt_SRSA","log_SRL","log_RTD"))  # order panels

REL_SAP_rootSS<- ggplot(data=Tab_metrics_rootSSRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Saprotrophs richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootSS


##Root SAPRO_PCOA2    
names(dataF)
Root_PCo2_Saprotrophs.RIRS2 <-lm(PCo2_Sapro ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootSCRS2 <- calc.relimp(Root_PCo2_Saprotrophs.RIRS2, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootSCRS2@pmvd = metrics_rootSCRS2@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootSCRS2=data.frame(variables=names(metrics_rootSCRS2@pmvd),varexpl.percent=metrics_rootSCRS2@pmvd) # nueva tabla

Tab_metrics_rootSCRS2$variables <- factor(Tab_metrics_rootSCRS2$variables,
                                          levels=c("log_rootC","log_rootN", "rootshoot_ratio","sqrt_diameter","log_RTD","sqrt_SRSA","log_SRL")) # order panels

REL_SAP_rootSC2<- ggplot(data=Tab_metrics_rootSCRS2, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Saprotrophs composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_rootSC2

####PAthogen

#log shoot mass_patho (root:shoot)

shootmass_patho <-lm(log_shootmass ~ Pathogens2+ S_Patho +  PCo1_Patho  + PCo2_Patho, data = dataF)
metrics_shootmass_patho <- calc.relimp(shootmass_patho, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_shootmass_patho@pmvd = metrics_shootmass_patho@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_shootmass_patho=data.frame(variables=names(metrics_shootmass_patho@pmvd),varexpl.percent=metrics_shootmass_patho@pmvd) # nueva tabla

Tab_metrics_shootmass_patho$variables <- factor(Tab_metrics_shootmass_patho$variables,
                                                levels=c("Pathogens2", "S_Patho", "PCo1_Patho", "PCo2_Patho"))  # order panels

REL_shootmass_patho<- ggplot(data=Tab_metrics_shootmass_patho, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Shoot mass _patho', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_shootmass_patho


#Root PATHO_AB

Root_Pathogens.RIRS <-lm(Pathogens2 ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPabRS <- calc.relimp(Root_Pathogens.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPabRS@pmvd = metrics_rootPabRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPabRS=data.frame(variables=names(metrics_rootPabRS@pmvd),varexpl.percent=metrics_rootPabRS@pmvd) # nueva tabla

Tab_metrics_rootPabRS$variables <- factor(Tab_metrics_rootPabRS$variables,
                                          levels=c("log_rootC","log_rootN","log_RTD","sqrt_SRSA","log_SRL","rootshoot_ratio","sqrt_diameter"))  # order panels

REL_rootPAT_AB<- ggplot(data=Tab_metrics_rootPabRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Pathogens abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_AB


#RELAIMPO PATHO_S 

Root_S_Pathogens.RIRS <-lm(S_Patho ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPSRS <- calc.relimp(Root_S_Pathogens.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPSRS@pmvd <- metrics_rootPSRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPSRS <- data.frame(variables=names(metrics_rootPSRS@pmvd),varexpl.percent=metrics_rootPSRS@pmvd) # nueva tabla

Tab_metrics_rootPSRS$variables <- factor(Tab_metrics_rootPSRS$variables,
                                         levels=c("log_rootN","rootshoot_ratio","log_rootC","sqrt_SRSA","log_RTD","sqrt_diameter","log_SRL"))  # order panels

REL_rootPAT_PS<- ggplot(data=Tab_metrics_rootPSRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Pathogens richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_PS


#RELAIMPO PATHO_PCOA2 #### Including Root:shoot 
names(dataF)
Root_PCo2_Pathogens.RIRS <-lm(PCo2_Patho ~ log_rootN+log_rootC+rootshoot_ratio+sqrt_diameter+log_RTD+log_SRL+sqrt_SRSA,data = dataF)
metrics_rootPCRS <- calc.relimp(Root_PCo2_Pathogens.RIRS, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_rootPCRS@pmvd = metrics_rootPCRS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_rootPCRS=data.frame(variables=names(metrics_rootPCRS@pmvd),varexpl.percent=metrics_rootPCRS@pmvd) # nueva tabla

Tab_metrics_rootPCRS$variables <- factor(Tab_metrics_rootPCRS$variables,
                                         levels=c("log_rootN","sqrt_diameter","log_rootC","log_RTD","log_SRL","rootshoot_ratio","sqrt_SRSA"))  # order panels

REL_rootPAT_PC2<- ggplot(data=Tab_metrics_rootPCRS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Pathogens composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_rootPAT_PC2


## 
ggarrange (REL_rootPAT_AB, REL_rootPAT_PS, REL_rootPAT_PC2,
           REL_SAP_rootAB,REL_SAP_rootSS,REL_SAP_rootSC2, 
           REL_MUT_rootMutAB, REL_MUT_rootMS, REL_MUT_rootMC2,
           labels = c("A", "B","C","D","E","F","G", "H", "I"),
           common.legend = FALSE, align = c("hv"))

ggarrange (REL_shootmass_patho, REL_shootmass_sapro, REL_shootmass_mutual,
           REL_shootmass_patho, REL_shootmass_sapro, REL_shootmass_mutual,
           REL_shootmass_patho, REL_shootmass_sapro, REL_shootmass_mutual,
           labels = c("A", "B","C"),
           common.legend = FALSE, align = c("hv"), ncol =3, nrow = 3)

# RAREFACTION 
library(vegan)
dataESV <- read.csv("ESV_totalFungComposition_mantel.csv")  # la tabla donde esta la "total" fungal composition. The ESV selected
esv_total <-dataESV %>% dplyr::select(starts_with("esv"))

S_esv<-specnumber(esv_total)
(raremax <- min(rowSums(esv_total))) #rows are the samples ( here each plant species). Se suman each y cual es el valor minimo encontrado en una row
S_rare <- rarefy(esv_total, raremax)
plot(S_esv, S_rare, xlab = "Observed No. of ESV", ylab = "Rarefied No. of ESV")
abline(0, 1)
rarecurve(esv_total, step = 20, col = "red", cex = 0.6)#, sample = raremax)

#rarefacton solo con "esv1"
esv_total2 <-NP_total %>% dplyr::select(starts_with("esv1"))
S_esv2<-specnumber(esv_total2)
(raremax2 <- min(rowSums(esv_total2))) #rows are the samples ( here each plant species). Se suman each y cual es el valor minimo encontrado en una row
S_rare2 <- rarefy(esv_total2, raremax2)
plot(S_esv2, S_rare2, xlab = "Observed No. of ESV", ylab = "Rarefied No. of ESV")
abline(0, 1)
rarecurve(esv_total2, step = 20, col = "red", cex = 0.6)#, sample = raremax)

###################################################

###################################################

####### NO INCLUIDO EN EL PAPER
##############################################################
##### PATH ANALYSES (2)
#... aqui he hecho.. Drought->fungal->shoot/diameter


# 1. teniendo shoot mass and root diameter as the variables
# 2. Teneindo "leaf traits" and "root traits" as the variables??? 
names(dataF)


model1 <- ' #incluyendo soil properties               ### LAS FUNGAL ATRIBUTES IMPORTANTES PARA SHOOT MASS
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho+
                  Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro+
                  Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual
diameter ~~ Pathogens2+log_Pathogen+ S_Patho+ log_S_Patho+ PCo1_Patho+ PCo2_Patho+
                  Saprotrophs2+log_Sapro+ S_Sapro+ log_S_Sapro+ PCo1_Sapro+ PCo2_Sapro+
                  Mutualists2+log_Mutual+ S_Mutual+ log_S_Mutual+ PCo1_Mutual+ PCo2_Mutual   ## no sabemos quien afecta a quien
Pathogens2 ~ TtoValue
log_Pathogen ~TtoValue
S_Patho ~TtoValue
log_S_Patho ~TtoValue
PCo1_Patho ~TtoValue
PCo2_Patho ~TtoValue
Saprotrophs2 ~TtoValue
log_Sapro ~TtoValue
S_Sapro ~TtoValue
log_S_Sapro ~TtoValue
PCo1_Sapro ~TtoValue
PCo2_Sapro ~TtoValue
Mutualists2 ~TtoValue
log_Mutual ~TtoValue
S_Mutual ~TtoValue
log_S_Mutual ~TtoValue
PCo1_Mutual ~TtoValue
PCo2_Mutual ~TtoValue

diameter ~TtoValue
log_shootmass ~ TtoValue+diameter
CP2_soilprop ~TtoValue

#remove close to cero to improve model

'
fit.mod1 <- sem(model1, data=dataF)
varTable(fit.mod1)
#dataF <- dataF %>% mutate(Saprotrophs2 = Saprotrophs2/100) # para ajustarlo a los demas # si lo corrigo o no el fitmeasures da lo mismo 

summary(fit.mod1, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.mod1, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # rmsea 0.15 
######

modelA <- ' #incluyendo soil properties               ### LAS FUNGAL ATRIBUTES IMPORTANTES PARA SHOOT MASS
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  PCo2_Patho + PCo2_Sapro + log_S_Mutual 
diameter ~~ PCo2_Patho + PCo2_Sapro     ## no sabemos quien afecta a quien
PCo2_Patho ~ TtoValue
PCo2_Sapro ~TtoValue
log_S_Mutual ~TtoValue
diameter ~TtoValue
log_shootmass ~ TtoValue+diameter
CP2_soilprop ~TtoValue

#remove close to cero to improve model
log_shootmass ~~ 0* CP2_soilprop
log_S_Mutual ~~ 0* diameter
'
fit.modA <- sem(modelA, data=dataF)
varTable(fit.modA)
#dataF <- dataF %>% mutate(Saprotrophs2 = Saprotrophs2/100) # para ajustarlo a los demas # si lo corrigo o no el fitmeasures da lo mismo 

summary(fit.modA, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.modA, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # rmsea 0.15 

#### 
names(dataF)
modelB <- ' #incluyendo soil properties                ### LOS FUNGAL ATTRIBUTES IMPORTANTES PARA DIAMETER
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  + S_Patho + S_Sapro + S_Mutual 
diameter ~~ S_Patho + S_Sapro + S_Mutual      ## no sabemos quien afecta a quien
S_Patho ~ TtoValue
S_Sapro ~TtoValue
S_Mutual ~TtoValue
diameter ~TtoValue
log_shootmass ~ TtoValue+diameter
CP2_soilprop ~ TtoValue

#remove close to cero to improve model
log_shootmass ~ 0*CP2_soilprop
'
fit.modB <- sem(modelB, data=dataF)
varTable(fit.modB)
#dataF <- dataF %>% mutate(Saprotrophs2 = Saprotrophs2/100) # para ajustarlo a los demas # si lo corrigo o no el fitmeasures da lo mismo 

summary(fit.modB, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.modB, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # rmsea 0.4

######
modelC <- ' #incluyendo soil properties                ### LOS FUNGAL ATTRIBUTES IMPORTANTES PARA leaf TRAITS
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  + PCo2_Patho + S_Sapro + log_Mutual 
diameter ~~ + S_Sapro + log_Mutual      ## no sabemos quien afecta a quien
PCo2_Patho ~ TtoValue
S_Sapro ~TtoValue
log_Mutual ~TtoValue
diameter ~TtoValue
log_shootmass ~ TtoValue+diameter
CP2_soilprop ~ TtoValue

#remove close to cero to improve model
log_shootmass ~~ 0* CP2_soilprop
PCo2_Patho ~~ 0* diameter
'
fit.modC <- sem(modelC, data=dataF)
varTable(fit.modC)
#dataF <- dataF %>% mutate(Saprotrophs2 = Saprotrophs2/100) # para ajustarlo a los demas # si lo corrigo o no el fitmeasures da lo mismo 

summary(fit.modC, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.modC, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # rmsea 0.19
######

modelD <- ' #incluyendo soil properties                ### LOS FUNGAL ATTRIBUTES IMPORTANTES PARA ROOTS TRAITS
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  + log_Pathogen + Saprotrophs2 + S_Mutual 
diameter ~~ + log_Pathogen + Saprotrophs2 + S_Mutual      ## no sabemos quien afecta a quien
log_Pathogen ~ TtoValue
Saprotrophs2 ~TtoValue
S_Mutual ~TtoValue
diameter ~TtoValue
log_shootmass ~ TtoValue+diameter
CP2_soilprop ~ TtoValue

#remove close to cero to improve model

'
fit.modD <- sem(modelD, data=dataF)
varTable(fit.modD)
#dataF <- dataF %>% mutate(Saprotrophs2 = Saprotrophs2/100) # para ajustarlo a los demas # si lo corrigo o no el fitmeasures da lo mismo 

summary(fit.modD, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.modD, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # rmsea 0.19

###
modelE <- ' #incluyendo soil properties                ### LOS  ATTRIBUTES de todo IMPORTANTES PARA DROUGHT
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  + S_Patho +  log_Mutual 
diameter ~~ + PCo1_Patho +  log_Mutual     ## no sabemos quien afecta a quien
S_Patho ~ TtoValue
log_Mutual ~TtoValue
diameter ~TtoValue
log_shootmass ~ TtoValue+diameter
CP2_soilprop ~ TtoValue

#remove close to cero to improve model

'
fit.modE <- sem(modelE, data=dataF)
varTable(fit.modE)
#dataF <- dataF %>% mutate(Saprotrophs2 = Saprotrophs2/100) # para ajustarlo a los demas # si lo corrigo o no el fitmeasures da lo mismo 

summary(fit.modE, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.modE, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # rmsea 0.27

# si quito log_Sapro (mejor muy poco el RMSEA 0.21)

##############################

# En verdad DROUGHT no affecta mucho..es mas la especie... entonces.. ver como los PLANT TRAITS  affecta las microbial communities
# luego se puede discutir q estos plant traits estan also afectado por drought...

cor(dataF[12:43])
`names(dataF)
names(dataF)
mydata <- cbind(dataF$sqrt_diameter, dataF$sqrt_root, dataF$sqrt_SRSA, dataF$log_shootmass,
    dataF$Pathogens2, dataF$log_Pathogen, dataF$S_Patho, dataF$log_S_Patho,
    dataF$Saprotrophs2, dataF$log_Sapro, dataF$S_Sapro, dataF$log_S_Sapro,
    dataF$Mutualists2, dataF$log_Mutual, dataF$S_Mutual, dataF$log_S_Mutual)#, method = c('pearson'))

as.matrix(mydata)

mydata1<-c(dataF$sqrt_diameter, dataF$sqrt_root, dataF$sqrt_SRSA, dataF$log_shootmass,
      dataF$Pathogens2, dataF$log_Pathogen, dataF$S_Patho, dataF$log_S_Patho,
      dataF$Saprotrophs2, dataF$log_Sapro, dataF$S_Sapro, dataF$log_S_Sapro,
      dataF$Mutualists2, dataF$log_Mutual, dataF$S_Mutual, dataF$log_S_Mutual)
mydata2 <- c("sqrt_diameter", "sqrt_root", "sqrt_SRSA", data = dataF)

head(mydata2)
cor(mydata2, use = "complete.obs")

cor(dataF$sqrt_diameter, dataF$sqrt_root, dataF$sqrt_SRSA, method = "pearson")#dataF$sqrt_SRSA, #dataF$log_shootmass,
         dataF$Pathogens2, dataF$log_Pathogen, dataF$S_Patho, dataF$log_S_Patho,
         dataF$Saprotrophs2, dataF$log_Sapro, dataF$S_Sapro, dataF$log_S_Sapro,
         dataF$Mutualists2, dataF$log_Mutual, dataF$S_Mutual, dataF$log_S_Mutual, method = "pearson")

install.packages("Hmisc")
library(Hmisc)
library(corrplot)


#####################################################################################################



############################################################################################################
############################################################################################################

## OLD ANALYSES 

###########################################
### CORRELATIONS TRAITS -FUNGAL

dcontrol = subset(data, Tto == "control")
ddrought = subset(data, Tto == "drought")

#datacontrol <- data[ which(data$tto=='control'), ]  subset is better
#datadrought <- data[ which(data$tto=='drought'), ]

str(dcontrol)
names(data)

#control 

#AMF = sqrt_AMF   VS sqrt_rootmass, sqrt_SRSA, log_SRL  # mejor root y SRSA
      # sqrt_Samf VS sqrt_rootmass, sqrt_SRSA
      #AMF_comp_2 VS log_SRL                         ## no good
#drought = no esta correlacionado  con ningun root trait 

data_herblegcontrol = subset(dcontrol, PlantGroup == "Herb" | PlantGroup == "Legume")
data_grasses = subset(data, PlantGroup == "Grass")
View(data_herbleg)

amf1<-ggplot(data, aes(x=sqrt_AMF, y=sqrt_RDA))+
  geom_point(colour="black",size=3)+
  geom_smooth(method=lm, se=TRUE)+ # ahora decoration
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"))+
  theme( axis.line = element_line(size = 0.5, linetype = "solid"))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"))+
  labs(title="",x ="sqrt (AMF abundance)", y = "sqrt (Root diameter)")+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black") # Add axis line
  )
amf1

#graph Herbs and legumes correlated with diameter 
data_herbleg = subset(datacontrol, PlantGroup == "Herb" | PlantGroup == "Legume")
data_grasses = subset(data, PlantGroup == "Grass")
View(data_herbleg)

amf2<-ggplot(data, aes(x=AMF_comp_2, y=sqrt_RDA))+
  geom_point(colour="black",size=3)+
  geom_smooth(method=lm, se=TRUE)+ # ahora decoration
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"))+
  theme( axis.line = element_line(size = 0.5, linetype = "solid"))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"))+
  labs(title="",x ="AMF composition", y = "sqrt (Root diameter)")+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) # Add axis line
amf2

#pathogens
names(data)
data_herbleg = subset(data, PlantGroup == "Herb" | PlantGroup == "Legume")
data_grassdr = subset(ddrought, PlantGroup == "Grass")

# pathogen composition was linked to root diameter and soil properties in grasses (drought)
#saprotrophs negative con NO3 under control. 
names
s1<-ggplot(data, aes(x=sqrt_sapro, y=NO3.mg.l.))+
  geom_point(colour="black",size=3)+
  geom_smooth(method=lm, se=TRUE)+ # ahora decoration
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme( axis.line = element_line(size = 0.5, linetype = "solid"))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))+
  labs(title="",x ="sqrt (Saprotroph abundance)")+
    theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))+
    ylab(bquote(bold('NO3 ('*~mg~kg^-1*')')))# Add axis line
s1

ggarrange(amf1, s1,labels = c("A", "B"),
          common.legend = TRUE, align=c("hv")) # allign del mismo size horiz y vertic

#regresion for the R2 value
summary(lm(NO3.mg.l.~sqrt_sapro, data=data))
summary(lm(sqrt_SRSA~sqrt_AMF, data=dcontrol))
summary(lm(sqrt_rootmass~sqrt_AMF, data=dcontrol))
summary(lm(sqrt_RDA~sqrt_AMF, data=data))
summary(lm(sqrt_RDA~AMF_comp_2, data=data))
summary(lm(sqrt_RDA~Pathog_comp_1, data=data))

#PEARSON Correlation 

#regression for drought and for non-drought
datacontrol <- data[ which(data$tto=='control'), ]
datadrought <- data[ which(data$tto=='drought'), ]


#log_pathogen5, sqrt_sapro, sqrt_AMF (0.79), sqrt_generpathog,sqrt_spepathog,sqrt_gensapro
#             specialistsaprotro, sqrt_speAMF
#richness: spathogen, sqrt_S_sapro, sqrt_S_AMF (0.84), generalist_pathog_S, specialist_pathog_S
#         generalist_sapro_S, sqrt_spesapro_S, sqrt_specialist_AMF_S
#shoot mass. it is better log shoot mass

data$log_pathogen5 <- as.numeric(data$log_pathogen5)
cor(datacontrol$log_pathogen5,datacontrol$log_shootmass, method = c('pearson'))
cor(datacontrol$sqrt_RDA, datacontrol$specialistsaprotro, method = c('pearson'))
str(data)


#playing
datap= read.csv("trialdata.csv")
datacontrol <- datap[ which(datap$tto=='control'), ]
datadrought <- datap[ which(datap$tto=='drought'), ]
cor(datacontrol$pathogen,datacontrol$shootmass, method = c('pearson'))
cor.test(datacontrol$pathogen,datacontrol$shootmass, method = 'pearson')


###REGRESIONES USING RII  (ALL data)

#1. RII_pathog_11 ~ RII_SoilNgkg
#2. RII_saprotr_11 ~ Tot_RootBio_RII

dataRII <- read.csv("RII_fungal_root_traits.csv") 

library(ggplot2)

modelRII=lm(RII_saprotr_11 ~ Tot_RootBio_RII, data = dataRII)
summary(modelRII)

regreRII<-ggplot(dataRII, aes(x=RII_saprotr_11, y=Tot_RootBio_RII))+
  geom_point(colour="black",size=3)+
  geom_smooth(method=lm, se=TRUE)+ # ahora decoration
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme( axis.line = element_line(size = 0.5, linetype = "solid"))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black") # Add axis line
  )
regreRII

###REGRESIONES USING RII  (MEAN por SPECIES)

# lo he intentado con sapro y root biomass...pero el result es peor q con todos los 240 puntos :(

dataRIImean <- read.csv("RII_fungal_root_mean_species.csv") 

library(ggplot2)

modelRIImean=lm(RII_saprotr_11 ~ Tot_RootBio_RII, data = dataRIImean)
summary(modelRIImean)

regreRIImean<-ggplot(dataRIImean, aes(x=RII_saprotr_11, y=Tot_RootBio_RII))+
  geom_point(colour="black",size=3)+
  geom_smooth(method=lm, se=TRUE)+ # ahora decoration
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme( axis.line = element_line(size = 0.5, linetype = "solid"))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black") # Add axis line
  )
regreRIImean

#REGRESIONES y GRAFICO  con RII

data= read.csv("regresion_pearson_RII y RIIvsControl.csv") 

#Pearson correlation que dieron significativas 

#AvDiam_RII       	RII_saprotr_av1  	118	  -0.18	 0.0509
#AvDiam_RII       	RII_AMF_av       	 95	  -0.21	 0.0429
#Shoot_RII        	RII_AMF_11       	117	   -0.29	 0.0013
#RTD_RII          	RII_S_pathogen_11	118	   -0.16	 0.0895
#SRL_RII          	RII_AMF_av       	 95	   0.21	 0.0373
#SRSA_RII         	RII_AMF_av       	 95	    0.17	 0.0910
#Tot_RootBio_RII  	RII_saprotr_11   	117	   0.18	 0.0567
#RII_SoilNgkg     	RII_pathog_11    	116	   0.17	 0.0652
#RII_SoilNgkg      	RII_S_AMF_11     	 62	   0.26	 0.0382
#RII_SoilCgkg      	RII_S_AMF_11     	 62	   0.26	 0.0414

##### hacer tb RII fungal vs control Root traits 

					
#Sqrt_Avdiam_control    RII_saprotr_av1 
#Sqrt_Avdiam_control
#log_shootmass_control
#logRTD_control
#logSRL3_control
#sqrt_SRSA_control
#sqrt_TotRootB_control
#	log_rootN_control	Roots_C_gkg_control


r<-ggplot(data, aes(x=RII_SoilCgkg, y=RII_S_AMF_11))+
  geom_point(colour="black",size=3)+
  geom_smooth(method=lm, se=TRUE)+ # ahora decoration
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme( axis.line = element_line(size = 0.5, linetype = "solid"))+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black") # Add axis line
  )
r

##############################################################################
#############################   MANTEL TEST ############
#########################################
## PARA SABER SI HAY CORRELACION ENTRE ESTAS DOS MATRICES!!! 

dataESV <- read.csv("ESV_totalFungComposition_mantel.csv")
totESV<-dataESV[, grep(pattern="^esv", colnames(dataESV))]  # seleccionar las columnas que empiezan con "esv" 

dataF <- read.csv ("traits_fungal.csv") # aqui estan los root traits tb

names(dataF)
head(data)
# Composition (ESV) vs root traits. 
#dos formas. Usando vegdist ("jaccard) Marina or dist.binary (mejor results with the first)
#OJO HE QUITADO LA ### 232  #### que no habia datos en ROOT TRAITS
###OJO. La primera columna siempre debe llevar el ID. La celda 1A vacia

#root trait matriz

dataF$sqrt_diameter <- sqrt(dataF$diameter)
dataF$log_SRL <- log(dataF$SRL)
dataF$sqrt_SRSA <- sqrt(dataF$SRSA)
dataF$sqrt_rootmass <-sqrt(dataF$rootmass)
dataF$log_rootN <- log(dataF$rootN)

roottrait <-(cbind(dataF$sqrt_diameter, dataF$log_SRL, dataF$sqrt_SRSA, dataF$sqrt_rootmass, dataF$log_rootN, dataF$rootC)) #matriz root traits#
roottrait.dist <- vegdist(scale(roottrait), "euclid")

rawroottrait = cbind(dataF$diameter, dataF$SRL, dataF$SRSA, dataF$rootmass, dataF$rootN, dataF$rootC) #raw data
rawroottrait.dist = vegdist(scale(rawroottrait), "euclid")

# Total fungal composition_VEGDIST 
ESV.dist <- vegdist(totESV, method="bray") # Sin poner Jaccard es Bry Curtis

mantel(ESV.dist, roottrait.dist) #pearson
mantel(ESV.dist, roottrait.dist, method="spearman")  
mantel(ESV.dist, rawroottrait.dist, method="spearman")  
mantel(ESV.dist, rawroottrait.dist) #pearson

#Total Fungal compo_ JACCARD. dIST.BINARY  
ESVbin = binarize(totESV, threshold = 0)   #que sea binaria: presence-absence. 
ESV.jac.dist = dist.binary(ESVbin,  method = 1) #1=jaccard

mantel(ESV.jac.dist, roottrait.dist) #root traits raw here. Transform is pract the same
mantel(ESV.jac.dist, rawroottrait.dist)
mantel(ESV.jac.dist, roottrait.dist, method = "spear")
mantel(ESV.jac.dist, rawroottrait.dist, method = "spear") 


#Pathogen composition. "BRAY" Tengo plant species con CERO patho. quitar estos pot en root traits. 

#datap= read.csv("Roottraits_4_pathogenESV.csv")
#roottrait.pat=dist(cbind(datap$RDA, datap$SRL, datap$SRSA, datap$rootmass, datap$rootN_g, datap$rootC_g))
#roottrait.pat.dist <- vegdist(scale(roottrait.pat), "euclid")
#pathoESV= read.csv("pathogen_ESV_NOALL.csv")  # He quitado las q era ceros# #dejar solo las var a usar. No pot or species column.
#pathoESV.dist <- vegdist(pathoESV)  
#mantel(pathoESV.dist, roottrait.pat.dist) 
#mantel(pathoESV.dist, roottrait.pat.dist, method="spear")

#Patho composition "VEGDIST JACARD or bray"   # con ninguno da significativo
patho_allESV = read.csv("pathogen_ESV_ALL.csv")
pathoESV.dist <- vegdist(patho_allESV, method="bray")  # puedo usar Jacard ver Semchenko 2018
mantel(pathoESV.dist, roottrait.dist)                  # n.s
mantel(pathoESV.dist, roottrait.dist, method="spear")   #n.s

mantel(pathoESV.dist, rawroottrait.dist)                #n.s
mantel(pathoESV.dist, roottrait.dist, method="spear")   #n.s

#Pathogen compo_ JACCARD DIST.BINARY 
pathoESVbin = binarize(patho_allESV, threshold = 0)   
pathoESV.jac.dist = dist.binary(pathoESVbin,  method = 1) 
mantel(pathoESV.jac.dist, roottrait.dist)                   # significativo!!! 
mantel(pathoESV.jac.dist, roottrait.dist, method = "spear")  # p=0.06


#Sapro composition "VEGDIST JACARD or Bray"
saproESV= read.csv("sapro_ESV.csv")   
saproESV.dist <- vegdist(saproESV, method="bray") 
mantel(saproESV.dist, roottrait.dist) 
mantel(saproESV.dist, roottrait.dist, method="spear") #p=0.06

mantel(saproESV.dist, rawroottrait.dist) 
mantel(saproESV.dist, roottrait.dist, method="spear") #p=0.07

#sapro compo "Jaccard. Dist binary 
saproESVbin = binarize(saproESV, threshold = 0)   
saproESV.jac.dist = dist.binary(saproESVbin,  method = 1)
mantel(saproESV.jac.dist, roottrait.dist) 
mantel(pathoESV.jac.dist, roottrait.dist, method = "spear")  #significativo

#  AMF_compo  Vegdist Jaccard or bray (no diferencias), 
AMFESV= read.csv("AMF_ESV.csv")   
AMFESV.dist <- vegdist(AMFESV, method="bray") 
mantel(AMFESV.dist, roottrait.dist) 
mantel(AMFESV.dist, roottrait.dist, method="spear")

mantel(AMFESV.dist, rawroottrait.dist) # significativo
mantel(AMFESV.dist, rawroottrait.dist, method="spear") #significativo


#AMF composition _vegdist Jaccard 
AMFESVbin = binarize(AMFESV, threshold = 0)   
AMF.jac.dist = dist.binary(AMFESVbin,  method = 1) # 1 means Jaccard
mantel(AMF.jac.dist, roottrait.dist)
mantel(AMF.jac.dist, roottrait.dist, method = "spear")


#USING BOTH JACCARD. Root trait dist ALL is 1. TODO 1 ...MMMMM. results weird.. NO tiene como sentido 
rootTbin = binarize(roottrait, threshold = 0)
roottrait.jac.dist = dist.binary(rootTbin, method = 1)
mantel(AMF.jac.dist, rootT.jac.dist)

####################################
##############################################################


################## PHYLA- PHYLUM abundances 

#p__Ascomycota   #p__Basidiomycota   #p__Cercozoa	
#p__Chytridiomycota	p__Entomophthoromycota (casi cero)	p__Mortierellomycota	
#p__Mucoromycota	p__Rozellomycota	p__Zoopagomycota	
#c__Archaeosporomycetes	c__Glomeromycetes	c__Paraglomeromycetes
str(dataF)

dataF$Tto <- factor(data$Tto, levels=c("control", "drought"))

anova(lm(p__Ascomycota~species*Tto, data=data))
ggplot(data, aes(x=species, y=p__Ascomycota)) + 
  geom_boxplot()+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))


anova(lm(p__Basidiomycota~species*Tto, data=data))
ggplot(data, aes(x=species, y=p__Basidiomycota, fill=Tto)) +
  geom_boxplot()+
  scale_fill_manual(values=c("darkturquoise","firebrick"))+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

sqrt_p__Chytridiomycota=sqrt((data$p__Chytridiomycota))
anova(lm(sqrt_p__Chytridiomycota~species*Tto, data=data))
ggplot(data, aes(x=species, y=p__Chytridiomycota)) + 
  geom_boxplot()+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

sqrt_p__Mortierellomycota=sqrt((data$p__Mortierellomycota))
anova(lm(sqrt_p__Mortierellomycota~species*Tto, data=data))
ggplot(data, aes(x=species, y=p__Mortierellomycota, fill=Tto)) +
  geom_boxplot()+
  scale_fill_manual(values=c("darkturquoise","firebrick"))+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

sqrt_p__Mucoromycota=sqrt((data$p__Mucoromycota))
anova(lm(sqrt_p__Mucoromycota~species*Tto, data=data))
ggplot(data, aes(x=species, y=p__Mucoromycota, fill=Tto)) +
  geom_boxplot()+
  scale_fill_manual(values=c("darkturquoise","firebrick"))+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

sqrt_p__Rozellomycota=sqrt((data$p__Rozellomycota))
anova(lm(sqrt_p__Rozellomycota~species*Tto, data=data))
ggplot(data, aes(x=species, y=p__Rozellomycota, fill=Tto)) +
  scale_fill_manual(values=c("darkturquoise","firebrick"))+
  geom_boxplot()+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

sqrt_c__Archaeosporomycetes=sqrt((data$c__Archaeosporomycetes))
anova(lm(sqrt_c__Archaeosporomycetes~species*Tto, data=data))
ggplot(data, aes(x=species, y=c__Archaeosporomycetes, fill=Tto)) +
  scale_fill_manual(values=c("darkturquoise","firebrick"))+
  geom_boxplot()+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

sqrt_c__Glomeromycetes=sqrt((data$c__Glomeromycetes))
anova(lm(sqrt_c__Glomeromycetes~species*Tto, data=data))

sqrt_c__Paraglomeromycetes=sqrt((data$c__Paraglomeromycetes))
anova(lm(sqrt_c__Paraglomeromycetes~species*Tto, data=data)) # Infostat signif interaccion.heteroc correg
ggplot(data, aes(x=species, y=c__Paraglomeromycetes)) +
  scale_fill_manual(values=c("darkturquoise","firebrick"))+
  geom_boxplot()+
  theme(axis.text.x = element_text(face="bold", color="black",size=10, angle=90))

mean(data$c__Glomeromycetes, data=data)
aggregate(data$p__Mucoromycota,by=list(data$species),mean,na.rm=TRUE)


### PATH ANALYSES  ( He usado el modelo 6 ed la segunda option)

dfcontrol= read.csv("rawdata_control_fungal.csv")
dfdrought=read.csv("rawdata_drought_fungal.csv")
str(dfcontrol)

##############
######
# Option 1. HAcer un stepAIC para esocger las variables y luego el path analyses. (como en soil blocks)

# 
# 1. Select the variables que van a ir en el modelo 

# a)Pathogens..He puesto todas las de pathogenos para escoger las mejores predictores.
M1.control.path <- lm(log_shootmass ~ log_pathogen5+RAIZ_generalistpathogen1                     # abundance
                      +RAIZ_specialistpathogen1+S_pathog+generalist_pathog_S+specialist_pathog_S # richness
                      ++Pathg_comp_1+ Pathg_comp_2, data=dfcontrol)                              # composition
fitM1path <- stepAIC(M1.control.path, direction = 'backward')
 #b) saprotrophs
M1.control.sapr <- lm(log_shootmass ~ RAIZ_sapro+RAIZ_generalistsaprotro1+specialistsaprotro 
                      +RAIZ_S_sapro+generalist_sapro_S+RAIZ_specialist_sapro_S+
                        Sapro_Comp_1+Sapro_comp2, data=dfcontrol)
fitM1sapr <- stepAIC(M1.control.sapr, direction = 'backward')               

#c) AMF
M1.control.AMF <- lm(log_shootmass ~ AMF_comp_1+AMF_comp_2+RAIZ_AMF+RAIZ_specialistAMF+RAIZ_S_AMF
                     +RAIZ_specialist_AMF_S,data=dfcontrol)
fitM1AMF <- stepAIC(M1.control.AMF, direction = 'backward')           

#d) soil properties
M1.control.soil <- lm(log_shootmass ~ sqrtSoilN+logSoilC5+RAIZ_NO3.mg.l.+RAIZ_PO4.mg.l.
                      +CP1_soilqual2+CP2_soilqual2,data=dfcontrol)
fitM1soil <- stepAIC(M1.control.soil, direction = 'backward')   

#e) Root traits
M1.control.root <- lm(log_shootmass ~ log_rootN + log_RootC+ sqrt_TotRootB
                      +sqrt_avdiam+logRTD+logSRl.3+sqrt_SRSA,data=dfcontrol)
fitM1root <- stepAIC(M1.control.root, direction = 'backward')   

str(dfcontrol)

# 2) Path analyses  (ME HE QUEDADO CON EL MODELO 4step)

# a) este modelo 6 fue el mejor de abajo (usando como predictores el AIC despues de linear regresion
#cada variables against shoot mass.... 
# he comprobado con el stepAIC selection y cambian algunos...ver el nuevo Path analyses 

m1step <- ' #
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  log_pathogen5+RAIZ_generalistsaprotro1+RAIZ_AMF+CP2_soilqual2+logSRl.3
log_pathogen5 ~ CP2_soilqual2+ logSRl.3
RAIZ_generalistsaprotro1 ~ CP2_soilqual2+ logSRl.3
RAIZ_AMF ~ CP2_soilqual2+ logSRl.3
'
fit.f1step <- sem(m1step, data=dfcontrol)
summary(fit.f1step, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f1step, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC"))  # RMESA = 0.28
 
    # sin correlaciones entre fungal groups. RMSEA = 0.28. el resto ok
    # correlacion entre AMF y pathogen . RMESA = 0.31..idem
    # correlacion entre AMF y sapro . RMESA= 0.16..idem


# Usar RTD y sin corrleciones entre fungal groups

m2step <- ' #
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  log_pathogen5+RAIZ_generalistsaprotro1+RAIZ_AMF+CP2_soilqual2+logRTD
log_pathogen5 ~ CP2_soilqual2+ logRTD
RAIZ_generalistsaprotro1 ~ CP2_soilqual2+ logRTD
RAIZ_AMF ~ CP2_soilqual2+ logRTD
'
fit.f2step <- sem(m2step, data=dfcontrol)
summary(fit.f2step, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f2step, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) 

        # sin correlaciones entre fungal groups. RMSEA = 0.3. el resto ok
        # correlacion entre AMF y pathogen . RMESA = 0.34..idem
        # correlacion entre AMF y sapro . RMESA= 0.16..idem

# usar diameter y sin correlaciones entre grupos

m3step <- ' #
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  log_pathogen5+RAIZ_generalistsaprotro1+RAIZ_AMF+CP2_soilqual2+sqrt_avdiam
log_pathogen5 ~ CP2_soilqual2+ sqrt_avdiam
RAIZ_generalistsaprotro1 ~ CP2_soilqual2+ sqrt_avdiam
RAIZ_AMF ~ CP2_soilqual2+ sqrt_avdiam
RAIZ_AMF~~RAIZ_generalistsaprotro1
'
fit.f3step <- sem(m3step, data=dfcontrol)
summary(fit.f3step, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f3step, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) 

      # sin correlaciones entre fungal groups. RMSEA = 0.2. el resto ok
      # correlacion entre AMF y pathogen . RMESA = 0.31..idem
      # correlacion entre AMF y sapro . RMESA= 0.16..idem

#####
# probar
# a) el m1step pero con patogen composition 2.  RMSEA= 0.23
# b)el m1step con Pathogen comp y AMF 2 (no correl among Fung)  RMESA = 0 ( lo otro bien)  AIC = -412.12 
# b.1) si a;ado correlacion pathoge y AMF baja el AIC en comparacion con el anterior  AIC=-410.11
# b.2) si a;ado en lugar de la ant. corre AMF -sapro. AIC = -410.4

#### EL ELEGIDO es m4step

m4step <- ' #
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+RAIZ_generalistsaprotro1+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
RAIZ_generalistsaprotro1 ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.f4step <- sem(m4step, data=dfcontrol)
summary(fit.f4step, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f4step, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) 

#### hacer el mismo modelo pero para DROUGHT 

m4step.drou <- ' #
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+RAIZ_generalistsaprotro1+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
RAIZ_generalistsaprotro1 ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.f4step.drou <- sem(m4step.drou, data=dfdrought)       # data drought
summary(fit.f4step.drou, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f4step.drou, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) 

############
### Option 2. La que hice primero. Siguiendo el de Semchenko 2018. He hecho una linear regression 
# para cada variable y con la de menor AIC he hecho mis path analyses.
 
# El mejor seria el modelo 6  ###

#modelo escogido el numero 6# si incluyo correlaciones entre fgroups no FIT. 

modelfull6 <- ' #incluyendo soil properties
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+specialistsaprotro+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
specialistsaprotro ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.f6 <- sem(modelfull6, data=dfcontrol)
summary(fit.f6, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f6, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # fit ok!
#quitar pathogens para subir el R2. Lo hice pero no cambia practicamente nada los R2. 

####### He hecho el mismo modelo 6 pero poner saprotropgs rather than specialist saprotroph para el de control. 
### es mas parsimoniosos el modelo 6 (con specialist saprotrofos)
str(dfcontrol)

#same model (6) for drought
fit.f6d <- sem(modelfull6, data=dfdrought)
summary(fit.f6d, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f6d, c("chisq", "df", "pvalue", "cfi", "rmsea"))
# si en el modelo 6 quito pathogens..da lo mismo mas o menos


##############
### PATH ANALYSES PER PLANT FUNCTIONAL GROUP 

dfcontrol= read.csv("rawdata_control_fungal.csv")
dfdrought= read.csv("rawdata_drought_fungal.csv")
dfcontrol_grasses = subset(dfcontrol, Group == "Grass")
dfdrought_grasses = subset(dfdrought, Group =="Grass")

# He usado el modelo (m4step) de stepAIC

m4.grasses <- ' #
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+RAIZ_generalistsaprotro1+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
RAIZ_generalistsaprotro1 ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.f4grasses <- sem(m4.grasses, data=dfcontrol_grasses)  #data grasses
summary(fit.f4grasses, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.f4grasses, c("chisq", "df", "pvalue", "cfi", "rmsea","AIC")) # RMSEA = 0.1 .Lo usamos

######
# Usanod specialist saprotro...pero el stepAIC no escogio specialist sapro por eso voy a dejar el anterior
modelfullgrasses <- '  # ES EL MODELO 6 solo para grasses
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+specialistsaprotro+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
specialistsaprotro ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.fg <- sem(modelfullgrasses, data=dfcontrol_grasses)
summary(fit.fg, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.fg, c("chisq", "df", "pvalue", "cfi", "rmsea"))
# he hecho el de grasses sin "soil" y nada da significativo

#same model for drought
fit.fgd <- sem(modelfullgrasses, data=dfdrought_grasses)
summary(fit.fgd, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.fgd, c("chisq", "df", "pvalue", "cfi", "rmsea"))

#Herbs 

dfcontrol= read.csv("rawdata_control_fungal.csv")
dfcontrol_herbs = subset(dfcontrol, Group == "Herb")
dfdrought_herbs = subset(dfdrought, Group == "Herb")

#modelo completo (6) NO FIT. He quitado las correlaciones among groups 
modelfullherbs <- '  #  
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+specialistsaprotro+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
specialistsaprotro ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.fh <- sem(modelfullherbs, data=dfcontrol_herbs)
summary(fit.fh, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.fh, c("chisq", "df", "pvalue", "cfi", "rmsea"))

#modelo drought
fit.fhd <- sem(modelfullherbs, data=dfdrought_herbs)
summary(fit.fhd, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.fhd, c("chisq", "df", "pvalue", "cfi", "rmsea"))


#Legumes. (escogido el modelo legumes 2 porq el mod1 no fit RMESA=0.3)

dfcontrol= read.csv("rawdata_control_fungal.csv")
dfcontrol_legumes = subset(dfcontrol, Group == "Legume")
dfdrought_legumes = subset(dfdrought, Group == "Legume")

modelfulllegumes <- '  # ES EL MODELO 6 
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  Pathg_comp_2+specialistsaprotro+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
Pathg_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
specialistsaprotro ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
#fit.fl <- sem(modelfulllegumes, data=dfcontrol_legumes)
#summary(fit.fl, standardized=TRUE, fit.measures=T,rsq=T) 
#fitMeasures(fit.fl, c("chisq", "df", "pvalue", "cfi", "rmsea")) # no fit RMSEA# 

#el escogido es el que sigue#
modelfulllegumes2 <- '  # he quitado pathogens, pero tb podria ser sin sapro
#regressions ~, correlations ~~, latente =~  
log_shootmass ~  specialistsaprotro+AMF_comp_2+CP2_soilqual2+sqrt_avdiam
specialistsaprotro ~ CP2_soilqual2+ sqrt_avdiam
AMF_comp_2 ~ CP2_soilqual2+ sqrt_avdiam
'
fit.fl2 <- sem(modelfulllegumes2, data=dfcontrol_legumes)
summary(fit.fl2, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.fl2, c("chisq", "df", "pvalue", "cfi", "rmsea")) # fit RMSEA 0.08# 
# he hecho el modelo quitando soil p (tiene el AIC mas alto) pero el RMSEA da peor 0.1 

#modelo drought
fit.fl2d <- sem(modelfulllegumes2, data=dfdrought_legumes)
summary(fit.fl2d, standardized=TRUE, fit.measures=T,rsq=T) 
fitMeasures(fit.fl2d, c("chisq", "df", "pvalue", "cfi", "rmsea")) # RMSEA 0.5 

# he probado quitando sapro y dejando pathogens... es mejor el modelo con sapro!!! 

#### GRAPHS ABUNDANCE- RICHNESS etc

fdata= read.csv("rawdata_control_drought_fungal.csv") # NO BORRAR LOS OTROS (control or drought)
names(fdata)

#ABUNDANCE
fdata$TTO = factor(fdata$TTO, levels = c("well-watered", "drought")) # order panels


##PATHOGENS abundance
abpatho<-ggplot(fdata, aes(x=species, y=pathogenab, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ ylim(0,10)+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="white", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="Pathogens (%)", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
abpatho

#Saprotrofos abundance
absapro<-ggplot(fdata, aes(x=species, y=saproab, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ ylim(0,70)+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="white", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="Saprotrophs (%)", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
absapro
#AMF abundance
abamf<-ggplot(fdata, aes(x=species, y=AMFab, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ ylim(0,3)+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='Plant species', y="AMF (%)", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
abamf


names(fdata)
##PATHOGENS richness per plant species
ggplot(fdata, aes(x=species, y=Spatho, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ 
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='Plant species', y="Pathogens", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels

#Saprotrofos richness per plant species 
ggplot(fdata, aes(x=species, y=S_sapro, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ ylim(0,70)+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='Plant species', y="Saprotrophs ", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels

#AMF richness per plant species 
ggplot(fdata, aes(x=species, y=S_AMF, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ ylim(0,3)+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='Plant species', y="AMF", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels

## PLANT FUNCTIONAL GROUP  (mejor los graph de abund que de richness. Sobretodo AMF)
##PATHOGENS abund
names(fdata)
PFG_pathoab<-ggplot(fdata, aes(x=Group, y=pathogenab, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ 
  geom_point(pch = 21, position = position_jitterdodge())+
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=0, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="Pathogen", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
PFG_pathoab

#SAPRO AB
PFG_saproab<-ggplot(fdata, aes(x=Group, y=saproab, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ 
  geom_point(pch = 21, position = position_jitterdodge())+
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=0, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="Saprotroph", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) 
PFG_saproab

#AMF ab
PFG_amfab<-ggplot(fdata, aes(x=Group, y=AMFab, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ ylim(0,3)+
  geom_point(pch = 21, position = position_jitterdodge())+
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=0, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='Plant functional group', y="AMF", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) 
PFG_amfab

ggarrange(PFG_pathoab, PFG_saproab,PFG_amfab, labels = c("A", "B","C"),
          common.legend = TRUE, align=c("hv"),ncol = 1, nrow = 3) # allign graph
### NO3
str(fdata)
fdata$TTO = factor(fdata$TTO, levels = c("well-watered", "drought")) # order panels

ggplot(fdata, aes(x=species, y=NO3.mg.l., fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ 
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ Group, scales = "free", space = "free")+ #facet grid no repetida
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', fill="")+
  ylab(bquote(bold('NO3 ('*~mg~kg^-1*')')))+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) 

# GRAPH specialist vs generalist

datgs<-read.csv("specialist_generalist_ab.csv")
names(datgs)

datgs$TTO = factor(datgs$TTO, levels = c("well-watered", "drought")) # order panels


dpatho = subset(datgs, fgroup == "Pathogen")
dsapro = subset(datgs, fgroup == "Saprotroph")
dAMF = subset(datgs, fgroup == "AMF")

pat<-ggplot(dpatho, aes(x=plantgroup, y=relabund, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+ 
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ subgroup, scales = "free", space = "free")+
theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="white", size=11, angle=0, hjust = 0.5),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="Pathogens (%)", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
pat

sap<-ggplot(dsapro, aes(x=plantgroup, y=relabund, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ subgroup, scales = "free", space = "free")+
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="white", size=11, angle=0, hjust = 0.5),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="Saprotrophs (%)", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
sap

am<-ggplot(dAMF, aes(x=plantgroup, y=relabund, fill=TTO)) +
  geom_boxplot(aes(fill=TTO))+
  geom_point(pch = 21, position = position_jitterdodge())+
  facet_grid(. ~ subgroup, scales = "free", space = "free")+
  theme_bw()+
  scale_fill_manual(values=c("lightblue","darkgoldenrod2"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=0, hjust = 0.5),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='', y="AMF (%)", fill="")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels
am

ggarrange(pat,sap,am, labels = c("A", "B", "C", "D"),
          ncol = 1, nrow = 3,
          common.legend = TRUE, align=c("hv")) 

#######################



# RELAIMPO _Groemping  ALL TRAITS (Leaf- root- soil )



dataF <- read.csv ("traits_fungal.csv")
names(dataF)
Mutualists2.RI <-lm(Mutualists2 ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)

metrics_Mab <- calc.relimp(Mutualists2.RI, type = c("lmg","pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)

metrics_Mab@pmvd = metrics_Mab@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_Mab=data.frame(variables=names(metrics_Mab@pmvd),varexpl.percent=metrics_Mab@pmvd) # nueva tabla

Tab_metrics_Mab$variables <- factor(Tab_metrics_Mab$variables,
                                levels=c("SRSA","diameter","LDMC","soilC","soilN","pH", "rootmass",
                                         "WSA","shootmass","SRL","SLA","RTD"))  # order panels


REL_MUT_AB<- ggplot(data=Tab_metrics_Mab, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Mutualists abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_AB
#########

## relaimpo S_Mutual 
S_Mutual.RI <-lm(S_Mutual ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)

metrics_MS <- calc.relimp(S_Mutual.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
par(cex.axis = 0.8)
plot(metrics_MS)

metrics_MS@pmvd = metrics_MS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_MS=data.frame(variables=names(metrics_MS@pmvd),varexpl.percent=metrics_MS@pmvd) # nueva tabla

Tab_metrics_MS$variables <- factor(Tab_metrics_MS$variables,
                                levels=c("pH", "WSA","diameter","SRSA", "SRL","LDMC","SLA","shootmass","soilN",
                                         "soilC","RTD","rootmass"))  # order panels


REL_MUT_MS<-ggplot(data=Tab_metrics_MS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Mutualists richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_MS

#########

## relaimpo PCO1_Mutual

PCo1_Mutual.RI <-lm(PCo1_Mutual ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_MC <- calc.relimp(PCo1_Mutual.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_MC@pmvd = metrics_MC@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_MC=data.frame(variables=names(metrics_MC@pmvd),varexpl.percent=metrics_MC@pmvd) # nueva tabla

Tab_metrics_MC$variables <- factor(Tab_metrics_MC$variables,
                                   levels=c("rootmass","soilN","soilC","SLA","SRL","shootmass","WSA","LDMC","SRSA","diameter",
                                            "pH","RTD"))  # order panels

REL_MUT_MC<-ggplot(data=Tab_metrics_MC, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Mutualists composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_MUT_MC

############
#RELAIMPO SAPRO_AB

Saprotrophs.RI <-lm(Saprotrophs2 ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_Sab <- calc.relimp(Saprotrophs.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_Sab@pmvd = metrics_Sab@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_Sab=data.frame(variables=names(metrics_Sab@pmvd),varexpl.percent=metrics_Sab@pmvd) # nueva tabla

Tab_metrics_Sab$variables <- factor(Tab_metrics_Sab$variables,
                                levels=c("rootmass","SRL","pH","WSA","diameter","shootmass", "LDMC",
                                         "soilN","soilC","SLA","SRSA","RTD"))  # order panels

REL_SAP_AB<- ggplot(data=Tab_metrics_Sab, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Saprotrophs abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_AB

#RELAIMPO SAPRO_S

S_Saprotrophs.RI <-lm(S_Sapro ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_SS <- calc.relimp(S_Saprotrophs.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_SS@pmvd = metrics_SS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_SS=data.frame(variables=names(metrics_SS@pmvd),varexpl.percent=metrics_SS@pmvd) # nueva tabla

Tab_metrics_SS$variables <- factor(Tab_metrics_SS$variables,
                                   levels=c("WSA","SRSA","LDMC","shootmass","rootmass","SRL", "pH",
                                            "RTD","diameter","SLA","soilC","soilN"))  # order panels

REL_SAP_SS<- ggplot(data=Tab_metrics_SS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Saprotrophs richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_SS

#RELAIMPO SAPRO_PCOA1
names(dataF)
PCo1_Saprotrophs.RI <-lm(PCo1_Sapro ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_SC <- calc.relimp(PCo1_Saprotrophs.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_SC@pmvd = metrics_SC@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_SC=data.frame(variables=names(metrics_SC@pmvd),varexpl.percent=metrics_SC@pmvd) # nueva tabla

Tab_metrics_SC$variables <- factor(Tab_metrics_SC$variables,
                                   levels=c("SRL","shootmass","WSA","rootmass","RTD","pH", "soilC",
                                            "soilN","SRSA","diameter","SLA","LDMC"))  # order panels

REL_SAP_SC<- ggplot(data=Tab_metrics_SC, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Saprotrophs composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_SAP_SC

#RELAIMPO PATHO_AB

Pathogens.RI <-lm(Pathogens2 ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_Pab <- calc.relimp(Pathogens.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_Pab@pmvd = metrics_Pab@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_Pab=data.frame(variables=names(metrics_Pab@pmvd),varexpl.percent=metrics_Pab@pmvd) # nueva tabla

Tab_metrics_Pab$variables <- factor(Tab_metrics_Pab$variables,
                                    levels=c("soilC","SRL","SRSA","shootmass","diameter","SLA", "RTD",
                                             "rootmass","WSA","LDMC","pH","soilN"))  # order panels

REL_PAT_AB<- ggplot(data=Tab_metrics_Pab, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Pathogens abundance', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_PAT_AB

#RELAIMPO PATHO_S

S_Pathogens.RI <-lm(S_Patho ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_PS <- calc.relimp(S_Pathogens.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_PS@pmvd <- metrics_PS@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_PS <- data.frame(variables=names(metrics_PS@pmvd),varexpl.percent=metrics_PS@pmvd) # nueva tabla

Tab_metrics_PS$variables <- factor(Tab_metrics_PS$variables,
                                   levels=c("WSA","LDMC","SRSA","SRL","shootmass","diameter","SLA", "RTD",
                                            "soilC","soilN","rootmass","pH"))  # order panels

REL_PAT_PS<- ggplot(data=Tab_metrics_PS, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Pathogens richness', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_PAT_PS

#RELAIMPO PATHO_PCOA1
names(dataF)
PCo1_Pathogens.RI <-lm(PCo1_Patho ~shootmass+LDMC+SLA+rootmass+diameter+RTD+SRL+SRSA+soilN+soilC+WSA+pH, data = dataF)
metrics_PC <- calc.relimp(PCo1_Pathogens.RI, type = c("pmvd"), rela = TRUE)# Hay mas indices e.g "first", "last","lmg". etc 
#rela= TRUE (forced to percentages)
metrics_PC@pmvd = metrics_PC@pmvd*100  # para q me de la variance explicada en porcentaje de una vez
Tab_metrics_PC=data.frame(variables=names(metrics_PC@pmvd),varexpl.percent=metrics_PC@pmvd) # nueva tabla

Tab_metrics_PC$variables <- factor(Tab_metrics_PC$variables,
                                   levels=c("WSA","shootmass","rootmass","diameter","RTD","SRL", "pH",
                                            "SRSA","LDMC","soilN","soilC","SLA"))  # order panels

REL_PATHO_PC<- ggplot(data=Tab_metrics_PC, aes(x=variables, y=varexpl.percent)) +
  geom_bar(stat="identity", fill="black")+ coord_flip()+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=10, color="black"))+
  labs(x='Predictors of Pathogens composition', y="Relative importance (%)")+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  geom_text(x=1, y=4000, label="R2 = ###%")
REL_PATHO_PC

ggarrange (REL_PAT_AB, REL_PAT_PS, REL_PATHO_PC,REL_SAP_AB, REL_SAP_SS,REL_SAP_SC, REL_MUT_AB, REL_MUT_MS, REL_MUT_MC,
           labels = c("A", "B","C","D","E","F","G", "H", "I"),
           common.legend = FALSE, align = c("hv"))

##NMDS (10 species)

d10NMDS = read.csv("ESV_10species_grasslands.csv")
str(d10NMDS)


#Identify columns ()
data1 = d10NMDS [,6:1964] #solo las ESV o lo que quiero grafica
data2 = d10NMDS [,1:3]    # las variables descriptivas o factores 

#ordination by NMDS
NMDS <-metaMDS(data1, distance = "bray", k = 2) # k = axis number, here check the stres (low 0.2 good)

#graph
co=c("black","blue","darkgoldenrod1","darkgoldenrod4","darkolivegreen1","darkolivegreen4", 
     "firebrick1","firebrick4","gray48","gray13","aliceblue","cadetblue","darkgoldenrod1","darkgoldenrod4","darkolivegreen1","darkolivegreen4", 
     "firebrick1","firebrick4","gray48","gray13","gray48","gray13")
co=c("black","blue","darkgoldenrod1","darkgoldenrod4","darkolivegreen1","darkolivegreen4", 
     "firebrick1","firebrick4","gray48","gray13","aliceblue","cadetblue","darkgoldenrod1","darkgoldenrod4","darkolivegreen1","darkolivegreen4", 
     "firebrick1","firebrick4","gray48","gray13","gray48","gray13")

shape = c(1:20)

plot(NMDS$points,col=co[data2$SpecieTTO],pch = shape[data2$SpecieTTO],
     cex=1.2, main="Fungal community composition",  xlab = "axis 1", ylab = "axis 2")

plot(NMDS$points, col=co[data2$SpecieTTO],pch = shape[data2$SpecieTTO],
     cex=1.2, main="Fungal community composition",  xlab = "NMDS 1", ylab = "NMDS 2",
     xlim=c(-1.5, 1.0), ylim=c(-0.6, 0.6))

#Connect the points that belong to the same treatment with ordispider
ordihull(NMDS, groups = data2$SpecieTTO,  label = TRUE)

#Add legend
txt <- c("Grassland","Marsh")
legend('topleft', txt , pch=c(18,16),col=c("red","blue"),cex=1, bty = "y")

##################### PERMANOVA
#Bootstrapping and testing for differences between the groups
fit <- adonis(data1 ~ SpecieTTO, data=data2, permutations=999, method="bray")
fit


#####################
#Check assumption of homogeneity of multivariate dispersion

# here si variation in the specieTTO 1 has to be similar to the variation SpTTO2 
#higher than 0.05 is good.WE asume there is homogeneity 
distances_data <- vegdist(data1)
anova(betadisper(distances_data, data2$SpecieTTO))  ### We need to correct// how? 

### NMDS FUNGAL GROUPS COMPOSITION


