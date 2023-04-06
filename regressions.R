install.packages("tidyverse")
library(ggplot2)


data <- read.csv("mydatafungalt.csv") 


data.control = data[ which(data$tto=='control'), ]
data.drought = data[ which(data$tto=='drought'), ]

str(data)
str(data.control)


x <- data$pathogen
y = data$shootmass

#transformate variable

logpathogen=log((data$pathogen)+5)
sqrtpathogen=sqrt((data$pathogen))

onlycontrol <- data[ which(data$tto=='control'), ]
onlydrought <- data[ which(data$tto=='drought'), ]

mod1patho = lm(shootmass ~ pathogen, data = onlycontrol)
summary(mod1patho)
plot(mod1patho) 

mod2patho = lm(shootmass ~ pathogen, data = onlydrought)
summary(mod2patho)
plot(mod2patho) 

# no correlation pathogen -shoot
# el resto de las correlaciones en infostat. 

#heteroces 0
#he probado con sqrt, log, raw.. y aqui el mejor para cad agrupo. 

##PATHOGEN AB
# el mejor fue log pathogen vs log shoot mass. pero no significancias
#generalists, specialists. No significancia

pathogengroup <- ggplot(data) +
  aes(log_pathogen5, log_shootmass, shape = Tto) +
  geom_point(aes(colour = Tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab("Log Pathogen abundance") +
  ylab("Log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
pathogengroup



aa = lm(shootmass_RII ~ log_pathogen5, data = data)
summary(aa) 
plot(aa) # 

aa<-ggplot(data, aes(x=log_pathogen5, y=shootmass_RII))+
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
aa


#SAPROTROFOS AB

# for total sapro sqrt sapro vs log shoot
#for genrelist sqrt
#for specialist . no efect en AOV ab. y no efecto aqui. 


saproAB <- ggplot(data) +
  aes(log_gensapro, log_shootmass, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" log gen Sapro abundance") +
  ylab("Log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
saproAB




bb = lm(shootmass_RII ~ sqrt_sapro, data = data)
summary(aa) 
plot(aa) # 

bb<-ggplot(data, aes(x=sqrt_sapro, y=shootmass_RII))+
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
bb

#AMF AB

# reuerda que no hay generalists. 


AMFAB <- ggplot(data) +
  aes(sqrt_speAMF, log_shootmass, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" sqrt spe abundance") +
  ylab("Log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
AMFAB


bb = lm(shootmass_RII ~ sqrt_speAMF, data = data)
summary(bb) 
plot(aa) # 

bb<-ggplot(data, aes(x=log_AMF5, y=shootmass_RII))+
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
bb

#PATHOGEN RICHNESS 
log_Spathogen= as.numeric(data$log_Spatho)  # convertir variable de Factor a numerico
str(log_Spathogen)

pathoS <- ggplot(data) +
  aes(log_Spathogen, log_shootmass, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" specialist_pathog_S") +
  ylab("Log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 0.5, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
pathoS


#SAPRO and AMF RICHNESS

# usar "S_sapro"
#usar S_AMF 

saproS <- ggplot(data) +
  aes(S_AMF, log_shootmass, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" S_AMF") +
  ylab("Log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 0.5, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
saproS

#saproS and amfS vs shoot RII

#usar S_sapro
#usar S_AMF

rr = lm(shootmass_RII ~ S_AMF, data = data) #log y sqrt con S_sapro y sin nada mas cerc a signif) 
summary(rr) 
plot(rr) # 

aa<-ggplot(data, aes(x=S_AMF, y=shootmass_RII))+
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
aa

###########
############
###########ROOT TRAITS VS MICROBIAL ABUNDANCE #
#
#1. Pearson correlation para identificar unas

# DIAMETRO. Correl with sqrt_AMF, Sqrt Spe sapro, S sapro, S all fun 


RDA = lm(RDA_RII ~ sqrt_AMF, data = data) #log y sqrt con S_sapro y sin nada mas cerc a signif) 
summary(RDA) #p=0.1 r2=0.01 

RDA2 = lm(RDA_RII ~ sqrt_spesapro, data = data) #log y sqrt con S_sapro y sin nada mas cerc a signif) 
summary(RDA2)  #p=0.6, r2= -0.007

RDA3 = lm(RDA_RII ~ S_sapro, data = data) #log y sqrt con S_sapro y sin nada mas cerc a signif) 
summary(RDA3) #p = 0.8, R2 = -0.008. 

RDA4 = lm(RDA_RII ~ S_allfung, data = data) #log y sqrt con S_sapro y sin nada mas cerc a signif) 
summary(RDA4) #p=0.7. r2= -0.007

RDA5 = lm(RDA_RII ~ sqrt_spesapro, data = data) #log y sqrt con S_sapro y sin nada mas cerc a signif) 
summary(RDA5) #p=0.6. r2= -0.007

#SRL . con sqrt_AMF, specialistsaprotro , AMf comp1 

SRL = lm(SRL_RII ~ sqrt_AMF, data=data)
summary(SRL) #p=0.3, r2= -0.0009. graph feo. 

SRL2 = lm(SRL_RII ~ sqrt_spesapro, data=data)
summary(SRL2) #p=0.1, r2= 0.01

SRL3 = lm(SRL_RII ~ AMF_comp_1, data=data)
summary(SRL3)#p=0.2, r2= 0.004

##SRSA     RAIZ_AMF, sqrt spe sapro, spesapro S 

SRSA = lm(SRSA_RII ~ sqrt_AMF, data=data)
summary(SRSA) #p=0.1, r2= -0.005

SRSA2 = lm(SRSA_RII ~ sqrt_spesapro, data=data)
summary(SRSA2) #p=0.04, r2= 0.02..graph ok 

SRSA3 = lm(SRSA_RII ~ specialist_sapro_S, data=data)
summary(SRSA3)#p=0.0006, r2= 0.08

#Root mass     raiz amf, s_amf, sapro comp 2

rootmass = lm(rootmass_RII ~ sqrt_AMF, data=data)
summary(rootmass) #p=0.1, r2= 0.006

rootmass2 = lm(rootmass_RII ~ S_AMF, data=data)
summary(rootmass2) #p=0.09, r2= 0.01

rootmass3 = lm(rootmass_RII ~ Sapro_comp2, data=data)
summary(rootmass3)#p=0.008, r2= 0.01

####graph 

# specialist_sapro_S
aa<-ggplot(data, aes(x=log_spesapro_S, y=SRL_RII))+
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
aa

#graph Root trait vs microbial group. droughtr and control separate. 

#### RDA (sqrt_RDA) 

# sqrt_AMF 
# pathogen            ok ..pero no tendency
# spepathogen         horrib
# sapro               ok...pero no tendency
# sqrt_sapro          ok... pero no tendency   f=2.24 p=0.1, r2= 0.01 /  f=1.35, p= 0.2 r2= 0.003
# generalistsaprotro  ok..idem
# S_AMF,              feo
# S_pathog            feo
#log_Spathogen        no se ve nada..
# sqrt_Ssapro         ok...ver summary    #####  f=1.96, p=0.1, r2= 0.008 / f= 4.78, p = 0.03, r2=0.03
# generalist_sapro_S  ok . similar al anterior
#specialist_sapro_S   feo 
#Fcomp_PC1            #meehh

Fcomp_PC1= as.numeric(data$Fcompos_PC1)  # convertir variable de Factor a numerico

RDA1.cont=  lm(sqrt_RDA ~ sqrt_sapro, data=data.control)
summary(RDA1.cont)
RDA1.drou = lm(sqrt_RDA ~ sqrt_sapro, data=data.drought)
summary(RDA1.drou)

RDA1.0 <- ggplot(data) +
  aes(sqrt_sapro, sqrt_RDA, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" sqrt_sapro") +
  ylab("sqrt_RDA") +
  scale_y_continuous(breaks = seq(0 , 0.5, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
RDA1.0




##### SRL    log SRL

#sqrt_AMF......Feo
#sqrt_sapro..... good   Control = f=0.2, p= 0.6, r2= -0.006 / drought= f=3.97, p = 0.04, r2=0.02
#sqrt_pathogen...no
#Fcompos_PC1....ok      f=1.38, , p= 0.24 r2= 0.003  /    f= 0.001, p=0.96, r2= -0.008
#Pathog_comp_1....maybe
#Sapro_Comp_1.......ok  f=1.73, p= 0.19, r2=0.006 /        f= 0.71, p= 0.4 r2= -0.002
#AMF_comp_1.......no
#sqrt_Ssapro......ok    f=0.2, p= 0.6, r2= -0.006  /        f=3.97, p= 0.04, r2= 0.02
#sqrt_Spatho.....no
#sqrt_Samf........no

SRL1.cont=  lm(log_SRL ~ sqrt_sapro, data=data.control)
summary(SRL1.cont)
SRL1.drou = lm(log_SRL ~ sqrt_Ssapro, data=data.drought)
summary(SRL1.drou)

SRL1.0 <- ggplot(data) +
  aes(sqrt_Ssapro, log_SRL, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" sqrt_Ssapro") +
  ylab("log_SRL") +
  scale_y_continuous(breaks = seq(0 , 0.5, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
SRL1.0

###

##### SRSA    sqrt SRSA

#sqrt_AMF......  feo graph.....
#sqrt_sapro..... OK ......f =0.02, p= 0.8 r2= -0.008 / f=0.9, p= 0.33 rr2= -0.0006
#sqrt_pathogen...OK........f=1.51 p = 0.2 r2= 0.004  /  f= 1.48, p = 0.2 r2= 0.004
#Fcompos_PC1....OK....... f= 0.41 p = 0.51 r2= -0.005 / f= 0.48 p = 0.48 r2= -0.004
#Pathog_comp_1....OK        f=1.37, p= 0.2 r2=0.003 /  f= 0.87, p = 0.3, r2= -0.001
#Sapro_Comp_1.......OK      f=1.76, p= 0.18 r2= 0.006 /  f= 0.24 p= 0.6 r2= -0.006
#AMF_comp_1....... NOOOOOO 
#sqrt_Ssapro...... OK       f=0.22 p=0.6 r2= -0.006 / f= 0.98, p=0.32 r2= 0.008
#sqrt_Spatho..... NO
#sqrt_Samf........no 

SRSA1.cont=  lm(sqrt_SRSA ~ sqrt_Ssapro, data=data.control)
summary(SRSA1.cont)
SRSA1.drou = lm(sqrt_SRSA ~ sqrt_Ssapro, data=data.drought)
summary(SRSA1.drou)


SRSA1.0 <- ggplot(data) +
  aes(sqrt_Ssapro, sqrt_SRSA, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" sqrt_Ssapro") +
  ylab("sqrt_SRSA") +
  scale_y_continuous(breaks = seq(0 , 0.5, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
SRSA1.0


#####root mass    sqrt root mass

#sqrt_AMF......  .....no
#sqrt_sapro.....  ok          f=0.66 p= 0.41 r2= -0.002 / 0.006, p = 0.9 r2= -0.008
#sqrt_pathogen... ok          f=0.002 p=0.62 r2= -0.006 / 0.45 p= 0.50 r2= -0.0047
#Fcompos_PC1.... ok           f=0.05 p= 0.82 r2= -0.008 / f= 0.0005 p= 0.9 r2= -0.008
#Pathog_comp_1.... ok         
#Sapro_Comp_1......ok
#AMF_comp_1....... no 
#sqrt_Ssapro...... ok         f=2.12 p= 0.14 r2 = 0.009 / f=7.36 p=0.007 r2 = 0.05
#sqrt_Spatho..... no
#sqrt_Samf........no

rootmass1.cont=  lm(sqrt_rootmass ~ sqrt_Ssapro, data=data.control) #regresion
summary(rootmass1.cont)
rootmass1.drou = lm(sqrt_rootmass ~ sqrt_Ssapro, data=data.drought)
summary(rootmass1.drou)

rootmass1.0 <- ggplot(data) +
  aes(sqrt_Samf, sqrt_rootmass, shape = tto) +
  geom_point(aes(colour = tto), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab(" sqrt_Samf") +
  ylab("sqrt_rootmass") +
  scale_y_continuous(breaks = seq(0 , 0.5, 1)) + 
  scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
rootmass1.0



#####Root N     log rootN

rootN= as.numeric(data$rootN_g) 
sqrt_rootN= sqrt(rootN)
log_rootN=log(rootN)
hist(log_rootN) 

#sqrt_AMF......  .....no
#sqrt_sapro.....       ok     
#sqrt_pathogen...    no    
#Fcompos_PC1....      ok      
#Pathog_comp_1....         
#Sapro_Comp_1......
#AMF_comp_1....... 
#sqrt_Ssapro......   ok       
#sqrt_Spatho.....   no
#sqrt_Samf........  no



