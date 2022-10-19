################
# (0) PACKAGES #
################

library(readxl)
library(patchwork)
library(FD)
library(ggplot2)
library(fastDummies)
library(lme4)
library(MuMIn)
library(rsq)

#######################################################################
# (1) IMPORT DATASET (table of diet/killifish body size class         #
#                    from Supporting information in Arim et al. 2010) #
#######################################################################

Prey_SizeClass <- read_excel("Prey_SizeClass.xlsx")
Prey_SizeClass<-as.data.frame(Prey_SizeClass)
View(Prey_SizeClass)

#######################
# (2) MATRIX ASSEMBLY #
#######################

# Matrix Q: spp (Prey) / traits
M.Q<-Prey_SizeClass[,c("Prey_item", "TG_Arim", "Functional_id", "TG_Laufer")] # Prey items categorized as BASAL, PRIMARY CONSUMER or PREDATOR 
rownames(M.Q)<-M.Q[,"Prey_item"]
M.Q<-M.Q[,-1]

# Matrix R: Body size classes -BSC- / mean BS
M.R<-as.numeric(colnames(Prey_SizeClass[,5:24])) # Mean body size value for each killifish BSC 
names(M.R)<-colnames(Prey_SizeClass[,5:24])

M.R2<-cbind(M.R, M.R^2) # Mean body size value (simple and quadratic) for each killifish BSC 
colnames(M.R2)[1:2]<-c("BS","BS2")
dim(M.R2)

# Matrix L: Prey / BSC
M.L<-Prey_SizeClass[,5:24] # Prey abundance in each killifish BSC
rownames(M.L)<-rownames(M.Q) 
M.L<-t(M.L)

################
# (3) ANALYSES #
################

#-------------------------------------------------------------#
# (3.1) TROPHIC GROUPS REPRESENTATION IN KILLIFISH DIET -CWM- #
#-------------------------------------------------------------#

cwm<-functcomp(M.Q[-c(4, 16, 21, 39, 48:49, 51),], M.L[,-c(4, 16, 21, 39, 48:49, 51)], CWM.type = "all")
cwm
cwm_R<-data.frame(cbind(cwm,M.R))
cwm_R

mA0.1<-glm(cwm_R[,"TG_Arim_0"]~cwm_R[,"M.R"])
mA0.2<-glm(cwm_R[,"TG_Arim_0"]~cwm_R[,"M.R"] + I(cwm_R[,"M.R"]^2))
summary(mA0.1)
summary(mA0.2)
anova(mA0.1, mA0.2, test = "Chisq") # Models are equivalent: retain mA0.1
plotA0<-ggplot(cwm_R, aes(M.R, TG_Arim_0))+
               geom_point(size=5, shape=19, color="navy")+
               geom_smooth(aes(x=M.R, y=TG_Arim_0), method='lm', formula= y~x) +
               labs(x="Predator mean body size (mm)", y="Community weighted mean (CWM)", title="Basal prey items")+
               theme_classic(base_size = 15) +
               geom_hline(yintercept=0, linetype="dashed", 
                         color = "black", size=1.5) +
               theme(axis.text=element_text(size=16),
                    axis.title=element_text(size=16,face="bold"))+
               theme(plot.title = element_text(hjust = 0.5))

mA1.1<-glm(cwm_R[,"TG_Arim_1"]~cwm_R[,"M.R"])
mA1.2<-glm(cwm_R[,"TG_Arim_1"]~cwm_R[,"M.R"] + I(cwm_R[,"M.R"]^2))
summary(mA1.1)
summary(mA1.2)
anova(mA1.1, mA1.2, test = "Chisq") #  Models are equivalent: retain mA1.1
plotA1<-ggplot(cwm_R, aes(M.R, TG_Arim_1))+
               geom_point(size=5, shape=19, color="navy")+
               geom_smooth(aes(x=M.R, y=TG_Arim_1), method='lm', formula= y~x) +
               labs(x="Predator mean body size (mm)", y="", title="Primary consumers prey items")+
               theme_classic(base_size = 15) +
               theme(axis.text=element_text(size=16),
                    axis.title=element_text(size=16,face="bold"))+
               theme(plot.title = element_text(hjust = 0.5))

mA2.1<-glm(cwm_R[,"TG_Arim_2"]~cwm_R[,"M.R"])
mA2.2<-glm(cwm_R[,"TG_Arim_2"]~cwm_R[,"M.R"] + I(cwm_R[,"M.R"]^2))
summary(mA2.1)
summary(mA2.2)
anova(mA2.1, mA2.2, test = "Chisq")            #  Models are equivalent: retain mA2.1
plotA2<-ggplot(cwm_R, aes(M.R, TG_Arim_2))+
               geom_point(size=5, shape=19, color="navy")+
               geom_smooth(aes(x=M.R, y=TG_Arim_2), method='lm', formula= y~x) +
               labs(x="Predator mean body size (mm)", y="", title="Predator prey items")+
               theme_classic(base_size = 15) +
               theme(axis.text=element_text(size=16),
                    axis.title=element_text(size=16,face="bold")) +
               theme(plot.title = element_text(hjust = 0.5))

plotA0 + plotA1 + plotA2


#-----------------------------#
# (3.2) CATS REGRESSION MODEL #
#-----------------------------#

# (3.2.1) Estimates the body size of each prey item 
Animal_2008_dieta<-read_excel("Animal_2008_dieta.xlsx")
Animal_2008_dieta<-as.data.frame(Animal_2008_dieta)
View(Animal_2008_dieta)

BS.prey<-prey.bs.diet(M=Animal_2008_dieta, col.prey=3, col.bs=5, col.diet=9)
BS.prey<-BS.prey[-31,]

quantile.prey<-quantile(BS.prey[,"mean.bs"], c(0.1, 0.5, 0.9), na.rm = T) # Vector with the quantile 0.1, 0.5 and 0.9 from the mean body size of prey items

# (3.2.2) Prepares matrix M.Q 
M.Q2<-M.Q[-c(4, 16, 21, 39, 48:49, 51),] # Eliminates NA rows
M.Q22<-M.Q2[order(rownames(M.Q2)),]      # Sorts alphabetically the rows to match the order of the rows in the matrix BS.prey

# (3.2.3) Trophic group categories translated into dummy variables (necessary in function EXPANDE)
traits_dummy22<-dummy_cols(M.Q22[,1], ignore_na = TRUE, remove_selected_columns = TRUE)  
colnames(traits_dummy22)
traits_dummy22<-traits_dummy22[,-2] # It is necessary to remove one column (trophic group) to avoid overadjusment (column removed: PRIMARY CONSUMERS)
traits_dummy22

# (3.2.4) Prepares matrix M.L 
M.L22<-M.L[,-c(4, 16, 21, 39, 48:49, 51)] # Eliminates prey items with trophic group=NA 

M.L22<-M.L22[,order(colnames(M.L22))] # Sorts alphabetically the columns to match the order of the rows in the matrix BS.prey

colnames(M.L22)==rownames(M.Q22)    # Names match
colnames(M.L22)==rownames(BS.prey)  # Names match

# (3.2.5) Prepares the matrix to perform the CATS regression model (uses function expande) 
M5<-expande(L = M.L22, Q = cbind(traits_dummy22, BS.prey[,"mean.bs"]), R = (M.R2))
colnames(M5)<-c("B", "P", "mean.bs", "BS", "BS2", "B_BS", "B_BS2", "P_BS", "P_BS2", "mean.bs_BS", "mean.bs_BS2", "N")

M6<-cbind(M5, rep(1:44, 20)) # Includes a column with the id of prey items
colnames(M6)[13]<-"prey.id"

M7<-cbind(M6,rep(1:20,  each=44)) # Includes a column with the id of killifish body size class
colnames(M7)[14]<-"BSC.id"

M8<-cbind(M7, M7[,"B"]*M7[,"mean.bs"], M7[,"P"]*M7[,"mean.bs"]) # Includes the interaction between prey trophic group and prey body size
colnames(M8)[15:16]<-c("B_mean.bs", "P_mean.bs")

# (3.2.6) CATS regression (mixed) model
m7.expande.6<-glmer(N ~ mean.bs + BS + BS2 + B + P +  B_BS + B_BS2 + P_BS + P_BS2 + mean.bs_BS + mean.bs_BS2 + 
                      mean.bs:B + B_mean.bs:BS + P_mean.bs:BS +
                      (1|prey.id) + (1|BSC.id), 
                    offset=rep(apply(M.L22,2,sum)/sum(apply(M.L22,2,sum)),20), 
                    family = "poisson", data = M8)

summary(m7.expande.6)
rsq(m7.expande.6)
coeff.m7.expande.6<-coeffs(m7.expande.6) 


#------------------------#
# (3.3) CATS MODEL PLOTS #
#------------------------#

plot.basales<-ggplot() +
              xlim(9.9, 36.2) +
              geom_function(fun = function(x) (coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["10%"] + coeff.m7.expande.6["B"] + coeff.m7.expande.6["mean.bs:B"]*quantile.prey["10%"]) + 
                              (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["10%"] + coeff.m7.expande.6["B_BS"] + coeff.m7.expande.6["BS:B_mean.bs"]*quantile.prey["10%"])*x +
                              (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["10%"] + coeff.m7.expande.6["B_BS2"])*x^2, 
                            colour = "springgreen4", linetype = 2, size = 1) +
              
              geom_function(fun = function(x) (coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["50%"] + coeff.m7.expande.6["B"] + coeff.m7.expande.6["mean.bs:B"]*quantile.prey["50%"]) + 
                              (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["50%"] + coeff.m7.expande.6["B_BS"] + coeff.m7.expande.6["BS:B_mean.bs"]*quantile.prey["50%"])*x +
                              (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["50%"] + coeff.m7.expande.6["B_BS2"])*x^2,
                            colour = "springgreen4", linetype = 1, size = 1) +
              
              geom_function(fun = function(x) (coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["90%"] + coeff.m7.expande.6["B"] + coeff.m7.expande.6["mean.bs:B"]*quantile.prey["90%"]) + 
                              (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["90%"] + coeff.m7.expande.6["B_BS"] + coeff.m7.expande.6["BS:B_mean.bs"]*quantile.prey["90%"])*x +
                              (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["90%"] + coeff.m7.expande.6["B_BS2"])*x^2,
                            colour = "springgreen4", linetype = 2, size = 2) +
              labs(x="Predator mean body size (mm)", y="Selection coefficient", title="Basal prey items") +
              
              theme_classic(base_size = 20) +
              theme(axis.text=element_text(size=16),
                    axis.title=element_text(size=16,face="bold")) +
              theme(plot.title = element_text(hjust = 0.5, size=16)) +
              geom_hline(yintercept=0, size=1)

plot.primary.consumers<-ggplot() +
                        xlim(9.9, 36.2) +
                        geom_function(fun = function(x) coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["10%"] +  # Intercepto para grupo de referencia (PC)
                                        (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["10%"])*x +      # termino lineal para grupo de referencia (PC)
                                        (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["10%"])*x^2, 
                                      colour = "blue", linetype = 2, size = 1) +
                        
                        geom_function(fun = function(x) coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["50%"] + # Intercepto para grupo de referencia (PC)
                                        (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["50%"])*x +     # termino lineal para grupo de referencia (PC)
                                        (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["50%"])*x^2,
                                      colour = "blue", linetype = 1, size = 1) +
                        
                        geom_function(fun = function(x) coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["90%"] + # Intercepto para grupo de referencia (PC)
                                        (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["90%"])*x +     # termino lineal para grupo de referencia (PC)
                                        (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["90%"])*x^2,
                                      colour = "blue", linetype = 2, size = 2) +
                        labs(x="Predator mean body size (mm)", y="", title="Primary consumer prey items") +
                        
                        theme_classic(base_size = 20) +
                        theme(axis.text=element_text(size=16),
                              axis.title=element_text(size=16,face="bold")) +
                        theme(plot.title = element_text(hjust = 0.5, size=16)) +
                        geom_hline(yintercept=0, size=1)

plot.predators<-ggplot() +
                  xlim(9.9, 36.2) +
                  geom_function(fun = function(x) (coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["10%"] + coeff.m7.expande.6["P"]) + # Intercepto para grupo P
                                  (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["10%"] + coeff.m7.expande.6["P_BS"] + coeff.m7.expande.6["BS:P_mean.bs"]*quantile.prey["10%"])*x +     # Termino lineal para grupo P
                                  (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["10%"] + coeff.m7.expande.6["P_BS2"])*x^2, 
                                colour = "purple", linetype = 2, size = 1) +
                  
                  geom_function(fun = function(x) (coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["50%"] + coeff.m7.expande.6["P"]) + # Intercepto para grupo P
                                  (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["50%"] + coeff.m7.expande.6["P_BS"] + coeff.m7.expande.6["BS:P_mean.bs"]*quantile.prey["50%"])*x +     # Termino lineal para grupo P
                                  (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["50%"] + coeff.m7.expande.6["P_BS2"])*x^2,
                                colour = "purple", linetype = 1, size = 1) +
                  
                  geom_function(fun = function(x) (coeff.m7.expande.6["(Intercept)"] + coeff.m7.expande.6["mean.bs"]*quantile.prey["90%"] + coeff.m7.expande.6["P"]) + # Intercepto para grupo P
                                  (coeff.m7.expande.6["BS"] + coeff.m7.expande.6["mean.bs_BS"]*quantile.prey["90%"] + coeff.m7.expande.6["P_BS"] + coeff.m7.expande.6["BS:P_mean.bs"]*quantile.prey["90%"])*x +     # Termino lineal para grupo P
                                  (coeff.m7.expande.6["BS2"] + coeff.m7.expande.6["mean.bs_BS2"]*quantile.prey["90%"] + coeff.m7.expande.6["P_BS2"])*x^2,
                                colour = "purple", linetype = 2, size = 2) +
                  labs(x="Predator mean body size (mm)", y="", title="Predator prey items") +
                  
                  theme_classic(base_size = 20) +
                  theme(axis.text=element_text(size=16),
                        axis.title=element_text(size=16,face="bold")) +
                  theme(plot.title = element_text(hjust = 0.5, size=16)) +
                  geom_hline(yintercept=0, size=1)

plot.basales + plot.primary.consumers + plot.predators


##############################
# (4) SUPPORTING INFORMATION #
##############################

#-----------------------------------------#
# (4.1) PLOT: PREY BODY SIZE DISTRIBUTION #
#-----------------------------------------#
Animal_2008_dieta_BS_TG <- read_excel("Animal_2008_dieta_BS_TG.xlsx")
View(Animal_2008_dieta_BS_TG)
str(Animal_2008_dieta_BS_TG)
Animal_2008_dieta_BS_TG<- as.data.frame(Animal_2008_dieta_BS_TG)
Animal_2008_dieta_BS_TG[,"Trophic_group"]<-as.factor(Animal_2008_dieta_BS_TG[,"Trophic_group"])
str(Animal_2008_dieta_BS_TG)


bs.BASAL<-Animal_2008_dieta_BS_TG[which(Animal_2008_dieta_BS_TG[,"Trophic_group"]=="Basal"),]
quantile(bs.BASAL[,"Biovol"], probs = c(0.01, 0.99))
bs.PRIMARY.CONSUMERS<-Animal_2008_dieta_BS_TG[which(Animal_2008_dieta_BS_TG[,"Trophic_group"]=="Primary_consumer"),]
quantile(bs.PRIMARY.CONSUMERS[,"Biovol"], probs = c(0.01, 0.99))
bs.PREDATORS<-Animal_2008_dieta_BS_TG[which(Animal_2008_dieta_BS_TG[,"Trophic_group"]=="Predator"),]
quantile(bs.PREDATORS[,"Biovol"], probs = c(0.01, 0.99))

Animal_2008_dieta_BS_TG_quantil_0.01_0.99<-rbind(bs.BASAL[which(bs.BASAL[,"Biovol"]>0 & bs.BASAL[,"Biovol"]<100),],
                                                 bs.PRIMARY.CONSUMERS[which(bs.PRIMARY.CONSUMERS[,"Biovol"]>0.01 & 
                                                                              bs.PRIMARY.CONSUMERS[,"Biovol"]<240),],
                                                 bs.PREDATORS[which(bs.PREDATORS[,"Biovol"]>0.02 & 
                                                                      bs.PREDATORS[,"Biovol"]<100),])

ggplot(Animal_2008_dieta_BS_TG_quantil_0.01_0.99[which(Animal_2008_dieta_BS_TG_quantil_0.01_0.99[,"Trophic_group"]==c("Primary_consumer","Predator")),], 
       aes(x=Biovol, fill=Trophic_group)) +
  geom_histogram( binwidth = 10) +
  theme_bw() +
  theme(legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 15)) +
  xlab(bquote("Estimated prey body size"~(mm^3))) +
  ylab("Frecuency") +
  
  facet_wrap(~Trophic_group, labeller = as_labeller(c("Predator"="Predator", "Primary_consumer"="Primary consumer"))) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        axis.title.x.top = element_text(size=15))


#------------------------------------------------------------------------------#
# (4.2) PLOT: PREY TROPHIC GROUP DISTRIBUTION AMONG PREDATOR BODY SIZE CLASSES #
#------------------------------------------------------------------------------#

Dieta_Arim_2010 <- read_excel("Dieta_Arim_2010.xlsx")
Dieta_Arim_2010<-as.data.frame(Dieta_Arim_2010)
Dieta_Arim_2010[,2]<-as.numeric(Dieta_Arim_2010[,2])
Dieta_Arim_2010<-Dieta_Arim_2010[-c(4, 16, 21, 39, 48, 49, 51),]
colnames(Dieta_Arim_2010)[c(15,22)]<-c("16.9", "36.2")
Prey_TG_distribution<-prey.proportion.BSC(M=Dieta_Arim_2010, TG.in=2, BSC.in=3:22)
BSC.labels<-colnames(Dieta_Arim_2010)[c(3:22)]

# Abundance of each TG in each BSC
plot.prey.abundance.BSC<-ggplot(Prey_TG_distribution, aes(fill=TG, y=Abundance, x=BSC)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Predator body size class (mm)") +
  labs(fill="Trophic group") +
  scale_fill_discrete(labels=c("Basal", "Primary consumer", "Predator"), 
                      type= c("springgreen4","blue", "purple")) +
  theme(axis.text.x = element_text(size = 14, angle = 60, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  theme(legend.text=element_text(size=10), 
        legend.title =element_text(size=12), legend.position = "none" ) +
  scale_x_discrete(limit=BSC.labels)

# Richness of each TG in each BSC
plot.prey.richness.BSC<-ggplot(Prey_TG_distribution, aes(fill=TG, y=Richness, x=BSC)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  xlab("Predator body size class (mm)") +
  labs(fill="Trophic group") +
  scale_fill_discrete(labels=c("Basal", "Primary consumer", "Predator"), 
                      type= c("springgreen4","blue", "purple")) +
  theme(axis.text.x = element_text(size = 14, angle = 60, hjust=1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  theme(legend.text=element_text(size=10), 
        legend.title =element_text(size=12), legend.position = c(0.2, 0.83)) +
  scale_x_discrete(limit=BSC.labels)


plot.prey.richness.BSC + plot.prey.abundance.BSC # 6x12


#---------------#
# END OF SCRIPT #
#---------------#


