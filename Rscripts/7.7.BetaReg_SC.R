# Beta regression on Selection Coefficients
library(betareg)
library(colorspace)
library(reshape2)
library(ggplot2)
source("Rscripts/BetaEffectSize.R")
source("Rscripts/label_scientific.R")
colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")

sccol<-c("#7071FF", "#3F4CE1")


# Read the data
SC<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F,row.names = 1)
SC<-SC[,c("pos", "EstSC")]
#Combine the shaped data and SC
betar<-read.csv("Output/BetaReg/BetaRegData.csv", stringsAsFactors = F, row.names = 1)
betaDF<-merge(betar, SC, by="pos")

################
#2. run Beta Regression -apply the same model as MF for comparison
mod1 <- betareg(EstSC ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2+NS1+NS2++NS3+NS4A+NS5A+NS5B+Shape, data = betaDF[betaDF$Stop == 0,])
AIC(mod1) #-91163.86
summary(mod1)
#              Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -6.538e+00  2.935e-02 -222.779  < 2e-16 ***
#t           -2.388e-01  3.402e-02   -7.019 2.23e-12 ***
#c           -5.089e-01  3.082e-02  -16.510  < 2e-16 ***
#g           -5.807e-01  3.380e-02  -17.184  < 2e-16 ***
#CpG          6.734e-02  1.173e-02    5.740 9.44e-09 ***
#Nonsyn       7.296e-01  2.803e-02   26.028  < 2e-16 ***
#bigAAChange  1.302e-01  1.050e-02   12.397  < 2e-16 ***
#Core         6.012e-02  2.223e-02    2.704 0.006846 ** 
#E1          -9.207e-02  2.132e-02   -4.318 1.57e-05 ***
#HVR1        -4.627e-01  6.069e-02   -7.623 2.48e-14 ***
#E2          -1.086e-01  1.847e-02   -5.879 4.13e-09 ***
#NS1         -1.727e-01  3.294e-02   -5.242 1.58e-07 ***
#NS2         -1.401e-01  2.094e-02   -6.690 2.23e-11 ***
#NS3         -6.983e-03  1.621e-02   -0.431 0.666722    
#NS4A         4.163e-05  3.274e-02    0.001 0.998986    
#NS5A        -5.772e-02  1.734e-02   -3.328 0.000876 ***
#NS5B         5.316e-02  1.814e-02    2.931 0.003380 ** 
#Shape       -1.003e-02  1.324e-02   -0.758 0.448483    
#t:Nonsyn     1.675e-01  3.610e-02    4.639 3.51e-06 ***
#c:Nonsyn     1.710e-01  3.307e-02    5.170 2.35e-07 ***
#g:Nonsyn    -8.965e-02  3.599e-02   -2.491 0.012750 *  


#remove NS3,Shape
mod2 <- update(mod1, ~. -NS3 -Shape)
AIC(mod2) # -91167.02
summary(mod2)

mod3 <- update(mod2, ~.-NS4A)
AIC(mod3) #-91168.98
summary(mod3)

#The best model from MF
mod3.2<-update(mod2, ~.-NS5A)
AIC(mod3.2)  #-91151.65  Not the best model

AIC(mod1,mod2,mod3,mod3.2)
#       df       AIC
#mod1   22 -91163.86
#mod2   20 -91167.02
#mod3   19 -91168.98
#mod3.2 19 -91151.65

sum3.2<-summary(mod3.2)

#A different model from MF is the best model (SC=MF+NS5A)
## Calculate effect size and save the the best model
effectsc<-BetaEffectSize(mod3)
write.csv(effectsc,"Output/BetaReg/BetaReg_SC_effectSize.csv")



# Create a combined effect size plot 

##effect size for MF & SC
emf<-read.csv("Output/BetaReg/BetaReg_MF_effectsize_SCmodel.csv",stringsAsFactors = F)
esc<-read.csv("Output/BetaReg/BetaReg_SC_effectSize.csv",stringsAsFactors = F, row.names = 1)

emf<-emf[,c("factor","Effect")]
colnames(emf)[2]<-"MF"
esc<-esc[,c("factor","Effect")]
colnames(esc)[2]<-"SC"

effects<-merge(emf, esc, by="factor", all =T)
effects<-effects[effects$factor!="(Intercept)",]
effects$factor<-factor(effects$factor, levels=emf$factor[-1])
#effects$percent<-effects$Effect*100

effects2<-melt(effects, id.var="factor")
colnames(effects2)[2:3]<-c("Type","Effect")

ggplot(effects2, aes(factor,Effect*100, group=Type, fill=Type, color=Type)) +
    geom_hline(yintercept = 0, size=0.5, color="gray40")+
    scale_fill_manual(values=c(paste0(colors2[4],"99"),paste0(sccol[1],"E6")), labels=c("MF","SC"))+
    scale_color_manual(values=c(colors2[4],sccol[2]))+
    geom_bar(stat="identity",position=position_dodge(width=.8), width = 1.2,  size=0.3)+
    theme_classic() +
    theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+

    labs(x="", y="Estimated effects (%)")+
    theme(legend.title = element_blank())

ggsave("Output/BetaReg/BetaReg.effects.MF-SC.pdf", width = 9.6, height = 5)


