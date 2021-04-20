#GLM/Beta regression preparation:
library(tidyverse)
library(zoo)
library(purrr)
library(MASS)
library(betareg)
library(miscTools)
source("Rscripts/baseRscript.R")

####################
SC<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F,row.names = 1)
SC<-SC[,c("pos", "EstSC")]
betar<-read.csv("Output/BetaReg/BetaReg.Data.Shape.csv", stringsAsFactors = F, row.names = 1)
betaDF<-merge(betar, SC, by="pos")

################
#2. run Beta Regression
mod.g1 <- betareg(EstSC ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1+NS2+NS4A+NS4B+NS5A+NS5B+Shape, data = betaDF[betaDF$Stop == 0,])
AIC(mod.g1) #-91173.68
summary(mod.g1)

#Coefficients (mean model with logit link):
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -6.545147   0.027430 -238.611  < 2e-16 ***
#t           -0.238879   0.034023   -7.021 2.20e-12 ***
#c           -0.508807   0.030824  -16.507  < 2e-16 ***
#g           -0.580692   0.033799  -17.181  < 2e-16 ***
#CpG          0.067652   0.011731    5.767 8.07e-09 ***
#Nonsyn       0.729289   0.028034   26.014  < 2e-16 ***
#bigAAChange  0.129948   0.010500   12.376  < 2e-16 ***
#Core         0.067545   0.018899    3.574 0.000352 ***
#E1          -0.085109   0.018466   -4.609 4.05e-06 ***
#HVR1        -0.455809   0.059808   -7.621 2.51e-14 ***
#E2          -0.101633   0.015102   -6.730 1.70e-11 ***
#NS1         -0.165872   0.031266   -5.305 1.13e-07 ***
#NS2         -0.133262   0.018183   -7.329 2.32e-13 ***
#NS4A         0.006899   0.031056    0.222 0.824208    
#NS4B         0.006849   0.016215    0.422 0.672766    
#NS5A        -0.050702   0.013623   -3.722 0.000198 ***
#NS5B         0.059375   0.014535    4.085 4.41e-05 ***
#Shape       -0.010908   0.013235   -0.824 0.409827    
#t:Nonsyn     0.167953   0.036103    4.652 3.29e-06 ***
#c:Nonsyn     0.171423   0.033075    5.183 2.19e-07 ***
#g:Nonsyn    -0.089150   0.035996   -2.477 0.013261 *  
    

#remove:    
mod.g2 <- update(mod.g1, ~. -NS4A - NS4B -Shape)
AIC(mod.g2) # 91178.69
summary(mod.g2)

#remove:    
mod.g3 <- update(mod.g1, ~. -NS4A  -Shape)
AIC(mod.g3) 

mod.g3 <- update(mod.g1, ~. -NS4B  -Shape)
AIC(mod.g3)

#remove the ns genes: NS3, NS4A, and MS5A
modI <- update(mod.g2, ~.  -g:Nonsyn  -t:Nonsyn)
AIC(modI) 

AIC(mod.g1,mod.g2,mod.g3)
#       df       AIC
#mod.g1 22 -91173.68
#mod.g2 19 -91178.69
#mod.g3 20 -91176.73

sum.g2<-summary(mod.g2)

modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output/BetaReg/BetaReg_SC_best_model.csv")


