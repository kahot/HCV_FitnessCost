#GLM/Beta regression preparation:
library(tidyverse)
library(zoo)
library(purrr)
library(MASS)
library(betareg)
library(miscTools)
source("Rscripts/baseRscript.R")

####################
SC<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F,row.names = 1)
SC<-SC[,c("pos", "EstSC")]
betar<-read.csv("Output1A/GLM/HCV.Ts.csv", stringsAsFactors = F, row.names = 1)
betaDF<-merge(betar, SC, by="pos")

################
#2. run Beta Regression
mod.g1 <- betareg(EstSC ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1+NS2+NS4A+NS4B+NS5A+NS5B+Shape, data = betaDF[betaDF$Stop == 0,])
AIC(mod.g1) #-91229.13
summary(mod.g1)

#Coefficients (mean model with logit link):
#Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -6.550994   0.027454 -238.621  < 2e-16 ***
#t           -0.235585   0.034034   -6.922 4.45e-12 ***
#c           -0.506702   0.030845  -16.427  < 2e-16 ***
#g           -0.577419   0.033811  -17.078  < 2e-16 ***
#CpG          0.067643   0.011734    5.765 8.17e-09 ***
#Nonsyn       0.729123   0.028060   25.985  < 2e-16 ***
#bigAAChange  0.129986   0.010501   12.378  < 2e-16 ***
#Core         0.067515   0.018900    3.572 0.000354 ***
#E1          -0.085084   0.018466   -4.608 4.07e-06 ***
#HVR1        -0.455732   0.059815   -7.619 2.56e-14 ***
#E2          -0.101629   0.015104   -6.729 1.71e-11 ***
#NS1         -0.165861   0.031262   -5.305 1.12e-07 ***
#NS2         -0.133246   0.018183   -7.328 2.34e-13 ***
#NS4A         0.006985   0.031053    0.225 0.822039    
#NS4B         0.006870   0.016216    0.424 0.671796    
#NS5A        -0.050651   0.013624   -3.718 0.000201 ***
#NS5B         0.059417   0.014536    4.088 4.36e-05 ***
#Shape       -0.010881   0.013235   -0.822 0.410986    
#t:Nonsyn     0.168183   0.036114    4.657 3.21e-06 ***
#c:Nonsyn     0.171550   0.033097    5.183 2.18e-07 ***
#g:Nonsyn    -0.088894   0.036008   -2.469 0.013559 * 

#remove:    
mod.g2 <- update(mod.g1, ~. -NS4A - NS4B -Shape)
AIC(mod.g2) # -91234.15
summary(mod.g2)
#            Estimate Std. Error  z value Pr(>|z|)
#(Intercept) -6.54972    0.02697 -242.851  < 2e-16 ***
#t           -0.23580    0.03403   -6.928 4.26e-12 ***
#c           -0.50666    0.03085  -16.425  < 2e-16 ***
#g           -0.57724    0.03381  -17.072  < 2e-16 ***
#CpG          0.06773    0.01173    5.773 7.80e-09 ***
#Nonsyn       0.72923    0.02806   25.988  < 2e-16 ***
#bigAAChange  0.12979    0.01049   12.367  < 2e-16 ***
#Core         0.05933    0.01694    3.503  0.00046 ***
#E1          -0.08776    0.01779   -4.933 8.12e-07 ***
#HVR1        -0.45704    0.05959   -7.670 1.72e-14 ***
#E2          -0.10398    0.01428   -7.284 3.24e-13 ***
#NS1         -0.16719    0.03082   -5.424 5.82e-08 ***
#NS2         -0.13456    0.01742   -7.725 1.12e-14 ***
#NS5A        -0.05346    0.01271   -4.206 2.60e-05 ***
#NS5B         0.05581    0.01363    4.094 4.24e-05 ***
#t:Nonsyn     0.16851    0.03611    4.666 3.07e-06 ***
#c:Nonsyn     0.17162    0.03310    5.185 2.16e-07 ***
#g:Nonsyn    -0.08908    0.03601   -2.474  0.01337 *  


#remove:    
mod.g3 <- update(mod.g1, ~. -NS4A  -Shape)
AIC(mod.g3) #-91232.37

mod.g3 <- update(mod.g1, ~. -NS4B  -Shape)
AIC(mod.g3) #-91232.19

#remove the ns genes: NS3, NS4A, and MS5A
modI <- update(mod.g2, ~.  -g:Nonsyn  -t:Nonsyn)
AIC(modI) #-91177.53

AIC(mod.g1,mod.g2,mod.g3)
#       df       AIC
#mod.g1 22 -91229.13
#mod.g2 19 -91234.15
#mod.g3 20 -91232.19

#sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)

modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output1A/GLM/BetaReg_SC_best_model.csv")



### plot

set.seed(123)
plot(mod.g2, which=1:4, type="pearson")
plot(mod.g2, which=5, type="deviance",sub.caption="")
plot(mod.g2, which=1, type="deviance",sub.caption="")
