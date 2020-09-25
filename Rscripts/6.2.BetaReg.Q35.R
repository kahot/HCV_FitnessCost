#GLM/Beta regression preparation:
library(tidyverse)
library(zoo)
library(purrr)
library(MASS)
library(betareg)
library(miscTools)
source("Rscripts/baseRscript.R")

####################

#Read data file
betar<-read.csv("Output1A/BetaReg/BetaReg.Data.Shape.csv", stringsAsFactors = F, row.names = 1)

#Run Beta Regression
mod.g1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2+NS3+NS4A+NS5A+NS5B+Shape, data = betar[betar$Stop == 0,])
AIC(mod.g1) #-75866.15
summary(mod.g1)
#using the new data = same results!
#              Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.4737177  0.0229838 -194.647  < 2e-16 ***
#t1           0.1157091  0.0208632    5.546 2.92e-08 ***
#c1          -0.5323171  0.0205090  -25.955  < 2e-16 ***
#g1          -0.7220995  0.0229500  -31.464  < 2e-16 ***
#CpG         -0.0725060  0.0134435   -5.393 6.91e-08 ***
#Nonsyn1     -0.7585679  0.0223505  -33.940  < 2e-16 ***
#bigAAChange -0.1249080  0.0140548   -8.887  < 2e-16 ***
#Core        -0.2152590  0.0240976   -8.933  < 2e-16 ***
#E1           0.0390887  0.0227549    1.718   0.0858 .  
#HVR1         0.2373264  0.0482918    4.914 8.90e-07 ***
#E2           0.0859036  0.0197335    4.353 1.34e-05 ***
#NS1          0.0714742  0.0330309    2.164   0.0305 *  
#NS2          0.0911630  0.0218428    4.174 3.00e-05 ***
#NS3          0.0035818  0.0179292    0.200   0.8417    
#NS4A        -0.0542810  0.0372726   -1.456   0.1453    
#NS5A        -0.0003067  0.0190949   -0.016   0.9872    
#NS5B        -0.1615474  0.0207384   -7.790 6.71e-15 ***
#t1:Nonsyn1  -0.1860138  0.0268162   -6.937 4.02e-12 ***
#c1:Nonsyn1  -0.1431568  0.0277536   -5.158 2.49e-07 ***
#g1:Nonsyn1   0.0929057  0.0286821    3.239   0.0012 ** 


#remove the ns genes: NS3, and NS5A        
mod.g2 <- update(mod.g1, ~. -NS3 - NS5A -Shape)
AIC(mod.g2) # -75870.06  *** BEST MODEL
summary(mod.g2)
#(Intercept) -4.47210    0.01836 -243.547  < 2e-16 ***
#t            0.11574    0.02086    5.548 2.89e-08 ***
#c           -0.53233    0.02050  -25.962  < 2e-16 ***
#g           -0.72216    0.02295  -31.467  < 2e-16 ***
#CpG         -0.07246    0.01343   -5.395 6.86e-08 ***
#Nonsyn      -0.75852    0.02235  -33.941  < 2e-16 ***
#bigAAChange -0.12491    0.01404   -8.894  < 2e-16 ***
#Core        -0.21688    0.01979  -10.960  < 2e-16 ***
#E1           0.03746    0.01816    2.063  0.03916 *  
#HVR1         0.23570    0.04633    5.088 3.63e-07 ***
#E2           0.08428    0.01419    5.940 2.85e-09 ***
#NS1          0.06986    0.03006    2.324  0.02013 *  
#NS2          0.08954    0.01700    5.266 1.40e-07 ***
#NS4A        -0.05591    0.03467   -1.612  0.10689    
#NS5B        -0.16317    0.01556  -10.484  < 2e-16 ***
#t:Nonsyn    -0.18612    0.02681   -6.942 3.87e-12 ***
#c:Nonsyn    -0.14322    0.02775   -5.161 2.46e-07 ***
#g:Nonsyn     0.09295    0.02868    3.241  0.00119 ** 


#remove the ns genes: NS3, NS4A, and MS5A
mod2 <- update(mod.g1, ~. -NS3 -NS4A - NS5A)
AIC(mod2) #-75869.38
summary(mod2)
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.47368    0.01834 -243.971  < 2e-16 ***
#t            0.11525    0.02086    5.524 3.32e-08 ***
#c           -0.53291    0.02051  -25.987  < 2e-16 ***
#g           -0.72290    0.02295  -31.498  < 2e-16 ***
#CpG         -0.07259    0.01343   -5.403 6.54e-08 ***
#Nonsyn      -0.75895    0.02235  -33.958  < 2e-16 ***
#bigAAChange -0.12471    0.01404   -8.879  < 2e-16 ***
#Core        -0.21476    0.01975  -10.873  < 2e-16 ***
#E1           0.03963    0.01812    2.187  0.02878 *  
#HVR1         0.23781    0.04632    5.134 2.83e-07 ***
#E2           0.08641    0.01414    6.113 9.78e-10 ***
#NS1          0.07202    0.03004    2.398  0.01649 *  
#NS2          0.09168    0.01696    5.405 6.47e-08 ***
#NS5B        -0.16105    0.01552  -10.380  < 2e-16 ***
#t:Nonsyn    -0.18616    0.02682   -6.942 3.86e-12 ***
#c:Nonsyn    -0.14228    0.02775   -5.127 2.94e-07 ***
#g:Nonsyn     0.09324    0.02868    3.251  0.00115 ** 

#Phi coefficients (precision model with identity link):
#        Estimate Std. Error z value Pr(>|z|)    
#(phi)  1145.80      18.81   60.92   <2e-16 ***
#
#Type of estimator: ML (maximum likelihood)
#Log-likelihood: 3.795e+04 on 19 Df
#Pseudo R-squared: 0.6125
#Number of iterations: 34 (BFGS) + 3 (Fisher scoring)         

AIC(mod.g1,mod.g2,mod2)
#         df       AIC
#mod.g1   21 -75866.15
#mod.g2   19 -75870.06
#mod.g2.2 18 -75869.38

sum.g1<-summary(mod.g1)
sum.g2<-summary(mod.g2)

modcoef1<-sum.g1$coefficients
coef1<-modcoef1[[1]]
write.csv(coef1,"Output1A/GLM/BetaReg_mod.g1_Ts.Q35.csv")
modcoef2<-sum.g2$coefficients
coef2<-modcoef2[[1]]
write.csv(coef2,"Output1A/BetaReg/BetaReg_mod.g2_Ts.Q35.csv")


### Calculate the effects:

model<-read.csv("Output1A/BetaReg/BetaReg_mod.g2_Ts.Q35.csv", stringsAsFactors = F, row.names = 1)
model$Effect<-''
for (g in 1:length(row.names(model)) ){
        if (g==1){
                model$Effect[1]<- exp(model[1,g])
        }
        else{
                model$Effect[g]<- (((exp(model[1,1] + model$Estimate[g]) - exp(model[1,1])) /exp(model[1,1])))#add estimate % column
        }
}

print(model)
write.csv(model, "Output1A/BetaReg/BetaReg_BestModel.csv")
