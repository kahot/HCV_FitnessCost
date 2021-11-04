# Run Beta regression
library(betareg)

#Read data file
betar<-read.csv("Output/BetaReg/BetaRegData.csv", stringsAsFactors = F, row.names = 1)

#Run Beta Regression
mod1 <- betareg(mean ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                        Core +E1 +HVR1++E2 +NS1 +NS2+NS3+NS4A+NS5A+NS5B+Shape, data = betar[betar$Stop == 0,])
AIC(mod1) #-75864.86
summary(mod1)
#using the new data = same results!
#              Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -4.473459   0.022983 -194.643  < 2e-16 ***
#t            0.115463   0.020862    5.535 3.12e-08 ***
#c           -0.532611   0.020508  -25.970  < 2e-16 ***
#g           -0.722394   0.022949  -31.478  < 2e-16 ***
#CpG         -0.072600   0.013443   -5.401 6.64e-08 ***
#Nonsyn      -0.758651   0.022349  -33.945  < 2e-16 ***
#bigAAChange -0.124901   0.014054   -8.887  < 2e-16 ***
#Core        -0.223432   0.025866   -8.638  < 2e-16 ***
#E1           0.037405   0.022843    1.637  0.10154    
#HVR1         0.237318   0.048289    4.914 8.90e-07 ***
#E2           0.084702   0.019787    4.281 1.86e-05 ***
#NS1          0.071473   0.033029    2.164  0.03047 *  
#NS2          0.091152   0.021842    4.173 3.00e-05 ***
#NS3          0.002000   0.018022    0.111  0.91164    
#NS4A        -0.054267   0.037271   -1.456  0.14538    
#NS5A        -0.002029   0.019209   -0.106  0.91590    
#NS5B        -0.164191   0.020978   -7.827 5.00e-15 ***
#Shape        0.012517   0.014699    0.852  0.39444    
#t:Nonsyn    -0.185900   0.026814   -6.933 4.12e-12 ***
#c:Nonsyn    -0.143042   0.027752   -5.154 2.55e-07 ***
#g:Nonsyn     0.092746   0.028682    3.234  0.00122 **

#remove the ns genes: NS3, and NS5A        
mod2 <- update(mod1, ~. -NS3 - NS5A -Shape)
AIC(mod2) # -75870.06  *** BEST MODEL
summary(mod2)
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


#remove NS4A
mod3 <- update(mod2, ~.-NS4A)
AIC(mod3) #-75869.38
summary(mod3)
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


AIC(mod1,mod2,mod3)
#         df       AIC
#mod.g1   22 -75864.86
#mod.g2   19 -75870.06
#mod.g2.2 18 -75869.38

#Save the best model
sum2<-summary(mod2)
modcoef2<-sum2$coefficients
model<-modcoef2[[1]]
write.csv(model,"Output/BetaReg/BetaReg_MF_bestModel.csv")


### Calculate the effects:
source("Rscripts/BetaEffectSize.R")

model<-read.csv("Output/BetaReg/BetaReg_MF_bestModel.csv", stringsAsFactors = F, row.names = 1)
effect<-BetaEffectSizeDF(model)
effect$factor<-rownames(effect)
effect$factor[effect$factor=='t']<-"T"
effect$factor[effect$factor=='c']<-"C"
effect$factor[effect$factor=='g']<-"G"
effect$factor[effect$factor=='t:Nonsyn']<-"T:Nonsyn"
effect$factor[effect$factor=='c:Nonsyn']<-"C:Nonsyn"
effect$factor[effect$factor=='g:Nonsyn']<-"G:Nonsyn"
effect$factor[effect$factor=='bigAAChange']<-"AAChange"

write.csv(effect, "Output/BetaReg/BetaReg_MF_effectsize.csv", row.names = F)

