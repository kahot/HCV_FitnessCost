# Calculate the effect sizes from beta regression results
# 'result' = betareg result (the output of betareg(), e.g. mod1)
BetaEffectSize<-function(result){
        res<-summary(result)
        model<-data.frame(res$coefficients[[1]])
        model$EstValue<-''
        model$Effect<-''
        for (i in 1:length(row.names(model)) ){
                if (i==1) {model$EstValue[1]<- exp(model[1,i])
                model$Effect[1]<-0}
                
                else{   model$EstValue[i]<-exp(model[1,1] + model$Estimate[i])
                model$Effect[i]<- (((exp(model[1,1] + model$Estimate[i]) - exp(model[1,1])) /exp(model[1,1]))) }
        }
        
        model$factor<-rownames(model)
        factors<-model$factor
        if ("t" %in% factors) model$factor[model$factor=='t']<-"T"
        if ("c" %in% factors) model$factor[model$factor=='c']<-"C"
        if ("g" %in% factors) model$factor[model$factor=='g']<-"G"
        if ("t:Nonsyn" %in% factors) model$factor[model$factor=='t:Nonsyn']<-"T:Nonsyn"
        if ("c:Nonsyn" %in% factors) model$factor[model$factor=='c:Nonsyn']<-"C:Nonsyn"
        if ("g:Nonsyn" %in% factors) model$factor[model$factor=='g:Nonsyn']<-"G:Nonsyn"
        if ("bigAAChange" %in% factors) model$factor[model$factor=='bigAAChange']<-"AAChange"
        
        return(model)
}

