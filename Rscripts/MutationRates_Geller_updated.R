#calculate mutation frequency based on the nucleotide from Geller's paper:
require(reshape2)
library(colorspace)
library(ggplot2)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")
colors2<-qualitative_hcl(6, palette="Dark3")


geller<-read.csv("Data/Geller2016dataset.csv",na.strings=c("","NA"),stringsAsFactors=FALSE)
#remove N sites
geller<-geller[geller$Sequence!="N",]

mean<-mean(geller$Site_average_mutation_frequency, na.rm = T) # 0.00002435308 
formatC(mean/23.3, format = "e", digits = 2) #1.02e-05


#Calculate mutation frequencies for each nucleotide subsitiution type for each line

mutationrates<-list()
mutfreq<-list()
for (i in c("Line1","Line2", "Line3")){
        id<-i
        name<-paste0(id,"data")
        DF<-geller[geller[paste0(id,"_valid_case")]==1,]
        
        dt<-DF[,c("H77_site","Sequence","codon_site")]
        dt$depth<-DF[, (paste0(id,"_sequence_depth"))]
        dt$to_U<-(DF[, (paste0(id,"_substitutions_for_U"))]/DF[, paste0(id,"_sequence_depth")])
        dt$to_A<-(DF[, (paste0(id,"_substitutions_for_A"))]/DF[, paste0(id,"_sequence_depth")])
        dt$to_C<-(DF[, (paste0(id,"_substitutions_for_C"))]/DF[, paste0(id,"_sequence_depth")])
        dt$to_G<-(DF[, (paste0(id,"_substitutions_for_G"))]/DF[, paste0(id,"_sequence_depth")])
        
        mutfreq[[i]]<-dt
        mutationrates[[i]]<-aggregate(dt[,4:7], by=list(dt$Sequence), mean, na.action = na.omit)
}


# the order is A, C, G, U
mut.rates<-lapply(mutationrates,function(x) x=x[,2:5])

library(abind) #abind comines multi-dim. arrays.
rates.matrix<-abind(mut.rates,along=3)
Geller.mut.freq<-apply(rates.matrix,c(1,2), mean)
Geller.mut.rates<-Geller.mut.freq/23.3
Geller.mut.rates<-data.frame(Geller.mut.rates)
Geller.mut.rates$nuc<-c("A","C","G","U")

Geller.mut.rates2<-melt(Geller.mut.rates, "nuc")
Geller.mut.rates2$mutations<-paste0(Geller.mut.rates2$nuc, gsub("to_", '', Geller.mut.rates2$variable))

Geller.mut.rates2<-Geller.mut.rates2[,c(4,3)]
colnames(Geller.mut.rates2)<-c("mutation","mut.rate")
write.csv(Geller.mut.rates2,"Data/Geller.mutation.rates_update.csv")



### calculate SE for mutation rates
#muttypes<-Geller.mut.rates2$mutation
MUT<-list()
for (i in 1:3){
        df<-mutfreq[[i]]

        mm<-melt(df, id.vars=c("H77_site","Sequence","depth"), measure.vars=c("to_U","to_A","to_C","to_G"))
        mm$mutation<-paste0(mm$Sequence, gsub("to_", '', mm$variable))
        mm<-mm[,c("mutation",'value',"depth")]
        MUT[[i]]<-mm
        
        
}

MM<-do.call(rbind, MUT)
MM$mut.rate<-MM$value/23.3
mutDF<-aggregate(MM$mut.rate, by=list(MM$mutation),mean)
colnames(mutDF)<-c("Mutation","MutRate")
mutationnames<-mutDF$Mutation

se<-data.frame(Mutation=mutationnames)
for (i in 1:nrow(se)){
        df<-MM[MM$mutation==mutationnames[i],]
        se$se[i]<-sqrt(mean(df$mut.rate)*(1-mean(df$mut.rate))/sum(df$depth))
}

mutDF$SE<-se$se
mutDF$CI<-mutDF$SE*1.96

write.csv(mutDF,"Output/Geller/Geller.MutRates.Summary_updated.csv")

#Plot MutRates with CI
#remove non-mutation combination
mutDF2<-mutDF[!(mutDF$Mutation=="AA"|mutDF$Mutation=="CC"|mutDF$Mutation=="UU"|mutDF$Mutation=="GG"),]

#se1
ggplot(mutDF2, aes(x=Mutation, y=MutRate))+
        geom_errorbar(aes(ymin=pmax(MutRate-SE,0), ymax=MutRate+SE), width=.2, size=.5, color="gray70")+
        geom_point(size =2, color="blue")+
        
        theme(axis.title.x=element_blank())+ylab("Estimated mutation Rate ± SE")+
        theme_bw()+labs(x="")

ggsave("Output/Geller/MutRates_all_SE.pdf", width = 5,height = 3)


ggplot(mutDF, aes(x=Mutation, y=MutRate))+
        geom_errorbar(aes(ymin=pmax(MutRate-CI,0), ymax=MutRate+CI), width=.3, size=.5, color="gray60")+
        geom_point(size =2, color="blue")+
        
        theme(axis.title.x=element_blank())+ylab("Estimated mutation Rate ± 95% CI")+
        theme_bw()+
        labs(x="")
ggsave("Output/Geller/MutRates_all_95CI.pdf", width = 5,height = 3)


#transition only
mutTs<-mutDF[mutDF$Mutation=="AG"|mutDF$Mutation=="GA"|mutDF$Mutation=="CU"|mutDF$Mutation=="UC",]
ggplot(mutTs, aes(x=Mutation, y=MutRate))+
        geom_errorbar(aes(ymin=pmax(MutRate-SE,0), ymax=MutRate+SE), width=.2, size=.5, color="#0000FF66")+
        geom_point(size = 2.5, color="blue")+
        
        theme(axis.title.x=element_blank())+ylab("Estimated mutation rate ± SE")+
        theme_bw()+
        labs(x="")
ggsave("Output/Geller/TsMutRates_SE.pdf", width = 4,height = 3)
p<-list()
ggplot(mutTs, aes(x=Mutation, y=MutRate))+
        geom_errorbar(data=mutTs, aes(ymin=pmax(MutRate-CI,0), ymax=MutRate+CI), width=.2, size=.5, color="#0000FF66")+
        scale_y_continuous(trans = 'log10', labels=label_scientific2)+
        geom_point(size =2.5, color="blue")+
        ylab("Estimated mutation rate ± 95% CI")+
        theme_bw()+theme(panel.grid.major.x = element_blank())+
        scale_x_discrete(labels=c(expression("A" %->% "G"),expression("C" %->% "U"),expression("G" %->% "A"),expression("U" %->% "C")))+ 
        theme(axis.text.x=element_text(color=1, size=12))+xlab('') 
ggsave("Output/Geller/TsMutRates_95CI.pdf", width = 4,height = 3)

MM2<-MM[MM$mutation=="AG"|MM$mutation=="GA"|MM$mutation=="CU"|MM$mutation=="UC",]
ggplot()+
        scale_y_continuous(trans = 'log10', labels=label_scientific2)+
        geom_boxplot(data=MM2, aes(x=mutation, y=mut.rate),outlier.alpha = 0.3, color="gray60",fill=paste0(colors2[5],"66"))+
        geom_errorbar(data=mutTs, aes(x=Mutation, y=mean, ymin=pmax(MutRate-CI,0), ymax=MutRate+CI), width=.2, size=.5, color=colors2[5])+
        geom_point(data=mutTs, aes(x=Mutation, y=MutRate),size =2.5, color="blue")+
        theme(axis.title.x=element_blank())+ylab("Estimated mutation rate ± 95% CI")+
        theme_bw()+
        xlab('')+
        scale_x_discrete(breaks=c("AG","CU","GA","UC"),labels=c(expression(A%->%G),expression(C%->%"T"),expression(G%->%A),expression("T"%->%C)))+
        theme(axis.text.x =element_text(color=1, size=12))
ggsave("Output/Geller/TsMutRates_CI_boxplot.pdf", width =4.5, heigh=4)

library(gridExtra)
pdf(paste0("Output/MutRates.pdf"), width = 10, height = 4)
do.call(grid.arrange, c(p, ncol=2))
dev.off()
ggplot()+
        scale_y_continuous(trans = 'log10', labels=label_scientific2)+
        geom_boxplot(data=MM2, aes(x=mutation, y=mut.rate),outlier.alpha = 0.3, color="gray60",fill=paste0(colors2[5],"66"))+
        
        geom_point(data=mutTs, aes(x=Mutation, y=MutRate),size =2.5, color="blue")+
        geom_errorbar(data=mutTs, aes(x=Mutation, y=MutRate, ymin=pmax(MutRate-SE,0), ymax=MutRate+SE), width=.2, size=.5, color="#0000FF66")+
        ylab("Estimated mutation rate ± SE")+
        theme_bw()+
        xlab('')+
        scale_x_discrete(breaks=c("AG","CU","GA","UC"),labels=c(expression(A%->%G),expression(C%->%"T"),expression(G%->%A),expression("T"%->%C)))+
        theme(axis.text.x =element_text(color=1, size=12))

ggsave("Output/Geller/MutRates_SE.pdf", width =4.5, heigh=4)



#### Wilcoxon Test #####

# AG vs. CT
wilcox.test(MM$value[MM$mutation=="AG"],MM$value[MM$mutation=="CU"], "greater")
#W = 28590000, p-value < 2.2e-16
# AG cs GA
wilcox.test(MM$value[MM$mutation=="AG"],MM$value[MM$mutation=="GA"], "greater")
#W = 28483000, p-value < 2.2e-16
wilcox.test(MM$value[MM$mutation=="UC"],MM$value[MM$mutation=="CU"], "greater")
#W = 30031000, p-value < 2.2e-16
wilcox.test(MM$value[MM$mutation=="UC"],MM$value[MM$mutation=="GA"], "greater")
#W = 3e+07, p-value < 2.2e-16

re1<-wilcox.test(MM$value[MM$mutation=="AG"],MM$value[MM$mutation=="CU"], "greater")

re2<-wilcox.test(MM$value[MM$mutation=="AG"],MM$value[MM$mutation=="GA"], "greater")
re3<-wilcox.test(MM$value[MM$mutation=="UC"],MM$value[MM$mutation=="CU"], "greater")
re4<-wilcox.test(MM$value[MM$mutation=="UC"],MM$value[MM$mutation=="GA"], "greater")
re1[3]
re2[3]
re3[3]
re4[3]

#########################################
#Create Mutation Freqeuncy Table

#divide the main df into 3 lines:

line1<-geller[geller$Line1_valid_case==1,c(1:3,11:20)]
line2<-geller[geller$Line2_valid_case==1,c(1:3,21:30)]
line3<-geller[geller$Line3_valid_case==1,c(1:3,31:40)]
lines<-c("line1","line2","line3")
CSV<-list()
for (j in 1:3 ){
        df<-get(lines[j])
        colnames(df)[5:12]<-c("depth", "u","c","g","a","ins","del","N" )
        df$maj<-(df$depth-df$u-df$c-df$g-df$a-df$N)
        
        for (i in 1:nrow(df)){
             nuc<-tolower(df$Sequence[i])
             df[i,nuc]<-df$maj[i]
        }
        CSV[[j]]<-df
        names(CSV)[j]<-lines[j]
        #write.csv(df,paste0("Output/Geller/",j,".csv"))
}

transition2<-function(nuc){
        if (is.na(nuc)|nuc=="-"){ 
                return (NA)
                next}
        else if (nuc=="a") return("g")
        else if (nuc=="g") return("a")
        else if (nuc=="c") return("u")
        else if (nuc=="u") return("c")
        else if (nuc=='n') return('n')
}


TypeOfSite<-c()
gell<-geller
gell$Sequence<-tolower(gell$Sequence)
gell$Sequence[gell$Sequence=="u"]<-"t"
TypeOfSite<-c()
for (codon in 1:(nrow(gell)/3)) {
        positions <- c(codon*3-2,codon*3-1, codon*3)
        WTcodon <- gell$Sequence[positions]  
        if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
                WTcodon<-c('n','n','n')
                mutant1codon<-c('n','n','n')
                mutant2codon<-c('n','n','n')
                mutant3codon<-c('n','n','n')}
        else{                        
                mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  
                mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
                mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
        }
        
        
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
}
        
gell$Type<-TypeOfSite


MF<-list()
for (i in 1:3){
        df<-CSV[[i]]
        df$Sequence<-tolower(df$Sequence)
        df$transition<-unlist(sapply(df$Sequence, function(x) transition2(x)))
        for (k in 1:nrow(df)){ df$freq.Ts[k]<-df[k,paste0(df$transition[k])]/df$depth[k]}
        write.csv(df,paste0("Output/Geller/Overview_",lines[i],".csv"))
        MF[[i]]<-df
        names(MF)[i]<-names(CSV[i])
}

l1<-MF[[1]]
l1<-l1[,c("H77_site", "freq.Ts")]
l2<-MF[[2]]
l2<-l2[,c("H77_site", "freq.Ts")]
l3<-MF[[3]]
l3<-l3[,c("H77_site", "freq.Ts")]


mut<-merge(l1, l2, by="H77_site", all=T)
mut<-merge(mut, l3, by="H77_site", all=T)
colnames(mut)[2:4]<-c("Line1","Line2","Line3")
mut$mean<-apply(mut[,2:4],1,mean, na.rm=T)
MutFreq<-merge(mut,gell[,c("H77_site","Sequence","Protein","Type")],by="H77_site", all=T)
colnames(MutFreq)<-c("pos","Line1","Line2","Line3","mean","ref","Protein", "Type")

write.csv(MutFreq, "Output/Geller/Geller_MutFrequency.csv")

###
#Syn vs. nonsyn

mf<-MutFreq[MutFreq$pos>381 & MutFreq$pos<=9336,]
nrow(mf[is.na(mf$mean),])

ov<-read.csv("Output/H77_Overview.csv",row.names = 1)
colnames(ov)[2:3]<-c("refH77","typeH77")

mf<-merge(mf, ov, by="pos", all.x = T)

mean(mf$mean[mf$makesCpG==1],na.rm=T) #0.0002717275
mean(mf$mean[mf$makesCpG==0],na.rm=T) #0.0001564748

mean(mf$mean[mf$bigAAChange==1],na.rm=T) #0.0001881563
mean(mf$mean[mf$bigAAChange==0],na.rm=T) #0.0001597525

write.csv(mf,"Output/Geller/Geller_mf_overview.csv")


###
#Run beta regression on in vitro data:
#Format the data:
nucord <- c("a", "t", "c", "g")
dfform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
nuc<-as.matrix(data.frame(a="",t="",c="",g=""))

dat<-mf
for (i in 1:nrow(dat)){
        atcg <- c(0,0,0,0)
        atcg[which(nucord == dat[i,]$ref)] <- 1
        nuc<-miscTools::insertRow(nuc,i,atcg)
        nonsyn <- as.numeric(regexpr("nonsyn",dat[i,]$typeH77) > 0)
        stop <- as.numeric(regexpr("stop",dat[i,]$typeH77) > 0)
        syn<-as.numeric(regexpr("^syn",dat[i,]$typeH77) > 0)
        new<-c(atcg,syn,nonsyn,stop)
        dfform<-miscTools::insertRow(dfform,i,new)
}

BRData<-cbind(mf[,c("pos","makesCpG","bigAAChange","mean")],dfform[1:nrow(dat),])

colnames(BRData)[2]<-"CpG"

#Add genes
genes<-read.csv("Data/HCV_annotations2.csv")
genenames<-genes$Gene[2:12]
gene<-c()
for (i in 2:12){
        gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)
BR2<-BRData
BRData<-merge(BRData, g, by ="pos", all.x=T)

for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(BRData)
        BRData[,n+1]<-0
        colnames(BRData)[n+1]<-gname
        BRData[BRData$gene==i,n+1]<-1
}
co<-which(colnames(BRData)=="NS1(P7)" )
colnames(BRData)[co]<-"NS1"
#Shape data
rna<-read.csv("Data/RNAStructure_Conserved2.csv", stringsAsFactors = F)
rnaList<-list()
for (i in 1:nrow(rna)){
        pos<-rna$Start[i]:rna$End[i]
        Shape<-rep(1, times=rna$End[i]-rna$Start[i]+1)
        rnaList[[i]]<-cbind(pos, Shape)
}

RnaShape<-data.frame(do.call(rbind,rnaList))
BRData<-merge(BRData,RnaShape, by="pos", all.x = T)
BRData$Shape[is.na(BRData$Shape)]<-0

which(colnames(BRData)=="gene")
BRData<-BRData[,-12]
BRData<-BRData[!is.na(BRData$mean),]
write.csv(BRData, "Output/Geller/geller_BRData.csv")
BRData<-read.csv("Output/Geller/geller_BRData.csv", row.names = 1, stringsAsFactors = F)

library(betareg)
df<-BRData[BRData$Stop == 0,]
#https://stackoverflow.com/questions/26385617/proportion-modeling-betareg-errors
#f y also assumes the extremes 0 and 1, a useful transformation in practice is (y · (n − 1) + 0.5)/n where n is the sample size (Smithson and Verkuilen 2006).

y.transf.betareg <- function(y){
        n.obs <- sum(!is.na(y))
        (y * (n.obs - 1) + 0.5) / n.obs
}



m1<-betareg(y.transf.betareg(mean) ~ t + c + g + CpG  + t:Nonsyn + c:Nonsyn + g:Nonsyn + Nonsyn +bigAAChange + 
                    Core +E1 +HVR1++E2 +NS1 +NS2+NS3+NS4A+NS5A+NS5B+Shape, data = BRData[BRData$Stop == 0,])

summary(m1)
AIC(m1) #-130092.8

m2<- betareg(y.transf.betareg(mean) ~ t + c + g +bigAAChange + 
                     Core +E1+NS1 +NS2+NS5B, data = BRData[BRData$Stop == 0,])

summary(m2)
AIC(m2) # -130105.3
#            Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -8.00720    0.01562 -512.500  < 2e-16 ***
#t           -0.04384    0.01774   -2.471 0.013475 *  
#c           -0.59383    0.01846  -32.163  < 2e-16 ***
#g           -0.71589    0.01924  -37.205  < 2e-16 ***
#bigAAChange  0.03719    0.01321    2.815 0.004877 ** 
#Core        -0.06044    0.03035   -1.991 0.046432 *  
#E1          -0.07366    0.02745   -2.684 0.007281 ** 
#NS1         -0.19609    0.04903   -3.999 6.36e-05 ***
#NS2         -0.19376    0.02714   -7.138 9.46e-13 ***
#NS5B        -0.05956    0.01701   -3.502 0.000461 ***


m3<- update(m2, ~. -Core)
m4<-update(m3, ~. -NS5B )
m5<-update(m4, ~. -E1 )
m6<-update(m5, ~. -bigAAChange )

AIC(m2,m3,m4,m5,m6)
#   df       AIC
#m2 11 -130105.3
#m3 10 -130103.2
#m4  9 -130094.5
#m5  8 -130092.1
#m6  7 -130085.7

source("Rscripts/BetaEffectSize.R")

effects<-BetaEffectSize(m2)
write.csv(effects, "Output/Geller/beta_reg_effects.csv")




#Calculate the effect size:
#result=glm or betareg result
BetaEffectSize<-function(result){
        res<-summary(result)
        model<-data.frame(res$coefficients[[1]])
        model$Effect<-''
        for (i in 1:length(row.names(model)) ){
                if (i==1) model$Effect[1]<- exp(model[1,i])

                else{
                        model$Effect[i]<- (((exp(model[1,1] + model$Estimate[i]) - exp(model[1,1])) /exp(model[1,1]))) }
        }
        return(model)
}




#####
counts<-data.frame(nt=c("a","c","g","t"))
nucs<-c("a","c","g","t")
for (i in 1:4){
        counts$syn[i]<-nrow(mf[mf$ref==nucs[i]&mf$Type=="syn",])
        counts$nonsyn[i]<-nrow(mf[mf$ref==nucs[i]&mf$Type=="nonsyn",])
}

table(mf[,c("ref","Type")])
counts
#  nt syn nonsyn
#1  a 630   1158
#2  c 867   1719
#3  g 792   1556
#4  t 719   1188

### Look at C->T mutation
mfC<-mf[mf$ref=="c",]
nrow(mfC[is.na(mfC$mean),]) #27 NA
mfC<-mfC[,c("pos","mean","Type")]
mfC<-mfC[order(mfC$mean),]
mfC<-mfC[!is.na(mfC$mean),]
summary(mfC$mean)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     
#0.000000 0.000029 0.000061 0.0001021 0.000196 0.005135   

mfC$quartile <- with(mfC, cut(mean, 
                                breaks=quantile(mean, probs=seq(0,1, by=0.25), na.rm=TRUE), 
                                include.lowest=TRUE, labels=c("1","2","3","4")))

write.csv(mfC,"Output/Geller/mf_C.csv")
mfC.q<-aggregate(mfC$mean, by=list(mfC$quartile), mean)
colnames(mfC.q)<-c("Q","mf")
mfC.q$mr<-mfC.q$mf/23.3



mfC1<-mfC[mfC$Type==""]


### where are the LM and HM sites?
Geller<-read.csv("Data/Geller2016dataset.csv",na.strings=c("","NA"),stringsAsFactors=FALSE)
plot(Geller$H77_site[Geller$LM_site==1], Geller$Site_average_mutation_frequency[Geller$LM_site==1], ylim=c(0,0.00001), pch=16)
plot(Geller$H77_site[Geller$HM_site==1], Geller$Site_average_mutation_frequency[Geller$HM_site==1], pch=16)
