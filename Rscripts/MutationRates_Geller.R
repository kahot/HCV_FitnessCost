#calculate mutation frequency based on the nucleotide from Geller's paper:
library(reshape2)
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


#Calculate mutation frequencies for each nucleotide substitution type for each line

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
        mutationrates[[i]]<-aggregate(dt[,5:8], by=list(dt$Sequence), mean, na.action = na.omit)
}


MUT<-list()
for (i in 1:3){
        df<-mutfreq[[i]]
        
        mm<-melt(df, id.vars=c("H77_site","Sequence","depth"), measure.vars=c("to_U","to_A","to_C","to_G"))
        mm$mutation<-paste0(mm$Sequence, gsub("to_", '', mm$variable))
        mm<-mm[,c("H77_site","mutation",'value',"depth")]
        mm$Line<-names(mutfreq)[i]
        MUT[[i]]<-mm
}

MM<-do.call(rbind, MUT)
MM$mut.rate<-MM$value/23.3
mutDF<-aggregate(MM$mut.rate, by=list(MM$mutation),mean)
colnames(mutDF)<-c("Mutation","MutRate")
mutationnames<-mutDF$Mutation
write.csv(mutDF,"Output/Geller/Geller.MutRates.Summary_updated.csv")


## estimate CIs for transition mutations

#first insert NA into the invalid sites. 
Depths<-data.frame(pos=DF$H77_site)
Lines<-list()
lns<-c("Line1","Line2", "Line3")
for (i in 1:3){
    id<-lns[i]
    name<-paste0(id,"data")
    line.df<-DF
    line.df<-line.df[,c(1,3, grep(id, colnames(line.df)))]
    line.df[line.df[paste0(id,"_valid_case")]==0,4:8]<-NA
    line.df[line.df[paste0(id,"_valid_case")]==0,]
    
    #create a depth matrix
    Depths[,id]<-line.df[, (paste0(id,"_sequence_depth"))]
    
    #create a count matrix
    line.df<-line.df[,c(1,2,5:8)]
    colnames(line.df)[3:6] <-gsub(paste0(id,"_substitutions_for_"), 'to',colnames(line.df[3:6]))
    Lines[[i]]<-line.df
    names(Lines)[i]<-id
}

meta<-Lines[[1]][1:2]
ToA<-cbind(meta, data.frame(sapply(Lines,"[[","toA")))
ToC<-cbind(meta, data.frame(sapply(Lines,"[[","toC")))
ToG<-cbind(meta, data.frame(sapply(Lines,"[[","toG")))
ToU<-cbind(meta, data.frame(sapply(Lines,"[[","toU")))

ToA$Mutation<-paste0(ToA$Sequence, "A")
ToC$Mutation<-paste0(ToA$Sequence, "C")
ToG$Mutation<-paste0(ToA$Sequence, "G")
ToU$Mutation<-paste0(ToA$Sequence, "U")

#extract transition mutations from the 4 dataframes and combine
muts<-c()
muts<-rbind(muts, ToA[ToA$Mutation=="GA",])
muts<-rbind(muts, ToC[ToC$Mutation=="UC",])
muts<-rbind(muts, ToG[ToG$Mutation=="AG",])
muts<-rbind(muts, ToU[ToU$Mutation=="CU",])
muts<-muts[order(muts$H77_site),]

Depths<-Depths[Depths$pos%in% muts$H77_site,]

#CI estimates using F-distributions
ci_low<-data.frame(pos=muts$H77_site)
ci_up<-data.frame(pos=muts$H77_site)

##calculate F values
#n1= n(total read #) -x (mutations#)+1
muts$n1<-Depths$Line1-muts$Line1+1
muts$n2<-Depths$Line2-muts$Line2+1
muts$n3<-Depths$Line3-muts$Line3+1
muts$w<-apply(muts[3:5],1, function(x) length(x[!is.na(x)]))

p=0.05
muts$f1.1<-apply(muts[,c("Line1","n1")], 1, function(x) ifelse(x["Line1"]==0, 0, qf(p/2, x["n1"]*2, x["Line1"]*2, lower.tail=FALSE)))
muts$f2.1<-qf(p/2, (muts$Line1+1)*2, (Depths$Line1-muts$Line1)*2, lower.tail=FALSE)

muts$f1.2<-apply(muts[,c("Line2","n2")], 1, function(x) ifelse(x["Line2"]==0, 0, qf(p/2, x["n2"]*2, x["Line2"]*2, lower.tail=FALSE)))
muts$f2.2<-qf(p/2, (muts$Line2+1)*2, (Depths$Line2-muts$Line2)*2, lower.tail=FALSE)

muts$f1.3<-apply(muts[,c("Line3","n3")], 1, function(x) ifelse(x["Line3"]==0, 0, qf(p/2, x["n3"]*2, x["Line3"]*2, lower.tail=FALSE)))
muts$f2.3<-qf(p/2, (muts$Line3+1)*2, (Depths$Line3-muts$Line3)*2, lower.tail=FALSE)

ci_low[,"Line1"]<-apply(muts[,c("Line1","n1","f1.1")], 1, function(x) ifelse(x["Line1"]==0, 0,x["Line1"]/(x["Line1"]+x["n1"]*x["f1.1"]))) 
ci_low[,"Line2"]<-apply(muts[,c("Line2","n2","f1.2")], 1, function(x) ifelse(x["Line2"]==0, 0,x["Line2"]/(x["Line2"]+x["n2"]*x["f1.2"]))) 
ci_low[,"Line3"]<-apply(muts[,c("Line3","n3","f1.3")], 1, function(x) ifelse(x["Line3"]==0, 0,x["Line3"]/(x["Line3"]+x["n3"]*x["f1.3"]))) 

ci_up[,"Line1"]<-((muts$Line1+1)*muts$f2.1)/((Depths$Line1-muts$Line1)+(muts$Line1+1)*muts$f2.1)
ci_up[,"Line2"]<-((muts$Line2+1)*muts$f2.2)/((Depths$Line2-muts$Line2)+(muts$Line2+1)*muts$f2.2)
ci_up[,"Line3"]<-((muts$Line3+1)*muts$f2.3)/((Depths$Line3-muts$Line3)+(muts$Line3+1)*muts$f2.3)

#CI esitmation for each position across populations
df.CI<-data.frame(pos=muts$H77_site, Mutation=muts$Mutation)

#rescaling factor for each site. Assing an equal weight for all samples at each position
Depths$wi2<-(1/muts$w)^2
Depths$W<-1/muts$w
r1=apply(Depths, 1, function(x) (sum(1/x[2:4], na.rm=T)*x["wi2"])^0.5)
r2=apply(Depths, 1, function(x) sum(1/sqrt(x[2:4]), na.rm=T)*x["W"])

R=r1/r2

#mut rates at each site for transition mutations
ts<-c("AG","GA","UC","CU")
MM2<-MM[MM$mutation %in% ts,]
MM_sites<-aggregate(MM2$mut.rate, by=list(MM2$H77_site),mean, na.rm=T)

MM_sites<-MM_sites[order(MM_sites$Group.1),]
MM_sites<-MM_sites[MM_sites$Group.1 %in% muts$H77_site,]

#Mean (equal weights estimation) w/o rescaling factors
df.CI[,"lower"]<-rowMeans(ci_low[2:4], na.rm=T)/23.3
df.CI[,"upper"]<-rowMeans(ci_up[2:4], na.rm=T)/23.3
#CI adjusted by the rescaling factor
df.CI[,"lower.R"]<-(MM_sites$x-df.CI$lower)*R
df.CI[,"upper.R"]<-(df.CI$upper-MM_sites$x)*R

#length(MM_sites$Group.1[MM_sites$Group.1==df.CI$pos])
#length(MM_sites$Group.1[MM_sites$Group.1==muts$H77_site])

df.CI$M.rates<-MM_sites$x
df.CI$diff<-MM_sites$x-df.CI$lower

CI.means<-aggregate(df.CI[3:6], by=list(df.CI$Mutation), mean, na.rm=T)
CI.means

CI.means<-merge(mutDF, CI.means, by="Group.1")
ggplot()+
    scale_y_continuous(trans = 'log10', labels=label_scientific2)+
    geom_boxplot(data=MM2, aes(x=mutation, y=mut.rate),outlier.alpha = 0.3, color="gray60",fill=paste0(colors2[5],"66"), width=0.5)+
    geom_errorbar(data=CI.means, aes(x=Group.1, y=x, ymin=x-lower.R, ymax=x+upper.R), width=0.1, size=.3, color="#507FB9")+
    geom_point(data=CI.means, aes(x=Group.1, y=x),size =1.8, color="blue")+
    theme(axis.title.x=element_blank())+ylab("Estimated mutation rate ± 95% CI")+
    theme_classic()+
    xlab('')+
    scale_x_discrete(breaks=c("AG","CU","GA","UC"),labels=c(expression(A%->%G),expression(C%->%"T"),expression(G%->%A),expression("T"%->%C)))+
    theme(axis.text.x =element_text(color=1, size=12))
ggsave("Output/Geller/MutRates_withCIestimates.pdf", width =4.5, heigh=4)

ggplot()+
    scale_y_continuous(trans = 'log10', labels=label_scientific2)+
    geom_boxplot(data=MM2, aes(x=mutation, y=mut.rate),outlier.alpha = 0.3, color="gray60",fill=paste0(colors2[5],"66"), width=0.5)+
    geom_errorbar(data=CI.means, aes(x=Group.1, y=x, ymin=lower, ymax=upper), width=0.1, size=.3, color="#507FB9")+
    geom_point(data=CI.means, aes(x=Group.1, y=x),size =1.8, color="blue")+
    theme(axis.title.x=element_blank())+ylab("Estimated mutation rate ± 95% CI")+
    theme_classic()+
    xlab('')+
    scale_x_discrete(breaks=c("AG","CU","GA","UC"),labels=c(expression(A%->%G),expression(C%->%"T"),expression(G%->%A),expression("T"%->%C)))+
    theme(axis.text.x =element_text(color=1, size=12))
ggsave("Output/Geller/MutRates_withCIestimates_withoutRescaling.pdf", width =4.5, heigh=4)


write.csv(CI.means,"Output/Geller/Geller.MutRates.CIs.csv")


#estimate CI with Agresi-Coull method
DF<-geller[!(geller$Line1_valid_case==0&geller$Line2_valid_case==0&geller$Line3_valid_case==0),]
ci<-data.frame(DF$H77_site)

#first insert NA into the invalid sites. 
Depths<-data.frame(pos=DF$H77_site)
Lines<-list()
lns<-c("Line1","Line2", "Line3")
for (i in 1:3){
    id<-lns[i]
    name<-paste0(id,"data")
    line.df<-DF
    line.df<-line.df[,c(1,3, grep(id, colnames(line.df)))]
    line.df[line.df[paste0(id,"_valid_case")]==0,4:8]<-NA
    line.df[line.df[paste0(id,"_valid_case")]==0,]
    
    #create a depth matrix
    Depths[,id]<-line.df[, (paste0(id,"_sequence_depth"))]
    
    #create a count matrix
    line.df<-line.df[,c(1,2,5:8)]
    colnames(line.df)[3:6] <-gsub(paste0(id,"_substitutions_for_"), 'to',colnames(line.df[3:6]))
    Lines[[i]]<-line.df
    names(Lines)[i]<-id
}

meta<-Lines[[1]][1:2]
ToA<-cbind(meta, data.frame(sapply(Lines,"[[","toA")))
ToC<-cbind(meta, data.frame(sapply(Lines,"[[","toC")))
ToG<-cbind(meta, data.frame(sapply(Lines,"[[","toG")))
ToU<-cbind(meta, data.frame(sapply(Lines,"[[","toU")))

ToA$Mutation<-paste0(ToA$Sequence, "A")
ToC$Mutation<-paste0(ToA$Sequence, "C")
ToG$Mutation<-paste0(ToA$Sequence, "G")
ToU$Mutation<-paste0(ToA$Sequence, "U")

muts<-list()
muts[[1]]<-ToA
muts[[2]]<-ToC
muts[[3]]<-ToG
muts[[4]]<-ToU

CIs_nuc<-list()
for (j in 1:4){
    df<-muts[[j]]
    ci_low<-data.frame(df$H77_site)
    ci_up<-data.frame(df$H77_site)
    for (i in 1: nrow(df)){
        if  (is.na(df$Line1[i])) ci_low$lower1[i]=NA; ci_up$upper1[i]=NA
        if (!is.na(df$Line1[i])) {
            #if (df$Line1[i]==0) cis$lower1[i]<-0; cis$lower1[i]<-0
            ci_estimates<-binom.confint(df$Line1[i], Depths$Line1[i], conf.level=0.95, methods="agresti-coull")
            ci_low$lower1[i]<-ci_estimates$lower
            ci_up$upper1[i]<-ci_estimates$upper
        }
        if  (is.na(df$Line2[i])) ci_low$lower2[i]=NA; ci_up$upper2[i]=NA
        if (!is.na(df$Line2[i])) {
            ci_estimates<-binom.confint(df$Line2[i], Depths$Line2[i], conf.level=0.95, methods="agresti-coull")
            ci_low$lower2[i]<-ci_estimates$lower
            ci_up$upper2[i]<-ci_estimates$upper
        }
        if  (is.na(df$Line3[i])) ci_low$lower3[i]=NA; ci_up$upper3[i]=NA
        if (!is.na(df$Line3[i])) {
            ci_estimates<-binom.confint(df$Line3[i], Depths$Line3[i], conf.level=0.95, methods="agresti-coull")
            ci_low$lower3[i]<-ci_estimates$lower
            ci_up$upper3[i]<-ci_estimates$upper
        }
    }
    cis<-data.frame(df$H77_site)
    cis$lower<-rowMeans(ci_low[2:4], na.rm=T)
    cis$upper<-rowMeans(ci_up[2:4], na.rm=T)
    cis$Mutations<-df$Mutation
    CIs_nuc[[j]]<-cis    
}

CI.mean<-data.frame()
for ( i in 1:4){
    df<-CIs_nuc[[i]]
    ave<-aggregate(df[,2:3], by=list(df$Mutations), mean, na.rm=T)  
    CI.mean<-rbind(CI.mean, ave)
}
ts<-c("AG","GA","UC","CU")
CI.Ts<-CI.mean[CI.mean$Group.1 %in% ts,]
colnames(CI.Ts)[1]<-"Mutation"
CI.Ts$lower<-CI.Ts$lower/23.3
CI.Ts$upper<-CI.Ts$upper/23.3
CI.Ts$lower[CI.Ts$lower<=0]<-0.00000000001
mutTs2_new<-merge(CI.Ts, mutTs2[,1:2], by="Mutation")


ggplot()+
    scale_y_continuous(trans = 'log10', labels=label_scientific2)+
    geom_boxplot(data=MM2, aes(x=mutation, y=mut.rate),outlier.alpha = 0.3, color="gray60",fill=paste0(colors2[5],"66"), width=0.5)+
    geom_errorbar(data=mutTs2_new, aes(x=Mutation, y=MutRate, ymin=lower, ymax=upper), width=0.1, size=.3, color="gray30")+
    geom_point(data=mutTs2, aes(x=Mutation, y=MutRate),size =1.8, color="blue")+
    theme(axis.title.x=element_blank())+ylab("Estimated mutation rate ± 95% CI")+
    theme_classic()+
    xlab('')+
    scale_x_discrete(breaks=c("AG","CU","GA","UC"),labels=c(expression(A%->%G),expression(C%->%"T"),expression(G%->%A),expression("T"%->%C)))+
    theme(axis.text.x =element_text(color=1, size=12))
ggsave("Output/Geller/TsMutRates_CIac.pdf", width =4.5, heigh=4)




    


#### Wilcoxon Test of mutation rates #####
mut<-c("AG","CU","GA","UC")
mcomb<-t(combn(mut,2))

WilcoxTest.mut<-data.frame(matrix(ncol=4,nrow=nrow(mcomb)))
colnames(WilcoxTest.mut)<-c("Mut1","Mut2","test","rawP")
for (i in 1:nrow(mcomb)) {
        vec1<-MM$value[MM$mutation==mcomb[i,1]]
        vec2<-MM$value[MM$mutation==mcomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
        
        WilcoxTest.mut$Mut1[i]<-mcomb[i,1]
        WilcoxTest.mut$Mut2[i]<-mcomb[i,2]
        WilcoxTest.mut$test[i]<-"greater"
        WilcoxTest.mut$rawP[i]<-result[[3]]
}   

WilcoxTest.mut2<-data.frame(matrix(ncol=4,nrow=nrow(mcomb)))
colnames(WilcoxTest.mut2)<-c("Mut1","Mut2","test","rawP")

for (i in 1:nrow(Ncomb)) {
        vec1<-MM$value[MM$mutation==mcomb[i,1]]
        vec2<-MM$value[MM$mutation==mcomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest.mut2$Mut1[i]<-mcomb[i,1]
        WilcoxTest.mut2$Mut2[i]<-mcomb[i,2]
        WilcoxTest.mut2$test[i]<-"less"
        WilcoxTest.mut2$rawP[i]<-result[[3]]
}   

WilcoxTestMut<-rbind(WilcoxTest.mut,WilcoxTest.mut2)
WilcoxTestMut<-Pcorrection(WilcoxTestMut)
write.csv(WilcoxTestMut, "Output/Geller/MutationRates_test.csv")


#########################################
#Create Mutation Frequency Overview

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

dat1<-mf[,c("mean","ref","Type")]

library(caret)
new<-dummyVars("~.", data=dat1)
brDF<-data.frame(predict(new, newdata = dat1))
colnames(brDF)[2:7]<-c("a","c","g","t","Nonsyn","Stop")
#genenm<-gsub("Protein","", colnames(brDF)[9:18])
#colnames(brDF)[9:18]<-genenm

brDF<-cbind(brDF[,1:7], mf[,c("pos","makesCpG", "bigAAChange")])
colnames(brDF)[which(colnames(brDF)=="makesCpG")]<-"CpG"

#Add genes
genes<-read.csv("Data/HCV_annotations2.csv")
gene<-c()
for (i in 2:12){
        gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)
BRData<-merge(brDF, g, by ="pos", all.x=T)

for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(BRData)
        BRData[,n+1]<-0
        colnames(BRData)[n+1]<-gname
        BRData[BRData$gene==i,n+1]<-1
}
colnames(BRData)[which(colnames(BRData)=="NS1(P7)")]<-"NS1"
#Shape data
rna<-read.csv("Data/RNAStructure_Conserved.csv", stringsAsFactors = F)
rnaList<-list()
for (i in 1:nrow(rna)){
        pos<-rna$Start[i]:rna$End[i]
        Shape<-rep(1, times=rna$End[i]-rna$Start[i]+1)
        rnaList[[i]]<-cbind(pos, Shape)
}

RnaShape<-data.frame(do.call(rbind,rnaList))
BRData<-merge(BRData,RnaShape, by="pos", all.x = T)
BRData$Shape[is.na(BRData$Shape)]<-0
BRData<-BRData[!is.na(BRData$mean),]
#write.csv(BRData, "Output/Geller/geller_BRData.csv")
#BRData<-read.csv("Output/Geller/geller_BRData.csv", row.names = 1, stringsAsFactors = F)

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
AIC(m1) #-129176.4

#Remove non-significant factors
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
#m2 11 -129184.2
#m3 10 -129182.4
#m4  9 -129174.5
#m5  8 -129172.6
#m6  7 -129167.0

source("Rscripts/BetaEffectSize.R")

effects<-BetaEffectSize(m2)
write.csv(effects, "Output/Geller/beta_reg_effects.csv")





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
