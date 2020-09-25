library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)
library(cowplot)
scCols<-c("#E16A86","#009ADE")
source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")
###########

#### Site by Site mut rate:
mf<-read.csv("Output1A/Geller/Geller_mf_overview.csv", row.names = 1, stringsAsFactors = F)
min(mf$mean[mf$mean>0], na.rm = T) #2.780489e-06
nrow(mf[mf$mean==0,]) #1044

#replace 0 with 1*10^-6
mf2<-mf
mf2$mean[mf2$mean==0]<-1*10^-6

mf2$mr<-mf2$mean/23.3

df<-read.csv("Output1A/Mutfreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names = 1 )
df<-df[,c(1,197:203)]

df.ind<-merge(df, mf2[,c("pos","mr")], by="pos")
df.ind$EstSC<-apply(df.ind[c("mr","mean")],1, function(x) EstimatedS(x[1],x[2]))
df.ind$EstSC<-as.numeric(df.ind$EstSC)
df<-df.ind
mean(df$EstSC, na.rm = T) # 0.002026144 vs. 0.002036205 (single mutrate for each nt)

depth<-read.csv("Output1A/ReadDepth_sum.csv", stringsAsFactors = F, row.names = 1)
df<-merge(df,depth, by="pos", all.x=T)
#write.csv(df,"Output1A/SelCoeff/SC.csv") #add the rad depth to SC files
df<-df[!is.na(df$EstSC),]


#### 
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-df$pos[nrow(df)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]

sc<-merge(df, genetable, by="pos")


#### Histogram comparison
SC<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
colnames(SC)[which(colnames(SC)=="EstSC")]<-"SC_original"
sc<-merge(sc, SC[,c("pos","SC_original")], by="pos")


scDF<-sc[,c("SC_original","EstSC")]
range(scDF$SC_original) #0.0001416349 0.0101180076
range(scDF$EstSC) # 1.63639e-06 2.19968e-01

colnames(scDF)<-c("Uniform", "Site-level")
scD<-melt(scDF)
colnames(scD)[1]<-c("SC")
library(gridExtra)

p1<-ggplot(scD[scD$SC=="Uniform",], aes(x=value))+ggtitle("Single mut rate")+geom_histogram(color="gray60", alpha = 0.5)+theme_bw()+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.000001,.5))+xlab("selection coefficient")
p2<-ggplot(scD[scD$SC=="Site-level",], aes(x=value))+ggtitle("Site-level mut rate")+geom_histogram(color="gray60", alpha = 0.5)+theme_bw()+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.000001,.5))+xlab("selection coefficient")


pdf("Output1A/SelCoeff_indivMutRates/SC_comparison_all.pdf", width = 4, height = 4)
grid.arrange(p1,p2,nrow=2)
dev.off()





##########################################
#Plot summary of 1) sel coef by type by nucleotide       
#1. 
k=1
transSC<-list()
for (i in c("a","t","c","g")) {
        for (type in c("syn","nonsyn")){
                datavector<-sc$EstSC[sc$Type==type & sc$ref==i]
                nt<-toupper(i)
                vname<-paste0(nt,".",type)
                dat<-data.frame(base=rep(nt,times=length(datavector)),
                                type=rep(type, times=length(datavector)), S.C.=datavector)
                transSC[[k]]<-dat
                names(transSC)[k]<-vname
                k=k+1
        }
}

scdata<-do.call(rbind, transSC)

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 
        
ggplot(scdata,aes(x=base, y=S.C., fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
        scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.0001,0.001,0.01,0.1), minor_breaks=ybreaks, labels=label_scientific2)+
        labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_fill_manual(values=col2_light[c(5,1)], ) + theme_bw()+
        theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")
ggsave("Output1A/SelCoeff/SC.byNT.indivMutRates.pdf", width = 5,height = 4)

summary<-data.frame(nt=c("a","c","g","t"))
summary$mean.inv[1]<-mean(sc$EstSC[sc$ref=="a"])
summary$mean.inv[2]<-mean(sc$EstSC[sc$ref=="c"])
summary$mean.inv[3]<-mean(sc$EstSC[sc$ref=="g"])
summary$mean.inv[4]<-mean(sc$EstSC[sc$ref=="t"])
#[1] 0.003035399
#[1] 0.001652175
#[1] 0.001381607
#[1] 0.0025262

summary$mean.sin[1]<-mean(SC$SC_original[SC$ref=="a"],na.rm=T)
summary$mean.sin[2]<-mean(SC$SC_original[SC$ref=="c"],na.rm=T)
summary$mean.sin[3]<-mean(SC$SC_original[SC$ref=="g"],na.rm=T)
summary$mean.sin[4]<-mean(SC$SC_original[SC$ref=="t"],na.rm=T)

summaryM<-melt(summary)

ggplot(summaryM, aes(x=nt, y=value, color=variable))+
        scale_color_manual(labels=c("Site-level mut rate","NT-level mut rate"), values=scCols)+
        geom_point(size=2)+theme_bw()+
        ylab("Mean selection coefficient")+xlab("")+
        scale_x_discrete(breaks=c("a","c","g","t"),labels=c(expression(A%->%G),expression(C%->%"T"),expression(G%->%A),expression("T"%->%C)))+
        theme(legend.title = element_blank())
ggsave("Output1A/SelCoeff/Mean.SC_indivRatevs.SingleRate.comparison.pdf", width = 5,height = 3.5)        

 
## Add average:
nt=c("a","c","g","t")
summary2<-data.frame(nt=c("A","C","G","T"))
for (i in 1:4){
        summary2$syn[i]<-mean(sc$EstSC[sc$ref==nt[i]&sc$Type=="syn"])
        summary2$nonsyn[i]<-mean(sc$EstSC[sc$ref==nt[i]&sc$Type=="nonsyn"])
}
summary2M<-melt(summary2)

ggplot()+
         geom_boxplot(data=scdata,aes(x=base, y=S.C., fill=factor(type)), outlier.alpha =.4, outlier.color = "gray60")+
        scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.0001,0.001,0.01,0.1), minor_breaks=ybreaks, labels=label_scientific2)+
        geom_point(data=summary2M, aes(x=nt, y=value, color=variable), size=3, position=position_dodge(width=.75))+
        labs(x="Nucleotide",y="Estimated selection coefficient")+
        scale_color_manual(values=colors2[c(5,1)] )+
        scale_fill_manual(values=col2_light[c(5,1)] ) + theme_bw()+
        theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10))+
        theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")
ggsave("Output1A/SelCoeff/SC.byNT_withMean_indivMutRates.pdf", width = 5,height = 4)

####
## wilcox test of SC by NT
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))
WilcoxTest.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt)<-c("NT1","NT2","test","P.value")

for (i in 1:nrow(Ncomb)) {
        vec1<-sc$EstSC[sc$ref==Ncomb[i,1]]
        vec2<-sc$EstSC[sc$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
        
        WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt$test[i]<-"greater"
        WilcoxTest.nt$P.value[i]<-result[[3]]
}   

WilcoxTest.nt2<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","P.value")

for (i in 1:nrow(Ncomb)) {
        vec1<-sc$EstSC[sc$ref==Ncomb[i,1]]
        vec2<-sc$EstSC[sc$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt2$test[i]<-"less"
        WilcoxTest.nt2$P.value[i]<-result[[3]]
}   

WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
write.csv(WilcoxTest.nt,"Output1A/SelCoeff/WilcoxTes_byNT_indivMutRates.csv")

