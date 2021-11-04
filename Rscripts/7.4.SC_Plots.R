library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)
library(cowplot)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")
###########

df<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
#coding regions only
df<-df[df$pos>=342,]


#### Calcualte the mean SC and GC_contents for each gene
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
write.csv(sc,"Output/SelCoeff/SC_gene.csv") 

#Average mean SC
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")

#SE
genenames<-genes$Gene[2:12]

SESC2<-data.frame(gene=genenames)
for (i in 1:11){
        df<-sc[sc$gene==genenames[i],]
        SESC2$se[i]<-sqrt(mean(df$EstSC)*(1-mean(df$EstSC))/mean(df$Depth))
}


AveSC$se<-SESC2$se
write.csv(AveSC, "Output/SelCoeff/SC_ave_se_by_gene.csv")

#######
##########################################
#Plot summary of 1) sel coef by type by nucleotide       
#2 no CPG creating mutations
sc3<-sc[sc$makesCpG==0,]
k=1
transSC<-list()
for (i in c("a","t","c","g")) {
    for (type in c("syn","nonsyn")){
        datavector<-sc3$EstSC[sc3$Type==type & sc3$ref==i]
        nt<-toupper(i)
        vname<-paste0(nt,".",type)
        dat<-data.frame(base=rep(nt,times=length(datavector)),
                        type=rep(type, times=length(datavector)), S.C.=datavector)
        transSC[[k]]<-dat
        names(transSC)[k]<-vname
        k=k+1
    }
}

scdata2<-do.call(rbind, transSC)

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 

ggplot(scdata2,aes(x=base, y=S.C., fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
    scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
    labs(x="Nucleotide",y="Estimated selection coefficient")+
    scale_fill_manual(values=colors2[c(5,1)]) + theme_bw()+
    theme(legend.title = element_blank())+theme(axis.text.x = element_text(size =10, color=1), axis.title.x = element_blank())+
    theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
    geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
    scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))

ggsave("Output/SelCoeff/SC.byNT_noCpG.pdf", width = 6,height = 4)



####
### Plot the effect size ###
effects<-read.csv("Output/BetaReg/BetaReg.effects.csv",stringsAsFactors = F)
effects$factor<-factor(effects$factor, levels=effects$factor[1:17])

effects$percent<-effects$percent*100
effects$percent_SC<-effects$percent_SC*100


##effect size for SC
ggplot(effects, aes(factor,percent_SC)) +
    geom_bar(stat="identity", color=colors2[5], fill=paste0(colors2[5],"CC"))+
    theme_classic() +
    theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")

ggsave("Output/SummaryFig.Filtered/BetaReg.effects.SC.pdf", width = 9, height = 6)


##effect size for MF & SC
colnames(effects)[4:5]<-c("MF","SC")
effects2<-melt(effects[,c(2,4,5)], id.var=c("factor"))
colnames(effects2)[2:3]<-c("Type","percent")

ggplot(effects2, aes(factor,percent, group=Type, fill=Type)) +
    scale_fill_manual(values=c(paste0(colors2[4],"99"),paste0(colors2[5],"CC")), labels=c("Mut freq","Sel coeff"))+
    scale_color_manual(values=c(colors2[4], colors2[5]))+
    geom_bar(stat="identity",color="gray",position=position_dodge(width=.5))+
    theme_test() +
    theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")+
    theme(legend.title = element_blank())

ggsave("Output/SummaryFig.Filtered/BetaReg.effects.MF-SC2.pdf", width = 9, height = 5)

