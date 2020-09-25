library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)
library(cowplot)
library(gridExtra)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")
###########

df<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
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
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")
SESC<-aggregate(sc$EstSC,by=list(sc$gene), std.error, na.rm = TRUE)
colnames(SESC)<-c("gene","se")

colors2<-qualitative_hcl(6, palette="Dark3")
scCols<-c("#E16A86","#009ADE")

## A

dt1<-sc[sc$ref=="a",]
A<-ggplot(dt1, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt1, Type=="nonsyn"),  color="gray60", alpha = 0.5,fill=scCols[1])+
        geom_histogram(data=subset(dt1, Type=="syn"),color="gray60", alpha = 0.5, fill=scCols[2])+
        scale_fill_manual(values=scCols)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("A")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt1$EstSC[dt1$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt1$EstSC[dt1$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)+
        theme(legend.position = "none")

dt3<-sc[sc$ref=="t",]
T_<-ggplot(dt3, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt3, Type=="nonsyn"),  color="gray60", alpha = 0.5,fill=scCols[1])+
        geom_histogram(data=subset(dt3, Type=="syn"),color="gray60", alpha = 0.5, fill=scCols[2])+
        scale_fill_manual(values=scCols)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("T")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt3$EstSC[dt3$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt3$EstSC[dt3$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)+
        theme(legend.position = "none")

dt4<-sc[sc$ref=="c" & sc$Type!='stop',]
mean1<-mean(dt4$EstSC[dt4$Type=="syn"], na.rm = T)
mean2<-mean(dt4$EstSC[dt4$Type=="nonsyn"], na.rm = T)

C<-ggplot(dt4, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt4, Type=="nonsyn"),  color="gray60", alpha = 0.5,fill=scCols[1])+
        geom_histogram(data=subset(dt4, Type=="syn"),color="gray60", alpha = 0.5, fill=scCols[2])+
        
        scale_fill_manual(values=scCols)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("C")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt4$EstSC[dt4$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt4$EstSC[dt4$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)+
        theme(legend.position = "none")


dt2<-sc[sc$ref=="g"&sc$Type!='stop',]

G<-ggplot(dt2, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt2, Type=="nonsyn"), aes(x=EstSC, fill=Type), color="gray60", alpha = 0.5)+
        geom_histogram(data=subset(dt2, Type=="syn"),aes(x=EstSC,fill=Type), color="gray60", alpha = 0.5)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        scale_fill_manual(name='', values=scCols,labels=c("Nonsyn","Syn"))+
        scale_color_manual(name='', values=scCols,labels=c("Nonsyn","Syn"))+
        theme_bw()+ylab("Count")+xlab("")+ggtitle("G")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt2$EstSC[dt2$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt2$EstSC[dt2$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)


plot_grid(A,T_,C,G, nrow = 1, ncol=4, rel_widths=c(1,1,1,1.5))
#grid.arrange(A,T_,C,G,  nrow = 1)
ggsave("Output1A/SelCoeff/SC_DFEplot.pdf", width=9, height=2.5)
