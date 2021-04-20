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

#add gene ingo
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

###
#SC means and SEs by Gene
genenames<-genes$Gene[2:12]
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")
scDF<-sc[,c("gene","EstSC","Depth")]
SESC<-data.frame(gene=genenames)
for (i in 1:11){
        df<-scDF[scDF$gene==genenames[i],]
        SESC$se[i]<-sqrt(mean(df$EstSC)*(1-mean(df$EstSC))/sum(df$Depth))
}


#AT content of each genes:
ATcontents<-list()
for (i in 1:11){
        genename<-genes$Gene[i+1]
        seq1<-as.character(sc$ref[sc$gene==genename])
        ATcontents[i]<-1-GC(seq1)
        names(ATcontents)[i]<-genename
        
}

ATs<-as.data.frame(do.call(rbind,ATcontents))
ATs$gene<-rownames(ATs)
colnames(ATs)[1]<-"AT"


SC_summary2<-merge(AveSC,SESC,by="gene")
SC_summary2<-merge(SC_summary2,ATs,by="gene")
write.csv(SC_summary2, "Output/SelCoeff/SC-AT-CpG_Summary_updated.csv")

#SC_summary2<-read.csv("Output/SelCoeff/SC-AT-CpG_Summary_updated.csv", stringsAsFactors = F, row.names = 1)
SC_summary2$gene<-factor(SC_summary2$gene, levels=c("Core", "E1", "HVR1", "E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ylabel<-expression(paste("Average cost (", italic("s"), ")"))
scaleFUN2<-function(x) sprintf("%.4f", x)
scPlot<-ggplot(SC_summary2,aes(gene,y=SC))+
        geom_point(color="royalblue", size=3)+
        geom_errorbar(aes(ymin=pmax(SC-se, 0), ymax=SC+se), width=.3,color="royalblue" )+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab(ylabel)+
        theme(panel.grid.major.x = element_blank())+
        scale_y_continuous(labels=scaleFUN2)

ATplot<-ggplot(SC_summary2,aes(x=gene,y=AT))+
        geom_point(color="#FF7A16", size=3.5)+
        theme_bw()+ylim(0.32,0.45)+
        theme(axis.title.x=element_blank())+ylab("AT content")+
        theme(panel.grid.major.x = element_blank())


plot_grid(scPlot, ATplot,nrow = 2, labels = c("", ""),
          rel_heights = c(1, 1))
ggsave(filename="Output/SelCoeff/SC-AT.pdf",width =7, height =4.5)


## Correlation:
results1<-cor.test(SC_summary2$AT,SC_summary2$SC, method = "spearman")
print(results1)
#      rho 
#0.2 
#p-value =  0.5576

