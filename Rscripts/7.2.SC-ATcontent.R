#Check the AT contents vs. SC values

library(ggplot2)
library(reshape2)
library(colorspace)
library(cowplot)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
###########

df<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
#coding regions only
df<-df[df$pos>=342,]

#add gene info
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
#SC means & SEs by Gene
genenames<-genes$Gene[2:12]
AveSC<-aggregate(sc$EstSC,by=list(sc$gene), mean, na.rm = TRUE)
colnames(AveSC)<-c("gene","SC")
se<-aggregate(sc$se,by=list(sc$gene), mean, na.rm = TRUE)
colnames(se)<-c("gene","SE")
AveSC$SE<-se$SE

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

SC_summary<-merge(AveSC,ATs,by="gene")
write.csv(SC_summary2, "Output/SelCoeff/SC-AT-Summary.csv")

SC_summary$gene<-factor(SC_summary$gene, levels=c("Core", "E1", "HVR1", "E2","NS1","NS2","NS3","NS4A","NS4B","NS5A","NS5B"))


ylabel<-expression(paste("Average cost (", italic("s"), ")"))
scaleFUN1<-function(x) sprintf("%.1f", x)
scPlot<-ggplot(SC_summary,aes(x=gene,y=SC))+
        geom_point(color="royalblue", size=3)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab(ylabel)+
        theme(panel.grid.major.x = element_blank())+
        #scale_y_continuous(labels=scaleFUN2)
        scale_y_continuous(label=label_scientific)

atPlot<-ggplot(SC_summary,aes(x=gene,y=AT*100))+
        geom_point(color="#FF7A16", size=3.5)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("% AT content")+
        theme(panel.grid.major.x = element_blank())+
        #scale_y_continuous(breaks=c(35,40,45),label=c("35.0", "40.0", "45.0"), limits = c(32,45))
        scale_y_continuous(labels=scaleFUN1,limits = c(32,45))

png("Output/SelCoeff/SC-AT.png",width=7, height=4.5, units="in",res=300)
ggdraw()+
        draw_plot(scPlot,x=0,y=0.5,width=1, height=0.5)+
        draw_plot(atPlot,x=0.05,y=0,width=0.95, height=0.5)+
        draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = 14,fontface = "plain")
dev.off()    



## Correlation:
results1<-cor.test(SC_summary$AT,SC_summary$SC, method = "spearman")
print(results1)
#      rho 
#0.2 
#p-value =  0.5576

