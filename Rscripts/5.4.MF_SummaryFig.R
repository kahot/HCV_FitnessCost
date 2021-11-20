# Create mutation frequency figures
library(reshape2)
library(colorspace)
library(ggplot2)
source("Rscripts/label_scientific.R")

colors2<-qualitative_hcl(6, palette="Dark3",)
col2_light1<-qualitative_hcl(6, palette="Set2")

## Create a summary plot 
#summary of mut freq and se by mutation type
TB3<-read.csv("Output/MutFreq/MF.Mean.SE.summary.csv", stringsAsFactors = F, row.names = 1)
#Reorganize the data frame
TB3<-TB3[-grep("Tv1|Tv2",TB3$type),]
TB3<-TB3[-grep("all.|All",TB3$type),]

TB3$Type<-''
TB3$Type[grep("syn",TB3$type)]<-"Syn"
TB3$Type[grep("ns",TB3$type)]<-"Nonsyn"
TB3$Type[grep("stop",TB3$type)]<-"Nonsense"
TB3$Type[TB3$Type=='']<-"All"
TB3$Mutation<-''
TB3$Mutation[grep("Ts",TB3$type)]<-"Transition"
TB3$Mutation[grep("tvs|Tvs",TB3$type)]<-"Tranversion"
TB3$Mutation[grep("Ts",TB3$type)]<-"Transition"

Summary<-TB3

#With 95% CI
#Nonsyn-tvs is too small for lower CI. Convert and make it to a dummy number to show the error bar.

Summary$CI_down[Summary$type=="tvs.ns"]<-0
Summary$CI_down[Summary$type=="tvs.stop"]<-0

#Order the xaxis
Summary$type<-factor(Summary$type, levels=c("Ts","Tvs","Ts.syn","tvs.syn","Ts.ns","tvs.ns","Ts.stop","tvs.stop"))
Summary$Type<-factor(Summary$Type, levels=c("All","Syn","Nonsyn","Nonsense"))

ggplot()+
    geom_errorbar(data=Summary, aes(x=Type,y=mean,ymin=CI_down , ymax=CI_up,fill=Mutation), 
                  position=position_dodge(width=.78),size=0.3, width=.2, color="gray30")+
    geom_point(data=Summary, aes(x=Type, y=mean, fill=Mutation),
               position=position_dodge(width=.78), color="gray30",shape=21, size=3)+
    theme_classic()+
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color=1))+
    theme(axis.title.x=element_blank())+ylab("Mean mutation freq. ± 95% C.I.")+
    scale_fill_manual(values=colors2[c(3,2)],  guide="none",labels=c("Transition","Transversion"))+
    geom_vline(xintercept = c(1:3)+0.5, color = "gray60", size=.4)+
    scale_y_continuous(trans = 'log10', labels=label_scientific, limits=c(0.0001,0.017))+
    annotate("segment", x = 0.8, xend = 1.2, y = 0.0115, yend = 0.0115, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.8, xend = 0.8, y = 0.0085, yend = 0.0115, colour = "gray60", size=0.3)+
    annotate("segment", x = 1.2, xend = 1.2, y = 0.0022, yend = 0.0115, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.0119, colour = "gray60", label="***", size=3)+
    annotate("segment", x = 1.8, xend = 2.8, y = 0.0145, yend = 0.0145, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 1.8, xend = 1.8, y = 0.012,  yend = 0.0145, colour  = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.8, xend = 2.8, y = 0.0075, yend = 0.0145, colour  = "#A5C0D9", size=0.3)+
    annotate("text", x = 2.3, y = 0.0148, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.85, xend = 3.8, y = 0.0145, yend = 0.0145, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.85, xend = 2.85,y = 0.0075, yend = 0.0145, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 3.8, xend = 3.8,  y = 0.0045, yend = 0.0145,   colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 3.3, y = 0.0148, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.2, xend = 3.2, y = 0.0062, yend = 0.0062, colour =  "#CEBCDE", size=0.3)+
    annotate("segment", x = 2.2, xend = 2.2, y = 0.0022, yend = 0.0062, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 3.2, xend = 3.2, y = 0.0011, yend = 0.0062, colour = "#CEBCDE", size=0.3)+
    annotate("text", x = 2.7, y = 0.0063, colour = "#CEBCDE", label="***", size=3)+
    annotate("text", x=4.3, y=0.015, color="gray20", label="Ts", size=4, hjust=0)+
    annotate("text", x=4.3, y=0.011, color="gray20", label="Tvs",size=4, hjust=0)+
    annotate("point",x=4.2, y=0.015, fill=colors2[3],shape=21, color="gray30",size=2.5)+
    annotate("point",x=4.2, y=0.011, fill=colors2[2],shape=21, color="gray30",size=2.5)
ggsave("Output/MutFreq/MF.byTyep_log_stats_CI.pdf", width = 4.8, height = 3.7)



#### Create a plot of mut freq. by gene and by type  

HCVFiles3<-list.files("Output/OverviewF/",pattern="overviewF.csv")
s<-length(HCVFiles3)
#read the transition mutation frequency file
Ts<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",row.names = 1,stringsAsFactors = F)
Ts<-Ts[Ts$pos>341&Ts$pos<8575,]
Ts<-Ts[Ts$Type!="stop",]

#Add gene information to Ts file
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
    gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[2:12]
mf1<-merge(Ts, genetable, by="pos", all.x=T )

#Read the SE file
se<-read.csv("Output/MutFreq/SEmatrix_Ts.csv",stringsAsFactors = F, row.names = 1)
# remove the sites not in  mf1
se<-se[se$pos %in% mf1$pos,]
#calculate mean per site
se$ave<-rowMeans(se[2:196],na.rm=T)
#attach the gene & type info
se<-merge(genetable, se[,c("pos","ave")], by="pos" )
se$Type<-mf1$Type

#Calculate the mean of mut. freq. and SE
GeneSummary<-aggregate(mf1$mean,by=list(mf1$gene, mf1$Type),FUN=mean)
colnames(GeneSummary)<-c("Gene","Type","Mean")

SE<-aggregate(se$ave,by=list(se$gene,se$Type),FUN=mean)
colnames(SE)<-c("Gene","Type", "SE")

GeneSummary$SE<-SE$SE

GeneSummary$CI<-1.96*GeneSummary$SE
write.csv(GeneSummary, "Output/MutFreq/MF.Summary.by.gene.by.type.csv")


GeneSummary$Gene<-factor(GeneSummary$Gene, levels=genenames)
GeneSummary$Type<-factor(GeneSummary$Type,levels=c("syn","nonsyn"))

mf2<-mf1[,c("pos","mean","gene","Type")]
mf2$gene<-factor(mf2$gene, levels=genenames)
mf2$Type<-factor(mf2$Type,levels=c("syn","nonsyn"))

ybreaks<- c(10^(-4),10^(-3),10^(-2),10^(-1)) 
ggplot()+
    geom_violin(data=mf2, aes(x=gene, y=mean, color=Type, fill=Type), size=.2,trim=F)+
    scale_y_continuous(trans="log10",breaks=ybreaks, labels=label_scientific)+
    geom_point(data=GeneSummary, aes(x=Gene, y=Mean, color=Type),size=1,
               position = position_dodge(width = 1))+
    geom_errorbar(data=GeneSummary, aes(x=Gene,y=Mean,ymin=pmax(Mean-CI, 0) , ymax=Mean+CI,fill=Type), 
                  position=position_dodge(1), size=0.2, width=.2, color="gray30")+
    xlab('')+ylab('Mutation frequency')+
    geom_point(data=GeneSummary, aes(x=Gene, y=Mean, fill=Type),size=1,color="gray40",
               position = position_dodge(width = 1))+
    theme_classic()+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = 1:10+0.5,  color = "gray80", size=.4)+
    scale_color_manual(values=colors2[c(5,1)])+
    scale_fill_manual(values=col2_light1[c(5,1)])+
    theme(legend.title = element_blank())+
    annotate("text", x=c(1,2,4:11), y=rep(0.0003, times=10), label="***",size=3.6, color="slategray3")

ggsave("Output/MutFreq/MF_violinePlot.byGene.pdf", width = 7, height = 5)







