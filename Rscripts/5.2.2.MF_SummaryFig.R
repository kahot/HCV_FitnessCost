library(plotrix)
library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(colorspace)
library(ggplot2)
source("Rscripts/label_scientific.R")

colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")
col2_light1<-qualitative_hcl(6, palette="Set2")
cold<-qualitative_hcl(5, "Cold")


#read in files
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq")
mf.files<-list()
for (i in 1:5){
        data<-read.csv(paste0("Output/MutFreq.filtered/Filtered.", files[i],".Q35.csv"),row.names = 1,stringsAsFactors = F)
        assign(paste0(files[i]),data)
        mf.files[[i]]<-data
        names(mf.files)[i]<-files[i]
}

#summary of mut freq and se by mutation type
TB3<-read.csv("Output/MutFreq.filtered/MF.Mean.SE.summary_updated.csv", stringsAsFactors = F, row.names = 1)


######
## Create a summary plot
Ts<-mf.files[[1]]
Ts<-Ts[Ts$pos>=342, ]

summary<-data.frame(Mutation=rep("Transition", times=4), Type=c("All", "Syn", "Nonsyn", "Nonsense"))
summary$Mean<-TB3$mean[1:4]
summary$SE<-TB3$se[1:4]

summary2<-data.frame(Mutation=rep("Transversion", times=4), Type=c("All", "Syn", "Nonsyn", "Nonsense"))
summary2$Mean[1]<-TB3$mean[TB3$type=="Tvs"]
summary2$Mean[2:4]<-TB3$mean[15:17]
summary2$SE[1]<-TB3$se[TB3$type=="Tvs"]
summary2$SE[2:4]<-TB3$se[15:17]
Summary<-rbind(summary, summary2)
Summary$Type<-factor(Summary$Type, levels=c("All", "Syn","Nonsyn","Nonsense"))


ggplot()+
    geom_rect(data=Summary[Summary$Mutation=="Transition",], aes(xmax=as.numeric(Type)-.35, xmin=as.numeric(Type)-.05, 
                                    ymax=Mean, ymin=0), fill=paste0(colors2[5],"CC")) + 
    scale_x_discrete(breaks=levels(Summary$Type))+    
    geom_rect(data=Summary[Summary$Mutation=="Transversion",], aes(xmax=as.numeric(Type)+.05, xmin=as.numeric(Type)+.35, 
                                                             ymax=Mean, ymin=0), fill=paste0(colors2[6],"CC")) + 
    scale_color_manual(values=paste0(colors2[c(5,6)],"E6"))+
    scale_fill_manual(values=paste0(colors2[c(5,6)],"E6"), labels=c("Transition","Transversion"))+
    geom_point(data=Summary, aes(x=Type,y=Mean,fill=Mutation, color=Mutation),position=position_dodge(.78),size=0.5)+
    geom_errorbar(data=Summary, aes(x=Type,y=Mean,ymin=pmax(Mean-SE, 0) , ymax=Mean+SE,fill=Mutation), position=position_dodge(.78), width=.2, color="gray30")+
    theme_linedraw()+
    theme(axis.title.x=element_blank())+ylab("Mean mutation frequency")+
    scale_y_continuous(trans = 'log10', labels=label_scientific, limits=c(0.0001,0.015))+
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1),axis.title.y = element_text(size=12))+
    geom_vline(xintercept = c(1:3)+0.5,  
               color = "gray60", size=.4)+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype=2, colour="gray60"),
          panel.grid.minor.y = element_line(linetype=2, colour="gray60"),
          legend.title = element_blank(),
          legend.text = element_text(size=11))+
    annotate("segment", x = 0.8, xend = 1.2, y = 0.0115, yend = 0.0115, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.8, xend = 0.8, y = 0.0085, yend = 0.0115, colour = "gray60", size=0.3)+
    annotate("segment", x = 1.2, xend = 1.2, y = 0.002, yend = 0.0115, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.0119, colour = "gray60", label="***", size=3)+
    annotate("segment", x = 1.8, xend = 2.8, y = 0.0117, yend = 0.0117, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 1.8, xend = 1.8, y = 0.009, yend = 0.0117, colour  = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.8, xend = 2.8, y = 0.005, yend = 0.0117, colour  = "#A5C0D9", size=0.3)+
    annotate("text", x = 2.3, y = 0.01199, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.85, xend = 3.8, y = 0.0117, yend = 0.0117, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.85, xend = 2.85, y = 0.005, yend = 0.0117, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 3.8, xend = 3.8, y = 0.003, yend = 0.0117,   colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 3.3, y = 0.01199, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.2, xend = 3.2, y = 0.004, yend = 0.004, colour =  "#CEBCDE", size=0.3)+
    annotate("segment", x = 2.2, xend = 2.2, y = 0.0015, yend = 0.004, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 3.2, xend = 3.2, y = 0.0005, yend = 0.004, colour = "#CEBCDE", size=0.3)+
    annotate("text", x = 2.8, y = 0.00425, colour = "#CEBCDE", label="***", size=3)
    
ggsave("Output/MutFreq.filtered/MeanMF.byTyep_log_stats.pdf", width = 5, height = 3.4)


#Non-log
ggplot()+
    geom_bar(data=Summary, aes(x=Type,y=Mean,fill=Mutation),position=position_dodge(.9), stat="identity",width=0.8)+
    scale_fill_manual(values=paste0(colors2[c(5,6)],"E6"), labels=c("Transition","Transversion"))+
    geom_errorbar(data=Summary, aes(x=Type,y=Mean,ymin=pmax(Mean-SE, 0) , ymax=Mean+SE,fill=Mutation), position=position_dodge(.9), width=.2, color="gray30")+
    theme_linedraw()+
    theme(axis.title.x=element_blank())+ylab("Mean mutation frequency")+
    ylim(0,0.01)+
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1),axis.title.y = element_text(size=12))+
    geom_vline(xintercept = c(1:3)+0.5,  
               color = "gray60", size=.4)+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype=2, colour="gray60"),
          panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=11))+
    annotate("segment", x = 0.8, xend = 1.2, y = 0.0085, yend = 0.0085, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.8, xend = 0.8, y = 0.006, yend = 0.0085, colour =  "gray60", size=0.3)+
    annotate("segment", x = 1.2, xend = 1.2, y = 0.002, yend = 0.0085, colour =  "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.0087, colour = "gray60", label="***", size=3)+
    annotate("segment", x = 1.8, xend = 2.75, y = 0.0095, yend = 0.0095, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 1.8, xend = 1.8, y = 0.0085,yend = 0.0095, colour =   "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.75, xend = 2.75, y = 0.005, yend = 0.0095, colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 2.3, y = 0.0097, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.85, xend = 3.75, y = 0.0095, yend = 0.0095, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.85, xend = 2.85, y = 0.005,  yend = 0.0095, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 3.75, xend = 3.75, y = 0.0026, yend = 0.0095, colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 3.3, y = 0.0097, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.25, xend = 3.2, y = 0.0036, yend = 0.0036, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 2.25, xend = 2.25, y = 0.0015, yend = 0.0036, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 3.2, xend = 3.2, y = 0.001, yend = 0.0036, colour = "#CEBCDE", size=0.3)+
    annotate("text", x = 2.75, y = 0.0039, colour = "#CEBCDE", label="***", size=3)

ggsave("Output/MutFreq.filtered/MeanMF.byTyep_stats.pdf", width = 5, height = 3.4)






####################################
#Plot summary of mutation frequency by type by nucleotide       

HCVFiles3<-list.files("Output/Overview3/",pattern="overview3.csv")
s<-length(HCVFiles3)
#exclude CpG creating mutations
Ts2<-Ts[Ts$makesCpG==0,]
k=1
transMF<-list()
for (i in c("a","t","c","g")) {
    for (type in c("syn","nonsyn")){
        datavector<-Ts$mean[Ts$Type==type & Ts$ref==i]
        nt<-toupper(i)
        vname<-paste0(nt,".",type)
        dat<-data.frame(base=rep(nt,times=length(datavector)),
                        type=rep(type, times=length(datavector)), MF=datavector)
        transMF[[k]]<-dat
        names(transMF)[k]<-vname
        k=k+1
    }
}

mfdata2<-do.call(rbind, transMF)

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 

ggplot(mfdata2,aes(x=base, y=MF, fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
    scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
    labs(x="",y="Mutation frequency")+
    scale_fill_manual(values=col2_light1[c(5,1)]) + theme_bw()+
    theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10, color=1))+
    theme(axis.text.y = element_text(size =10), axis.title.y= element_text(size =12), panel.grid.major.x=element_blank())+
    geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
    scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))

ggsave("Output/MutFreq.filtered/MF.byNT._noCpG.pdf", width = 6,height = 4)


####### Create plots of mut frq by gene
## By genes
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

depth<-read.csv("Output/ReadDepth_sum.csv",stringsAsFactors = F, row.names = 1)

mfs<-merge(Ts,depth, by="pos", all.x=T)

mf1<-mfs[!is.na(mfs$mean),]
mf1<-merge(mf1, genetable, by="pos", all.x=T )

SumMFGenes<-aggregate(mf1$mean,by=list(mf1$gene),FUN=mean)
colnames(SumMFGenes)<-c("Gene", "Mean")
meSE<-data.frame(gene=genenames)

for (i in 1:11){
    df<-mf1[mf1$gene==genenames[i],]
    meSE$SE[i]<-sqrt(mean(df$mean)*(1-mean(df$mean))/mean(df$Depth))
}

for (i in 1:11){
    df<-mf1[mf1$gene==genenames[i],]
    meSE$SE2[i]<-sqrt(mean(df$mean)*(1-mean(df$mean))/nrow(df))
}




sumG<-cbind(SumMFGenes, meSE$SE)
colnames(sumG)[3]<-"SE"
sumG$Gene<-factor(sumG$Gene, levels=genenames)

write.csv(sumG, "Output/MutFreq.filtered/MF_Summary_Table.by.gene.csv")

#
mf2<-mf1[,c("pos","mean","gene","Type")]
mf2$gene<-factor(mf2$gene, levels=genenames)

ybreaks<-c(seq(0.0005, 0.0009, by=0.0001),seq(0.001,0.009, by=0.001), seq(0.01, 0.04, by=0.01))
ylabel<-c(rep("", times=5),0.001, 0.002, 0.003,"", 0.005,"", 0.007, "","",0.01, 0.02, 0.03,"")
ggplot(sumG, aes(x=Gene, y=Mean))+
    scale_y_continuous(trans="log10",breaks = ybreaks,labels=ylabel)+
    geom_point(data=mf2, aes(x=gene, y=mean), position=position_jitter(width=0.2),stat = "identity", col=paste0(colors2[5],"4D"), size=0.2)+
    geom_point(color="blue3")+
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.25, position=position_dodge(width=0.3),color="blue3")+
    theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+
    ylab("Mutation frequency")+
    theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13), panel.grid.minor.y = element_blank())+
    geom_vline(xintercept = c(1:10)+0.5,  
               color = "gray70", size=.5)+
    theme(axis.title.x=element_blank())
ggsave(filename="Output/SummaryFig.Filtered/Ave.mf.byGene_updated.pdf", width = 7, height = 4)


######## Plot mut freq by Gene and by type (Syn & nonsyn)
#exclude stop mutations
mf1<-mf1[mf1$Type!="stop",]

SumMFGenes<-aggregate(mf1$mean,by=list(mf1$gene, mf1$Type),FUN=mean)
colnames(SumMFGenes)<-c("Gene","Type","Mean")

SE1<-data.frame(Gene=genenames, Type="syn")
SE2<-data.frame(Gene=genenames, Type="nonsyn")
for (i in 1:11){
    df1<-mf1[mf1$gene==genenames[i]&mf1$Type=="syn",]
    df2<-mf1[mf1$gene==genenames[i]&mf1$Type=="nonsyn",]
    
    SE1$SE[i]<-sqrt(mean(df1$mean)*(1-mean(df1$mean))/sum(df1$Depth))
    SE2$SE[i]<-sqrt(mean(df2$mean)*(1-mean(df2$mean))/sum(df2$Depth))
}
seSum<-rbind(SE2,SE1)
sumG2<-merge(SumMFGenes, seSum, by=c("Gene","Type"))
sumG2$Gene<-factor(sumG2$Gene, levels=genenames)
write.csv(sumG2, "Output/MutFreq.filtered/MF_Summary_Table.by.gene.byType.csv")

#ggplot(sumG2, aes(x=Gene, y=Mean, group=Type, color=Type))+
#    geom_point(position=position_dodge(width=0.3))+scale_color_manual(values=colors2[c(1,5)])+
#    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.4, position=position_dodge(width=0.3))+
#    theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Mutation frequency")+
#    theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
#    geom_vline(xintercept = c(1:10)+0.5,  
#               color = "gray70", size=.5)+
#    theme(axis.title.x=element_blank())
#ggsave(filename="Output/SummaryFig.Filtered/Ave.MF_by.gene_by.type_updated.pdf", width = 8.5, height = 5)

###make it into one figure
mf1$gene<-factor(mf1$gene, levels=genenames)


ggplot(sumG2, aes(x=Gene, y=Mean))+
    scale_y_continuous(trans="log10",breaks = ybreaks,labels=ylabel)+
    geom_point(data=mf1, aes(x=gene, y=mean, color=Type), position=position_jitterdodge(dodge.width=0.6, jitter.width=0.15),
               alpha=0.1,stat = "identity", size=0.2)+
    geom_point(data=sumG2, aes(x=Gene, y=Mean, group=Type, color=Type), color="gray60", position=position_dodge(width=0.6))+
    geom_errorbar(data=sumG2, aes(ymin=Mean-SE, ymax=Mean+SE,color=Type), width=.6, position=position_dodge(width=0.6))+
    theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+
    ylab("Mutation frequency")+
    theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13), panel.grid.minor.y = element_blank())+
    geom_vline(xintercept = c(1:10)+0.5,  
               color = "gray70", size=.5)+
    theme(axis.title.x=element_blank())+
    scale_color_manual(values=colors2[c(1,5)], labels=c("Nonsyn","Syn"))+
    theme(legend.title=element_blank())+
    guides(fill = guide_legend(override.aes = list(linetype = 0)),
           color = guide_legend(override.aes = list(linetype = 0)))+
    annotate("text", x=c(1,2,4:11), y=rep(0.035, times=10), label="*",size=4)

ggsave("Output/MutFreq.filtered//MutFreq.byGene.byType2.pdf", width=8, height = 4)       








