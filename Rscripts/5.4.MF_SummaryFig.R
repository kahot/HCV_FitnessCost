# Create mutation frequency figures
library(reshape2)
library(colorspace)
library(ggplot2)
source("Rscripts/label_scientific.R")

colors2<-qualitative_hcl(6, palette="Dark3")
col8<-qualitative_hcl(8, palette="Dark3")
col2_light1<-qualitative_hcl(6, palette="Set2")

## Create a summary plot 
#summary of mut freq and CI by mutation type
mf<-read.csv("Output/MutFreq/MF_summary.csv", stringsAsFactors = F, row.names = 1)
colnames(mf)<-c("Mutation","All","Syn","Nonsyn",'Nonsense')
mf<-mf[mf$Mutation=="Ts"|mf$Mutation=="Tvs",]

mfm<-melt(mf, id.vars="Mutation")
mfm$id<-paste0(mfm$Mutation,".",mfm$variable)
mfm$id<-gsub(".All","", mfm$id )
colnames(mfm)[2:3]<-c("Type","MF")

Sum<-read.csv("Output/MutFreq/CI_TypeofMutaions_summary.csv", stringsAsFactors = F, row.names = 1)
#consistent variable
Sum$Type[Sum$Type=="all"]<-"All"
Sum$Type[Sum$Type=="syn"]<-"Syn"
Sum$Type[Sum$Type=="nonsyn"]<-"Nonsyn"
Sum$Type[Sum$Type=="stop"]<-"Nonsense"
Sum$id<-paste0(Sum$Mutation,".",Sum$Type)
Sum$id<-gsub(".All","", Sum$id )

mfmF<-merge(mfm, sumF[,c(2,3,6)], id.vars="id")
mfmAC<-merge(mfm, sumAC[,c(2,3,6)], id.vars="id")
sumF2<-Sum[Sum$Method=="no.rescaling",]

#Order the xaxis
mfm$id<-factor(mfm$id, levels=c("Ts","Tvs","Ts.Syn","Tvs.Syn","Ts.Nonsyn","Tvs.Nonsyn","Ts.Nonsense","Tvs.Nonsense"))
mfm$Type<-factor(mfm$Type, levels=c("All","Syn","Nonsyn","Nonsense"))

mfmF<-merge(mfm, sumF[,c(2,3,6)], id.vars="id")
mfmAC<-merge(mfm, sumAC[,c(2,3,6)], id.vars="id")
mfmF2<-merge(mfm, sumF2[,c(2,3,6)], id.vars="id")

#CIs based on F-distributions with rescaling
ggplot()+
    geom_errorbar(data=mfmF, aes(x=Type,y=MF,ymin=MF-low, ymax=MF+up,fill=Mutation), 
                  position=position_dodge(width=.78),size=0.4, width=.5, color="gray30")+
    geom_point(data=mfmF, aes(x=Type,y=MF, fill=Mutation),
               position=position_dodge(width=.78), color="gray30",shape=21, size=3.2)+
    theme_classic()+
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color=1))+
    theme(axis.title.x=element_blank())+ylab("Mean mutation freq. ± 95% C.I.")+
    scale_fill_manual(values=colors2[c(3,2)],  guide="none",labels=c("Transition","Transversion"))+
    geom_vline(xintercept = c(1:3)+0.5, color = "gray60", size=.4)+
    scale_y_continuous(trans = 'log10', labels=label_scientific)+
    annotate("segment", x = 0.8, xend = 1.2, y = 0.01, yend = 0.01, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.8, xend = 0.8, y = 0.007, yend = 0.01, colour = "gray60", size=0.3)+
    annotate("segment", x = 1.2, xend = 1.2, y = 0.0015, yend = 0.01, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.0105, colour = "gray60", label="***", size=3)+
    annotate("segment", x = 1.8, xend = 2.8, y = 0.012, yend = 0.012, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 1.8, xend = 1.8, y = 0.01,  yend = 0.012, colour  = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.8, xend = 2.8, y = 0.0045, yend = 0.012, colour  = "#A5C0D9", size=0.3)+
    annotate("text", x = 2.3, y = 0.0125, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.85, xend = 3.8, y = 0.012, yend = 0.012, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.85, xend = 2.85,y = 0.0045, yend = 0.012, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 3.8, xend = 3.8,  y = 0.003, yend = 0.012,   colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 3.3, y = 0.0125, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.2, xend = 3.2, y = 0.004, yend = 0.004, colour =  "#CEBCDE", size=0.3)+
    annotate("segment", x = 2.2, xend = 2.2, y = 0.002, yend = 0.004, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 3.2, xend = 3.2, y = 0.001, yend = 0.004, colour = "#CEBCDE", size=0.3)+
    annotate("text", x = 2.7, y = 0.0042, colour = "#CEBCDE", label="***", size=3)+
    annotate("text", x=4.3, y=0.015, color="gray20", label="Ts", size=4, hjust=0)+
    annotate("text", x=4.3, y=0.011, color="gray20", label="Tvs",size=4, hjust=0)+
    annotate("point",x=4.2, y=0.015, fill=colors2[3],shape=21, color="gray30",size=2.5)+
    annotate("point",x=4.2, y=0.011, fill=colors2[2],shape=21, color="gray30",size=2.5)
ggsave("Output/MutFreq/MF.byTyep_log_stats_CI.pdf", width = 4.8, height = 3.7)


#Agresti–Coull
ggplot()+
    geom_errorbar(data=mfmAC, aes(x=Type,y=MF,ymin=MF-low, ymax=MF+up,fill=Mutation), 
                  position=position_dodge(width=.78),size=0.3, width=.2, color="gray30")+
    geom_point(data=mfmAC, aes(x=Type,y=MF, fill=Mutation),
               position=position_dodge(width=.78), color="gray30",shape=21, size=3)+
    theme_classic()+
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color=1))+
    theme(axis.title.x=element_blank())+ylab("Mean mutation freq. ± 95% C.I.")+
    scale_fill_manual(values=colors2[c(3,2)],  guide="none",labels=c("Transition","Transversion"))+
    geom_vline(xintercept = c(1:3)+0.5, color = "gray60", size=.4)+
    scale_y_continuous(trans = 'log10', labels=label_scientific)+
    annotate("segment", x = 0.8, xend = 1.2, y = 0.017, yend = 0.017, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.8, xend = 0.8, y = 0.014, yend = 0.017, colour = "gray60", size=0.3)+
    annotate("segment", x = 1.2, xend = 1.2, y = 0.005, yend = 0.017, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.018, colour = "gray60", label="***", size=3)+
    annotate("segment", x = 1.8, xend = 2.8, y = 0.026, yend = 0.026, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 1.8, xend = 1.8, y = 0.022, yend = 0.026, colour  = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.8, xend = 2.8, y = 0.012, yend = 0.026, colour  = "#A5C0D9", size=0.3)+
    annotate("text", x = 2.3, y = 0.027, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.85, xend = 3.8, y = 0.026, yend = 0.026, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.85, xend = 2.85,y = 0.012, yend = 0.026, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 3.8, xend = 3.8,  y = 0.008,  yend = 0.026,   colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 3.3, y = 0.027, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.2, xend = 3.2, y = 0.010, yend = 0.01, colour =  "#CEBCDE", size=0.3)+
    annotate("segment", x = 2.2, xend = 2.2, y = 0.006, yend = 0.01, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 3.2, xend = 3.2, y = 0.004, yend = 0.01, colour = "#CEBCDE", size=0.3)+
    annotate("text", x = 2.7, y = 0.0105, colour = "#CEBCDE", label="***", size=3)+
    annotate("text", x=4.3, y=0.015, color="gray20", label="Ts", size=4, hjust=0)+
    annotate("text", x=4.3, y=0.011, color="gray20", label="Tvs",size=4, hjust=0)+
    annotate("point",x=4.2, y=0.015, fill=colors2[3],shape=21, color="gray30",size=2.5)+
    annotate("point",x=4.2, y=0.011, fill=colors2[2],shape=21, color="gray30",size=2.5)
ggsave("Output/MutFreq/MF.byTyep_log_stats_CI_agresi-coull.pdf", width = 4.8, height = 3.7)

#F-distribution method without rescaling
ggplot()+
    geom_errorbar(data=mfmF2, aes(x=Type,y=MF,ymin=MF-low, ymax=MF+up,fill=Mutation), 
                  position=position_dodge(width=.78),size=0.3, width=.2, color="gray30")+
    geom_point(data=mfmF2, aes(x=Type,y=MF, fill=Mutation),
               position=position_dodge(width=.78), color="gray30",shape=21, size=3.2)+
    theme_classic()+
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color=1))+
    theme(axis.title.x=element_blank())+ylab("Mean mutation freq. ± 95% C.I.")+
    scale_fill_manual(values=colors2[c(3,2)],  guide="none",labels=c("Transition","Transversion"))+
    geom_vline(xintercept = c(1:3)+0.5, color = "gray60", size=.4)+
    scale_y_continuous(trans = 'log10', labels=label_scientific)+
    annotate("segment", x = 0.8, xend = 1.2, y = 0.017, yend = 0.017, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.8, xend = 0.8, y = 0.014, yend = 0.017, colour = "gray60", size=0.3)+
    annotate("segment", x = 1.2, xend = 1.2, y = 0.005, yend = 0.017, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.018, colour = "gray60", label="***", size=3)+
    annotate("segment", x = 1.8, xend = 2.8, y = 0.026, yend = 0.026, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 1.8, xend = 1.8, y = 0.022, yend = 0.026, colour  = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.8, xend = 2.8, y = 0.012, yend = 0.026, colour  = "#A5C0D9", size=0.3)+
    annotate("text", x = 2.3, y = 0.027, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.85, xend = 3.8, y = 0.026, yend = 0.026, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 2.85, xend = 2.85,y = 0.012, yend = 0.026, colour = "#A5C0D9", size=0.3)+
    annotate("segment", x = 3.8, xend = 3.8,  y = 0.008,  yend = 0.026,   colour = "#A5C0D9", size=0.3)+
    annotate("text", x = 3.3, y = 0.027, colour = "#A5C0D9", label="***", size=3)+
    annotate("segment", x = 2.2, xend = 3.2, y = 0.010, yend = 0.01, colour =  "#CEBCDE", size=0.3)+
    annotate("segment", x = 2.2, xend = 2.2, y = 0.006, yend = 0.01, colour = "#CEBCDE", size=0.3)+
    annotate("segment", x = 3.2, xend = 3.2, y = 0.004, yend = 0.01, colour = "#CEBCDE", size=0.3)+
    annotate("text", x = 2.7, y = 0.0105, colour = "#CEBCDE", label="***", size=3)+
    annotate("text", x=4.3, y=0.015, color="gray20", label="Ts", size=4, hjust=0)+
    annotate("text", x=4.3, y=0.011, color="gray20", label="Tvs",size=4, hjust=0)+
    annotate("point",x=4.2, y=0.015, fill=colors2[3],shape=21, color="gray30",size=2.5)+
    annotate("point",x=4.2, y=0.011, fill=colors2[2],shape=21, color="gray30",size=2.5)
ggsave("Output/MutFreq/MF.byTyep_log_stats_CI_FnoRescaling.pdf", width = 4.8, height = 3.7)




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

#Read the CI file
ci<-read.csv("Output/MutFreq/CI.Ts.byGene.csv",stringsAsFactors = F, row.names = 1)
ci$Type[ci$Type=="Nonsyn"]<-"nonsyn"
ci$Type[ci$Type=="Syn"]<-"syn"
ci$id<-paste0(ci$gene,".",ci$Type)

ci.ac<-read.csv("Output/MutFreq/CI.Ts.byGene_Agresti-coull.csv",stringsAsFactors = F, row.names = 1)
ci.ac$Type[ci.ac$Type=="Nonsyn"]<-"nonsyn"
ci.ac$Type[ci.ac$Type=="Syn"]<-"syn"
ci.ac$id<-paste0(ci.ac$gene,".",ci.ac$Type)

ci.noR<-read.csv("Output/MutFreq/CI.Ts.byGene_noRescaling.csv",stringsAsFactors = F, row.names = 1)
ci.noR$Type[ci.noR$Type=="Nonsyn"]<-"nonsyn"
ci.noR$Type[ci.noR$Type=="Syn"]<-"syn"
ci.noR$id<-paste0(ci.noR$gene,".",ci.noR$Type)

#Calculate the mean of mut. freq.
summary<-aggregate(mf1$mean,by=list(mf1$gene, mf1$Type),FUN=mean)
colnames(summary)<-c("Gene","Type","Mean")
summary$id<-paste0(summary$Gene,".",summary$Type)

GeneSummary<-merge(summary, ci[c("id","low","up")], by = "id")
colnames(GeneSummary)[5:6]<-c("CI.low","CI.up")

GeneSummary<-merge(GeneSummary, ci.ac[c("id","low","up")], by = "id")
colnames(GeneSummary)[7:8]<-c("CI.low.ac","CI.up.ac")

GeneSummary<-merge(GeneSummary, ci.noR[c("id","low","up")], by = "id")
colnames(GeneSummary)[9:10]<-c("CI.low.noR","CI.up.noR")

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
    geom_errorbar(data=GeneSummary, aes(x=Gene,y=Mean,ymin=Mean-CI.low , ymax=Mean+CI.up,fill=Type), 
                  position=position_dodge(1), size=0.2, width=.3, color="gray30")+
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

#with agresi-coull method
ggplot()+
    geom_violin(data=mf2, aes(x=gene, y=mean, color=Type, fill=Type), size=.2,trim=F)+
    scale_y_continuous(trans="log10",breaks=ybreaks, labels=label_scientific)+
    geom_point(data=GeneSummary, aes(x=Gene, y=Mean, color=Type),size=1,
               position = position_dodge(width = 1))+
    geom_errorbar(data=GeneSummary, aes(x=Gene,y=Mean,ymin=Mean-CI.low.ac , ymax=Mean+CI.up.ac,fill=Type), 
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
ggsave("Output/MutFreq/MF_violinePlot.byGene_AC.pdf", width = 7, height = 5)

#without rescaling factors
ggplot()+
    geom_violin(data=mf2, aes(x=gene, y=mean, color=Type, fill=Type), size=.2,trim=F)+
    scale_y_continuous(trans="log10",breaks=ybreaks, labels=label_scientific)+
    geom_point(data=GeneSummary, aes(x=Gene, y=Mean, color=Type),size=1,
               position = position_dodge(width = 1))+
    geom_errorbar(data=GeneSummary, aes(x=Gene,y=Mean,ymin=Mean-CI.low.noR , ymax=Mean+CI.up.noR,fill=Type), 
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
ggsave("Output/MutFreq/MF_violinePlot.byGene_noRescaling.pdf", width = 7, height = 5)








### Plot mutation freq. across the genome based on the mutation types (Fig 1C)

pdf("Output/MutFreq/MutFreq_acrossGenome.pdf",width=14,height=5)
maxnuc=mf1$pos[nrow(mf1)]
par(mar = c(4.5,5,1,1))

plot(mfs$pos[1:maxnuc],mfs$mean[1:maxnuc],
     log="y", ylab="Mutation frequency",cex.lab=1.4,
     yaxt="n", xlab="Genome position",xaxt='n',
     col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,0.1))
axis(1,at=c(seq(500,8500,by=1000)), labels=c(seq(500,8500,by=1000)))
eaxis(side = 2, at = 10^((-1):(-(4))), cex=2)
rect(genes$start[genes$Gene=="Core"],3.2*10^-4 , genes$end[genes$Gene=="Core"], 0.1, col = "#FFC0001A", border="#FFC0001A")
rect(genes$start[genes$Gene=="HVR1"],3.2*10^-4 , genes$end[genes$Gene=="HVR1"], 0.2, col = "#FFC0001A", border="#FFC0001A")


for(i in 2:4){abline(h = 1:10 * 10^(-i), col = "gray80")}

for (i in 1:maxnuc){
    if (is.na(mfs$Type[i])==T) next
    if (mfs$Type[i]=="stop") {points(mfs$pos[i],mfs$mean[i],pch=21,col='gray30',lwd=0.3, bg="black",cex=.4)
        next}
    if (mfs$Type[i]=="syn") {c=col8[6]}
    if (mfs$Type[i]=="nonsyn"&mfs$ref[i]%in%c("c","g")) {c=col8[7]}
    if (mfs$Type[i]=="nonsyn"&mfs$ref[i]%in%c("a","t")) {c=col8[1]}
    points(mfs$pos[i],mfs$mean[i],pch=21,col='gray30',lwd=0.3, bg=paste0(c,"99"),cex=.5)
}
ylow<-0.00025
for (j in 2:(nrow(genes)-1)){
    xleft<-genes$start[j]
    xright<-genes$start[j+1]
    
    if ((j==4|j==6|j==9)){
        rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
        text(xleft+80, 1.44*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        #mtext(paste0(genes$Gene[j]),side= 1, line=-0.1, at= xleft+80, col="black", cex=0.8)
    }
    else if (j==12){
        rect(xleft,ylow,genes$end[j],1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
        text(xleft+600,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
    }
    else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
        text(xright-(xright-xleft)/2,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}

roll100<-rollmean(mfs$mean, k=100, na.rm=T, align="center")
mfs$roll100<-c(rep(NA, times=50),roll100, rep(NA, times=49))
lines(roll100~pos,data=mfs, col="#001ade", lwd=1)
abline(v=genes$end, col="gray80", lwd=.5)
#Add legend
legpos=300; legposV=0.1
rect(legpos, 0.42*legposV, (legpos+1000), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.9,pch=21,bg=col8[6],col=1,cex=1)
text((legpos+150),legposV*0.9,"Syn",adj=0, cex=1)
points((legpos+100),legposV*0.74,pch=21,bg=col8[1],col=1,cex=1)
text((legpos+150),legposV*0.74,"Nonsyn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.6,pch=21,bg=col8[7],col=1,cex=1)
text((legpos+150),legposV*0.6,"Nonsyn, C/G",adj=0, cex=1)
points((legpos+100),legposV*0.49,pch=21,bg=1,col=1,cex=1)
text((legpos+150),legposV*0.49,"Nonsense",adj=0, cex=1)

box()
dev.off()
