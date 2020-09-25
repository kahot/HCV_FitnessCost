library(plotrix)
library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(colorspace)
library(ggplot2)
colors2<-qualitative_hcl(6, palette="Dark3")
library(dplyr)
source("Rscripts/baseRscript.R")


#Create the filtered mut frequency table of all samples 

HCVFiles_overview3<-list.files("Output1A/Overview3/",pattern="overview3.csv")
FilteredOverview2<-list()
for (i in 1:length(HCVFiles_overview3)){ 
        overviews2<-read.csv(paste0("Output1A/Overview3/",HCVFiles_overview3[i]),stringsAsFactors=FALSE, row.names=1)
        FilteredOverview2[[i]]<-overviews2
        names(FilteredOverview2)[i]<-substr(paste(HCVFiles_overview3[i]),start=1,stop=7)
}

##################################
MutFreq_Ts<-list()
MutFreq_tv1<-list()
MutFreq_tv2<-list()
MutFreq_tvs<-list()
MutFreq_all<-list()

for (i in 1:length(FilteredOverview2)){
        dat<-FilteredOverview2[[i]]
        filename<-names(FilteredOverview2)[i]
        
        MutFreq_Ts[[i]]<-dat[,c("pos","freq.Ts.ref")] 
        MutFreq_tv1[[i]]<-dat[,c("pos","freq.transv1.ref")] 
        MutFreq_tv2[[i]]<-dat[,c("pos","freq.transv2.ref")] 
        MutFreq_tvs[[i]]<-dat[,c("pos","freq.transv.ref")] 
        MutFreq_all[[i]]<-dat[,c("pos","freq.mutations.ref")] 
        
        names(MutFreq_Ts)[i]<-filename
        names(MutFreq_tv1)[i]<-filename
        names(MutFreq_tv2)[i]<-filename
        names(MutFreq_tvs)[i]<-filename
        names(MutFreq_all)[i]<-filename
        
        
}
#assign column names for the list
for (i in 1:length(MutFreq_Ts)) {
        colnames(MutFreq_Ts[[i]])<-c("pos",paste0(names(MutFreq_Ts[i])))
        colnames(MutFreq_tv1[[i]])<-c("pos",paste0(names(MutFreq_tv1[i])))
        colnames(MutFreq_tv2[[i]])<-c("pos",paste0(names(MutFreq_tv2[i])))
        colnames(MutFreq_tvs[[i]])<-c("pos",paste0(names(MutFreq_tvs[i])))
        colnames(MutFreq_all[[i]])<-c("pos",paste0(names(MutFreq_all[i])))
}

Ts<-MutFreq_Ts%>% purrr::reduce(full_join, by='pos')
Tv1.MutFreq<-MutFreq_tv1 %>% purrr::reduce(full_join, by='pos')
Tv2.MutFreq<-MutFreq_tv2 %>% purrr::reduce(full_join, by='pos')
Tvs.MutFreq<-MutFreq_tvs %>% purrr::reduce(full_join, by='pos')
AllMutFreq<-MutFreq_all %>% purrr::reduce(full_join, by='pos')


### all mut. freq with metadata ###

files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )
cnames<-c("",".tv1",".tv2",".tvs","")
mf.files<-list()
s<-length(FilteredOverview2)
M<-FilteredOverview2[[3]]
colnames(M)[33:35]<-c("MutAA","MutAA.tv1","MutAA.tv2")

for (i in 1:5) {
        dat<-get(files[i])
        dat$mean<-rowMeans(dat[2:(s+1)],na.rm=T, dims=)
        if (i==1|i==5) muttypes<-M[,c("pos","ref", "Type","WTAA","MutAA","makesCpG","bigAAChange")]
        if(i==2|i==3) muttypes<-M[,c("pos","ref", paste0("Type",cnames[i]),"WTAA",paste0("MutAA",cnames[i]), paste0("makesCpG",cnames[i]), paste0("bigAAChange",cnames[i]))]
        if(i==4) muttypes<-M[,c("pos","ref","Type.tv1","Type.tv2","WTAA","MutAA.tv1","MutAA.tv2","makesCpG.tv1","makesCpG.tv1","bigAAChange.tv1","bigAAChange.tv2")]
        
        dat2<-merge(dat,muttypes,by="pos")
        mf.files[[i]]<-dat2
        names(mf.files)[i]<-files[i]
        write.csv(dat2,paste0("Output1A/MutFreq.filtered/Filtered.",files[i],".Q35.csv"))
}



##############################
##############################

#Create the mut freq summary table by gene and by type:
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

#read in files
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq")
mf.files<-list()
for (i in 1:5){
        data<-read.csv(paste0("Output1A/MutFreq.filtered/Filtered.", files[i],".Q35.csv"),row.names = 1,stringsAsFactors = F)
        assign(paste0(files[i]),data)
        mf.files[[i]]<-data
        names(mf.files)[i]<-files[i]
}

#Summarize the mean and se for all types of mutations
# 1.Transitions:
Ts<-mf.files[[1]]
Ts<-Ts[Ts$pos>=342,]

#calculate SE with read depth from all files:
#create a corresponding matrix of read deapth for each Ts site

#Get the file names (SeqData files)
SeqDt<-list.files("Output1A/SeqDataQ35/",pattern="SeqData")

Reads.Ts<-list()
for (i in 1:length(SeqDt)){   
        #for (i in 1:1){
        id<-substr(paste(SeqDt[i]),start=9,stop=15)
        print(id)
        DF<-read.csv(paste0("Output1A/SeqDataQ35/",SeqDt[i]),stringsAsFactors=FALSE, row.names = 1)
        DF<-DF[DF$pos>=342 & DF$pos<=8575,]
        DF$TsReads<-NA
        #determine the read depth of Transition mutations at every site.
        for (k in 1:nrow(DF)){
                if (is.na(DF$MajNt[k])) next
                else if (DF$MajNt[k]!=DF$ref[k]) next
                else {DF$TsReads[k]<-DF[k,paste0(DF$ref[k])]+DF[k,paste0(DF$transition.ref[k])]}
        }
        
        Reads.Ts[[i]]<-DF[,c("pos","TsReads")]
        names(Reads.Ts)[i]<-id
}

for (i in 1:length(Reads.Ts)) {
        colnames(Reads.Ts[[i]])<-c("pos",paste0(names(Reads.Ts[i])))
}

Reads_Ts<-Reads.Ts%>% purrr::reduce(full_join, by='pos')
write.csv(Reads_Ts, "Output1A/ReadDepth_Transitions.csv")

## Total Reads
Reads<-list()
for (i in 1:length(SeqDt)){   
        #for (i in 1:1){
        id<-substr(paste(SeqDt[i]),start=9,stop=15)
        print(id)
        DF<-read.csv(paste0("Output1A/SeqDataQ35/",SeqDt[i]),stringsAsFactors=FALSE, row.names = 1)
        DF<-DF[DF$pos>=342 & DF$pos<=8575,]
        
        Reads[[i]]<-DF[,c("pos","TotalReads")]
        names(Reads)[i]<-id
}

for (i in 1:length(Reads)) {
        colnames(Reads[[i]])<-c("pos",paste0(names(Reads[i])))
}

Reads_Total<-Reads%>% purrr::reduce(full_join, by='pos')
write.csv(Reads_Total, "Output1A/ReadDepth_All.csv")

Reads_Total$Depth<-rowSums(Reads_Total[2:196], na.rm=T)


# Use general read deapth for now:
depth<-Reads_Total[,c("pos","Depth")]
write.csv(depth,"Output1A/ReadDepth_sum.csv")

#### Start from here if previously calculated 'depth'
depth<-read.csv("Output1A/ReadDepth_sum.csv",stringsAsFactors = F, row.names = 1)

tb3<-data.frame(type=c("Ts","Ts.syn","Ts.ns","Ts.stop", "Tv1","Tv1.syn","Tv1.ns","Tv1.stop","Tv2","Tv2.syn","Tv2.ns","Tv2.stop","Tvs","All" ))
files<-c("Ts","Tv1.MutFreq","Tv2.MutFreq","Tvs.MutFreq","AllMutFreq" )
## Add total read depth for each site to calculate SE
for (i in 1:length(files)){
        dt<-mf.files[[i]]
        #coding region only
        dt<-dt[dt$pos>341,]
        dt<-merge(dt,depth, by="pos", all.x = T)
        if (i<=3){
                colnames(dt)[199]<-"Type"
                k<-(i-1)*4+1
                tb3$mean[k]<-mean(dt$mean,na.rm=T)
                tb3$se[k]<-sqrt(tb3$mean[k]*(1-tb3$mean[k])/mean(dt$Depth)) 
                k<-k+1
                tb3$mean[k]<-mean(dt$mean[dt$Type=="syn"],na.rm=T)
                tb3$se[k]<-sqrt(tb3$mean[k]*(1-tb3$mean[k])/mean(dt$Depth[dt$Type=="syn"])) 
                k<-k+1
                tb3$mean[k]<-mean(dt$mean[dt$Type=="nonsyn"],na.rm=T)
                tb3$se[k]<-sqrt(tb3$mean[k]*(1-tb3$mean[k])/mean(dt$Depth[dt$Type=="nonsyn"])) 
                k<-k+1
                tb3$mean[k]<-mean(dt$mean[dt$Type=="stop"],na.rm=T)
                tb3$se[k]<-sqrt(tb3$mean[k]*(1-tb3$mean[k])/mean(dt$Depth[dt$Type=="stop"])) 
        }
        #Tvs
        if (i==4) {k=13
        tb3$mean[k]<-mean(dt$mean,na.rm=T)
        tb3$se[k]<-sqrt(tb3$mean[k]*(1-tb3$mean[k])/mean(dt$Depth)) }
        if (i==5) {k=14
        tb3$mean[k]<-mean(dt$mean,na.rm=T)
        tb3$se[k]<-sqrt(tb3$mean[k]*(1-tb3$mean[k])/mean(dt$Depth)) }
}

tb3$type<-as.character(tb3$type)

addtb3<-data.frame(type=c("tvs.syn","tvs.ns","tvs.stop", "all.syn","all.ns","all.stop"),
                   mean=c(mean(tb3$mean[6],tb3$mean[10]),mean(tb3$mean[7],tb3$mean[11]),mean(tb3$mean[8],tb3$mean[12]),
                          mean(tb3$mean[2],tb3$mean[6],tb3$mean[10]),mean(tb3$mean[3],tb3$mean[7],tb3$mean[11]),mean(tb3$mean[4],tb3$mean[8],tb3$mean[12])),
                   se  =c(mean(tb3$se[6],tb3$se[10]),mean(tb3$se[7],tb3$se[11]),mean(tb3$se[8],tb3$se[12]),
                          mean(tb3$se[2],tb3$se[6],tb3$se[10]),mean(tb3$se[3],tb3$se[7],tb3$se[11]),mean(tb3$se[4],tb3$se[8],tb3$se[12])))

TB3<-rbind(tb3, addtb3)

write.csv(TB3, "Output1A/MutFreq.filtered/MF.Mean.SE.summary_updated.csv")








######## Plot mut freq by Gene and by type
##
mf1<-mfs[!is.na(mfs$mean),]
mf1<-merge(mf1, genetable, by="pos", all.x=T )

#exclude stop mutations
mf1<-mf1[mf1$Type!="stop",]

SumMFGenes<-aggregate(mf1$mean,by=list(mf1$gene, mf1$Type),FUN=mean)
colnames(SumMFGenes)<-c("Gene","Type","Mean")

#scDF<-sc[,c("gene","EstSC","Depth")]
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
write.csv(sumG2, "Output1A/MutFreq.filtered/MF_Summary_Table.by.gene.byType.csv")


ggplot(sumG2, aes(x=Gene, y=Mean, group=Type, color=Type))+
        geom_point(position=position_dodge(width=0.3))+scale_color_manual(values=colors2[c(1,5)])+
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.4, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Mutation frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(axis.title.x=element_blank())
ggsave(filename="Output1A/SummaryFig.Filtered/Ave.MF_by.gene_by.type_updated.pdf", width = 8.5, height = 5)


###make it into one figure

sumG$Type<-"Mean"
sumG<-sumG[,c(1,4,2,3)]
sumG2<-rbind(sumG2,sumG)

ggplot(sumG2, aes(x=Gene, y=Mean))+
        scale_y_continuous(trans="log10",breaks = ybreaks,labels=ylabel)+
        geom_point(data=mf2, aes(x=gene, y=mean), position=position_jitter(width=0.2),stat = "identity", col=paste0(cold[2],"4D"), size=0.2)+
        geom_point(data=sumG2, aes(x=Gene, y=Mean, group=Type, color=Type), position=position_dodge(width=0.6))+
        geom_errorbar(data=sumG2, aes(ymin=Mean-SE, ymax=Mean+SE,color=Type), width=.6, position=position_dodge(width=0.6))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+
        ylab("Mutation frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13), panel.grid.minor.y = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  
                   color = "gray70", size=.5)+
        theme(axis.title.x=element_blank())+
        scale_color_manual(values=c("blue3",colors2[c(1,5)]), labels=c("Mean","Nonsyn mean","Syn mean"))+
        theme(legend.title=element_blank())+
        guides(fill = guide_legend(override.aes = list(linetype = 0)),
               color = guide_legend(override.aes = list(linetype = 0)))+
        annotate("text", x=c(1,2,4:11), y=rep(0.035, times=10), label="*",size=4)

ggsave("Output1A/SummaryFig.Filtered/MutFreq.byGene.byType.pdf", width=8, height = 4)       













##############################

### Average read depth 
Reads<-data.frame(Seq=names(FilteredOverview2))

for (i in 1:length(FilteredOverview2)){
        dat<-FilteredOverview2[[i]]
        dat<-dat[!is.na(dat$freq.Ts),]
        Reads$ave[i]<-mean(dat$TotalReads, na.rm=T) 
        Reads$median[i]<-median(dat$TotalReads, na.rm = T)
        Reads$max[i]<-max(dat$TotalReads, na.rm=T)
        Reads$min[i]<-min(dat$TotalReads, na.rm=T)
        
}

mean(Reads$ave)
mean(Reads$median) #5282.449
mean(Reads$max) #22005.66

sum(Reads$ave>=5000)

### Structural vs. non-structural genes mut. freq
TS<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)

st<-TS[TS$pos>341& TS$pos<=2579,]
nonst<-TS[TS$pos<=2579,]

mean(st$mean, na.rm=T) #0.004964672
mean(nonst$mean, na.rm = T)  #0.004872566

r1<-wilcox.test(st$mean,nonst$mean, alternative = "greater", paired = FALSE) 
r1[[3]]  #P=0.1054109 Not Significant

######
######
#Calculate CI
#1. use proportion SE and using z-score

Ts$n<-apply(Ts[,2:196], 1, FUN=function(x) sum(!is.na(x)))


Ts$SE<-sqrt(Ts[,"mean"]*(1-Ts[,"mean"])/Ts[,"n"]) 
Ts$CI_up<-Ts[,"mean"]+Ts[,"SE"]*1.96
Ts$CI_down<-Ts[,"mean"]-Ts[,"SE"]*1.96
Ts$CI<-Ts[,"SE"]*1.96

Ts2<-Ts[,c("pos","mean","SE","CI", "CI_up","CI_down")]
write.csv(Ts2,"Output1A/MutFreq.filtered/Ts_C_summary.csv")
ggplot(Ts2, aes(x=pos, y=mean))+
        geom_errorbar(aes(ymin=mean-CI, ymax=mean+CI), width=.2, size=.2, color="gray60")+
        geom_point(size =.3, color="blue")+

        theme(axis.title.x=element_blank())+ylab("Mutation frequency ± 95% CI")+
        theme_bw()+
        labs(x="")

        

ggplot(Ts2)+
        geom_point(aes(x=pos, y=mean), color="blue", size=.3)+
        geom_point(aes(x=pos, y=CI_up), size=.5, color="gray30", shape=45)+
        theme_bw()+
        theme(axis.title.x=element_blank())+ylab("Mutation frequency ± 95% CI")+
        labs(x="")


