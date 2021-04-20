#GLM/Beta regression preparation:
# RUnnig only the filtered Q35 Ts data (no Ts_zero and Ts-NA)

library(tidyverse)
library(zoo)
library(purrr)
library(DataCombine)
library(miscTools)


source("Rscripts/baseRscript.R")

## Run 10.2.Filterout_lowfrew.R first ##

#### For preparing BetaReg data formatting

nucord <- c("a", "t", "c", "g")
brform<-as.matrix(data.frame(a="",t="",c="",g="", Syn ="",Nonsyn="",Stop=""))
nuc<-as.matrix(data.frame(a="",t="",c="",g=""))

dat<-read.csv("Output/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names=1)
dat<-dat[dat$pos>=342,]
dat1<-dat[,c("pos","makesCpG","bigAAChange","mean")]

for (i in 1:nrow(dat)){
        atcg <- c(0,0,0,0)
        atcg[which(nucord == dat[i,]$ref)] <- 1
        nuc<-miscTools::insertRow(nuc,i,atcg)
        nonsyn <- as.numeric(regexpr("nonsyn",dat[i,]$Type) > 0)
        stop <- as.numeric(regexpr("stop",dat[i,]$Type) > 0)
        syn<-as.numeric(regexpr("^syn",dat[i,]$Type) > 0)
        new<-c(atcg,syn,nonsyn,stop)
        brform<-insertRow(brform,i,new)
}
BetaRegData<-cbind(dat1$pos,brform[1:nrow(dat1),])
BetaRegData<-cbind(BetaRegData,dat1[,2:4])
colnames(BetaRegData)[9]<-"CpG"
#write.csv(BetaRegData, "Output/BetaReg/BetaRegData_Ts_FilteredData.csv")


### addd gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene
gene<-c()
for (i in 2:12){
        gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)


colnames(BetaRegData)[1]<-"pos"
betar<-merge(BetaRegData, g, by ="pos")

for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(betar)
        betar[,n+1]<-0
        colnames(betar)[n+1]<-gname
        betar[betar$gene==i,n+1]<-1
}
co<-which(colnames(betar)=="NS1(P7)" )
colnames(betar)[co]<-"NS1"

write.csv(betar,"Output/BetaReg/BetaRegFull.Ts.FilteredData.csv")


## Add RNA structural information
DF<-read.csv("Output/BetaReg/BetaRegFull.Ts.FilteredData.csv", stringsAsFactors = F, row.names = 1)
colnames(DF)[1]<-"pos"
rna<-read.csv("Data/RNAStructure_Conserved2.csv", stringsAsFactors = F)
rnaList<-list()
for (i in 1:nrow(rna)){
        pos<-rna$Start[i]:rna$End[i]
        Shape<-rep(1, times=rna$End[i]-rna$Start[i]+1)
        rnaList[[i]]<-cbind(pos, Shape)
}

RnaShape<-data.frame(do.call(rbind,rnaList))
DF<-merge(DF,RnaShape, by="pos", all.x = T)
DF$Shape[is.na(DF$Shape)]<-0

write.csv(DF, "Output/BetaReg/BetaReg.Data.Shape.csv")




