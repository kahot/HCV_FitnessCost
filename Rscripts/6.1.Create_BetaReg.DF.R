#Beta regression preparation:

library(DataCombine)
library(miscTools)
library(caret)

#### Covert the information to the BetaReg data format

dat<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",row.names=1)
dat<-dat[dat$pos>=342,c("pos","ref","mean","Type","makesCpG", "bigAAChange")]

dat1<-dat[,c("mean","ref","Type")]

new<-dummyVars("~.", data=dat1)
brDF<-data.frame(predict(new, newdata = dat1))
colnames(brDF)[2:7]<-c("a","c","g","t","Nonsyn","Stop")
brDF<-cbind(brDF[,1:7], dat[,c("pos","makesCpG", "bigAAChange")])

colnames(brDF)[9]<-"CpG"


### add gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"

gene<-c()
for (i in 2:12){
        gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)

#add gene info as 2_12
brDF<-merge(brDF, g, by ="pos")

#add as one hot encoindg
for (i in 2:12){
        gname<-paste(genes$Gene[i])
        n<-ncol(brDF)
        brDF[,n+1]<-0
        colnames(brDF)[n+1]<-gname
        brDF[brDF$gene==i,n+1]<-1
}

# Add RNA structural information
rna<-read.csv("Data/RNAStructure_conserved.csv", stringsAsFactors = F)
rnaList<-list()
for (i in 1:nrow(rna)){
        pos<-rna$Start[i]:rna$End[i]
        Shape<-rep(1, times=rna$End[i]-rna$Start[i]+1)
        rnaList[[i]]<-cbind(pos, Shape)
}

RnaShape<-data.frame(do.call(rbind,rnaList))
brDF<-merge(brDF,RnaShape, by="pos", all.x = T)
brDF$Shape[is.na(brDF$Shape)]<-0

write.csv(brDF, "Output/BetaReg/BetaRegData.csv")



