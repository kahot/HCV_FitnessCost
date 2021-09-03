library(tidyverse)
library(zoo)
library(purrr)
library(colorspace)
source("Rscripts/baseRscript.R")

### Test of mut freq between mutation types
#1) transition mutations
Ts<-read.csv("Output/MutFreq.filtered/Filtered.Ts.Q35.csv", row.names = 1, stringsAsFactors = F)
Ts<-Ts[Ts$pos>=342, ]

r1<-wilcox.test(Ts$mean[Ts$Type=="syn"], Ts$mean[Ts$Type=="nonsyn"], alternative = "greater", paired = FALSE) 
r2<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Ts$mean[Ts$Type=="stop"], alternative = "greater", paired = FALSE) 
r1[[3]]  #P=0
r2[[3]]  #P= 1.933446e-41
 
#2) transversion mutations
## 2.1) Summary Stats
Tv1<-read.csv("Output/MutFreq.filtered/Filtered.Tv1.MutFreq.Q35.csv", row.names = 1, stringsAsFactors = F)
Tv1<-Tv1[Tv1$pos>=342, ]
Tv2<-read.csv("Output/MutFreq.filtered/Filtered.Tv2.MutFreq.Q35.csv", row.names = 1, stringsAsFactors = F)
Tv2<-Tv2[Tv2$pos>=342, ]

Syn<-c(Tv1$mean[Tv1$Type.tv1=="syn"],Tv2$mean[Tv2$Type.tv2=="syn"])
Nonsyn <- c(Tv1$mean[Tv1$Type.tv1=="nonsyn"],Tv2$mean[Tv2$Type.tv2=="nonsyn"])
Stop <- c(Tv1$mean[Tv1$Type.tv1=="stop"],Tv2$mean[Tv2$Type.tv2=="stop"])

r1<-wilcox.test(Syn, Nonsyn, alternative = "greater", paired = FALSE) 
r2<-wilcox.test(Nonsyn, Stop, alternative = "greater", paired = FALSE) 
wilcox.test(Nonsyn, Stop, alternative = "less", paired = FALSE) 

r1[[3]] #2.788528e-175
r2[[3]] #1

### Transition vs. transversion

tvs<-c(Tv1$mean, Tv2$mean)

r3<-wilcox.test(Ts$mean, tvs, alternative = "greater", paired = FALSE) 
r3[[3]]
#W = 124680000, p-value < 2.2e-16 (P=0)

r4<-wilcox.test(Ts$mean[Ts$Type=="syn"], Syn, alternative = "greater", paired = FALSE) 
r4[[3]]
#W = 7866400, p-value < 2.2e-16 (P=0)
r5<-wilcox.test(Ts$mean[Ts$Type=="nonsyn"], Nonsyn, alternative = "greater", paired = FALSE) 
r5[[3]]
#W = 62349000, p-value < 2.2e-16 (P=0)



##############
# Transition mutations: test by nucleotide and by gene
Ts2<-Ts[Ts$makesCpG==0,]

depth<-read.csv("Output/ReadDepth_sum.csv",stringsAsFactors = F, row.names = 1)



## Wilcox test by NT
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))

for (f in 1:2){
        if (f==1) {dat<-Ts;  fname <- "" }
        if (f==2) {dat<-Ts2; fname<-"_noCpG"}
        
        WilcoxTest.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt)<-c("NT1","NT2","test","P.value")
        
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$mean[dat$ref==Ncomb[i,1]]
                vec2<-dat$mean[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
                
                WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt$test[i]<-"less"
                WilcoxTest.nt$P.value[i]<-result[[3]]
        }   
        
        WilcoxTest.nt2<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","P.value")
        
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$mean[dat$ref==Ncomb[i,1]]
                vec2<-dat$mean[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
                
                WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt2$test[i]<-"greater"
                WilcoxTest.nt2$P.value[i]<-result[[3]]
        }   
        
        WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
        write.csv(WilcoxTest.nt, paste0("Output/SummaryStats/MF_WilcoxTestResults_byNT", fname,".csv"))
}



## Wilcoxon test by gene
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[2:12]
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
    gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector


#run Wilcoxin Test  
Gcomb<-t(combn(genenames,2))
WilcoxTest.gene<-data.frame(matrix(ncol=3,nrow=11))
colnames(WilcoxTest.gene)<-c("gene","test","P.value")


mf1<-merge(Ts,depth, by="pos", all.x=T)
mf1<-mf1[!is.na(mf1$mean),]
mf1<-merge(mf1, genetable, by="pos", all.x=T )

for (i in 1:11) {
    vec1<-mf1$mean[mf1$gene==genenames[1] & mf1$Type=="syn"]
    vec2<-mf1$mean[mf1$gene==genenames[i] & mf1$Type=="nonsyn"]
    result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
    
    WilcoxTest.gene$gene[i]<-genenames[(i+1)]
    WilcoxTest.gene$test[i]<-"greater"
    WilcoxTest.gene$P.value[i]<-result[[3]]
}   

write.csv(WilcoxTest.gene,"Output/SummaryStats/MF_WilcoxTestResults_byGene.synvsNonsyn.csv")











