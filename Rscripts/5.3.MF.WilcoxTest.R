library(tidyverse)
library(zoo)
library(purrr)
library(colorspace)
source("Rscripts/baseRscript.R")


Ts<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv", row.names = 1, stringsAsFactors = F)
Ts<-Ts[Ts$pos>=342, ]
Ts2<-Ts[Ts$makesCpG==0,]

depth<-read.csv("Output1A/ReadDepth_sum.csv",stringsAsFactors = F, row.names = 1)



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
        write.csv(WilcoxTest.nt, paste0("Output1A/SummaryStats/MF_WilcoxTestResults_byNT", fname,".csv"))
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

write.csv(WilcoxTest.gene,"Output1A/SummaryStats/MF_WilcoxTestResults_byGene.synvsNonsyn.csv")











