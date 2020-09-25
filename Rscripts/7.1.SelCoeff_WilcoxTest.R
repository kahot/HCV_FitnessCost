library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
source("Rscripts/baseRscript.R")

#SC<-list()
#fnames<-c("Ts", "Ts_NA", "Ts_zero")
mutrates<-read.csv("Output1A/Geller/Geller.MutRates.Summary_updated.csv", row.names = 1, stringsAsFactors = F)


df<-read.csv("Output1A/Mutfreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names = 1 )
df<-df[,c(1,197:203)]
        
df$TSmutrate[df$ref=="a"]<-mutrates$MutRate[mutrates$Mutation=="AG"]
df$TSmutrate[df$ref=="c"]<-mutrates$MutRate[mutrates$Mutation=="CU"]
df$TSmutrate[df$ref=="g"]<-mutrates$MutRate[mutrates$Mutation=="GA"]
df$TSmutrate[df$ref=="t"]<-mutrates$MutRate[mutrates$Mutation=="UC"]
df$EstSC<-""
for (j in 1:nrow(df)){
        df$EstSC[j] <- EstimatedS(df$TSmutrate[j],df$mean[j])
}
df$EstSC<-as.numeric(df$EstSC)
depth<-read.csv("Output1A/ReadDepth_sum.csv", stringsAsFactors = F, row.names = 1)
df<-merge(df, depth, by="pos")

write.csv(df,"Output1A/SelCoeff/SC.csv")

df<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
### Wilcoxon Test on SC  (use mean SC)
## Use A & T only for CpG Analysis
ty<-which(colnames(df)=="Type")
dat<-df
#fname=names(SC)[i]
k=1
TypeList<-list()
for (typeofsite in c("syn", "nonsyn","stop")){
        all<-dat$EstSC[dat[,ty]==typeofsite]
        allcpg<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==1]
        allnoncpg1<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==0 & dat$ref=="t"]
        allnoncpg2<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==0 & dat$ref=="a"]
        for (wtnt in c("a", "t", "c", "g")){
                selco<- dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt]
                sc_NonCpG<-dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt & dat$makesCpG==0]
                sc_CpG<-dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt & dat$makesCpG==1]
                vectorname<-paste0(typeofsite,"_",wtnt)
                assign(vectorname, selco)
                vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
                assign(vname1, sc_NonCpG)
                vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
                assign(vname2, sc_CpG)
                }
        typev<-paste0(typeofsite,"_all")
        assign(typev, all)
        TypeList[[k]]<-all
        names(TypeList)[k]<-typev
        k=k+1
        cpgv1<-paste0(typeofsite,"_allCpG")
        assign(cpgv1, allcpg)
        TypeList[[k]]<-allcpg
        names(TypeList)[k]<-cpgv1
        k=k+1
        cpgv0<-paste0(typeofsite,"_allnonCpG")
        allnoncpg<-c(allnoncpg1,allnoncpg2)
        assign(cpgv0, allnoncpg)
        TypeList[[k]]<-allnoncpg
        names(TypeList)[k]<-cpgv0
        k=k+1
        }

#wilcox.test 
wilcoxtest1<-data.frame("test"=matrix(nrow=4))
        
re1<-wilcox.test(syn_all,nonsyn_all, alternative = "less", paired = FALSE) 
wilcoxtest1$test[1]<-re1[[7]]
wilcoxtest1$P.value[1]<-re1[[3]]
        
re2<-wilcox.test(nonsyn_all,stop_all, alternative = "less", paired = FALSE) 
wilcoxtest1$test[2]<-re2[[7]]
wilcoxtest1$P.value[2]<-re2[[3]]
        
re3<-wilcox.test(syn_allCpG,syn_allnonCpG, alternative = "greater", paired = FALSE) 
wilcoxtest1$test[3]<-re3[[7]]
wilcoxtest1$P.value[3]<-re3[[3]]
        
re4<-wilcox.test(nonsyn_allCpG,nonsyn_allnonCpG, alternative = "greater", paired = FALSE) 
wilcoxtest1$test[4]<-re4[[7]]
wilcoxtest1$P.value[4]<-re4[[3]]
        
write.csv(wilcoxtest1,"Output1A/SelCoeff/WilcoxonResults_by_Type.csv")

Type.sc<-data.frame("mean"=matrix(nrow=9))
for (i in 1:9){
        rownames(Type.sc)[i]<-names(TypeList)[i]
        Type.sc$mean[i]<-mean(TypeList[[i]],na.rm=T)
        Type.sc$se[i]<-std.error(TypeList[[i]],na.rm=T)
}
write.csv(Type.sc,"Output1A/SelCoeff/SC_byType_Summary.csv")
      
  
## Test on NT by NT
S<-data.frame()
se<-data.frame()
S_CpG<-data.frame()
se_CpG<-data.frame()
S_nonCpG<-data.frame()
se_nonCpG<-data.frame()
for (typeofsite in c("syn", "nonsyn","stop")){
     for (wtnt in c("a", "t", "c", "g")){
        sc<- dat$EstSC[dat[,ty]==typeofsite & dat$ref==wtnt]
        S[typeofsite,wtnt]<-mean(sc[!is.na(sc)])
        se[typeofsite,wtnt]<-std.error(sc[!is.na(sc)])
        
        S_NonCpG<-dat$EstSC[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==0]
        S_nonCpG[typeofsite,wtnt]<-mean(S_NonCpG[!is.na(S_NonCpG)])
        se_nonCpG[typeofsite,wtnt]<-std.error(S_NonCpG[!is.na(S_NonCpG)])
        
        Sc_CpG<-dat$EstSC[dat$Type==typeofsite & dat$ref==wtnt & dat$makesCpG==1]
        S_CpG[typeofsite,wtnt]<-mean(Sc_CpG[!is.na(Sc_CpG)])
        se_CpG[typeofsite,wtnt]<-std.error(Sc_CpG[!is.na(Sc_CpG)])
        
        vectorname<-paste0(typeofsite,"_",wtnt)
        assign(vectorname, sc)
        vname1<<-paste0(typeofsite,"_",wtnt,"_noncpg")
        assign(vname1, S_NonCpG)
        vname2<<-paste0(typeofsite,"_",wtnt,"_cpg")
        assign(vname2, Sc_CpG)
     }
}
rownames(se)<-c("syn_se","nonsyn_se","stop_se")
rownames(S_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(se_nonCpG)<-c("syn_noncpg_se","nonsyn_noncpg_se","stop_noncpg_se")
rownames(S_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")
rownames(se_CpG)<-c("syn_cpg_se","nonsyn_cpg_se","stop_cpg_se")

SCbyType<-rbind(S,se,S_nonCpG,se_nonCpG,S_CpG,se_CpG)
SCbyType2<-data.frame(t(SCbyType))
write.csv(SCbyType2,"Output1A/SelCoeff/SC_byNt_byType.csv")

# Wilcoxin Test by nucleotide 
nuc.sc<-data.frame("syn.cpg"=matrix(nrow=4))
rownames(nuc.sc)<-c("a","t","c","g")


WilcoxTest.nt.sc<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.nt.sc)<-c("nt","test","P.value")
k=1
for (i in c("a","t","c","g")) {
    if (i=="a"|i=="t"){
        syncpg<-get(paste0("syn_",i,"_cpg"))
        synnoncpg<-get(paste0("syn_",i,"_noncpg"))
        nonsyncpg<-get(paste0("nonsyn_",i,"_cpg"))
        nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
        
        nuc.sc$syn.cpg[k]<-mean(syncpg, na.rm=T)
        nuc.sc$syn.cpg.se[k]<-std.error(syncpg, na.rm=T)
        nuc.sc$syn.ncpg[k]<-mean(synnoncpg, na.rm=T)
        nuc.sc$syn.ncpg.se[k]<-std.error(synnoncpg, na.rm=T)
        nuc.sc$ns.cpg[k]<-mean(nonsyncpg, na.rm=T)
        nuc.sc$ns.cpg.se[k]<-std.error(nonsyncpg, na.rm=T)
        nuc.sc$ns.ncpg[k]<-mean(nonsynnoncpg, na.rm=T)
        nuc.sc$ns.ncpg.se[k]<-std.error(nonsynnoncpg, na.rm=T)
        nuc.sc$stop.ncpg[k]<-NA
        nuc.sc$stop.ncpg.se[k]<-NA
        
        if (i=="a"){
                result1<-wilcox.test(syncpg, synnoncpg, alternative = "greater", paired = FALSE) 
                result2<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "greater", paired = FALSE)
                for (r in 1:2){
                        result<-get(paste0('result',r))
                        WilcoxTest.nt.sc$nt[r]<-i
                        WilcoxTest.nt.sc$test[r]<-result[[7]]
                        WilcoxTest.nt.sc$P.value[r]<-result[[3]]}
        }
        if (i=="t"){
                result3<-wilcox.test(syncpg, synnoncpg, alternative = "greater", paired = FALSE) 
                result4<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "greater", paired = FALSE)
                for (r in 3:4){
                        result<-get(paste0('result',r))
                        WilcoxTest.nt.sc$nt[r]<-i
                        WilcoxTest.nt.sc$test[r]<-result[[7]]
                        WilcoxTest.nt.sc$P.value[r]<-result[[3]]}
        }
    }
else {  synnoncpg<-get(paste0("syn_",i,"_noncpg"))
        nonsynnoncpg<-get(paste0("nonsyn_",i,"_noncpg"))
        stopnoncpg<-get(paste0("stop_",i,"_noncpg"))
        
        nuc.sc$syn.cpg[k]<-NA
        nuc.sc$syn.cpg.se[k]<-NA
        nuc.sc$syn.ncpg[k]<-mean(synnoncpg, na.rm=T)
        nuc.sc$syn.ncpg.se[k]<-std.error(synnoncpg, na.rm=T)
        nuc.sc$ns.cpg[k]<-NA
        nuc.sc$ns.cpg.se[k]<-NA
        nuc.sc$ns.ncpg[k]<-mean(nonsynnoncpg, na.rm=T)
        nuc.sc$ns.ncpg.se[k]<-std.error(nonsynnoncpg, na.rm=T)
        nuc.sc$stop.ncpg[k]<-mean(stopnoncpg, na.rm=T)
        nuc.sc$stop.ncpg.se[k]<-std.error(stopnoncpg, na.rm=T)
        
        if (i =="c") {
                result5<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "less", paired = FALSE) 
                WilcoxTest.nt.sc$nt[5]<-i
                WilcoxTest.nt.sc$test[5]<-result5[[7]]
                WilcoxTest.nt.sc$P.value[5]<-result5[[3]]     }           
        if (i =="g") {
                result6<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "less", paired = FALSE) 
                WilcoxTest.nt.sc$nt[6]<-i
                WilcoxTest.nt.sc$test[6]<-result6[[7]]
                WilcoxTest.nt.sc$P.value[6]<-result6[[3]]    }
        }
     k=k+1
} 
write.csv(WilcoxTest.nt.sc,"Output1A/SelCoeff/WilcoxTestResults_byNT.csv")
write.csv(nuc.sc,"Output1A/SelCoeff/SC_Summary_byNT.csv")

#####

## Wilcoxon test on SCs by gene 
df<-read.csv("Output1A/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
#coding regions only
df<-df[df$pos>=342,]
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[1:12]
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-df$pos[nrow(df)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]
sc<-merge(df, genetable, by="pos")

Gcomb<-t(combn(genenames[2:12],2))
WilcoxTest2.gene<-data.frame(matrix(ncol=4,nrow=nrow(Gcomb)))
colnames(WilcoxTest2.gene)<-c("gene1","gene2","test","P.value")

for (i in 1:nrow(Gcomb)) {
        vec1<-sc$EstSC[sc$gene==Gcomb[i,1]]
        vec2<-sc$EstSC[sc$gene==Gcomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest2.gene$gene1[i]<-Gcomb[i,1]
        WilcoxTest2.gene$gene2[i]<-Gcomb[i,2]
        WilcoxTest2.gene$test[i]<-"less"
        WilcoxTest2.gene$P.value[i]<-result[[3]]
}   

WilcoxTest2.gene2<-data.frame(matrix(ncol=4,nrow=nrow(Gcomb)))
colnames(WilcoxTest2.gene2)<-c("gene1","gene2","test","P.value")

for (i in 1:nrow(Gcomb)) {
        vec1<-sc$EstSC[sc$gene==Gcomb[i,1]]
        vec2<-sc$EstSC[sc$gene==Gcomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
        
        WilcoxTest2.gene2$gene1[i]<-Gcomb[i,1]
        WilcoxTest2.gene2$gene2[i]<-Gcomb[i,2]
        WilcoxTest2.gene2$test[i]<-"greater"
        WilcoxTest2.gene2$P.value[i]<-result[[3]]
}        

WilcoxTest2.gene<-rbind(WilcoxTest2.gene,WilcoxTest2.gene2)
write.csv(WilcoxTest2.gene,"Output1A/SelCoeff/SC_WilcoxTestResults_byGene.csv")


#####
## wilcox test of SC by NT
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))
WilcoxTest.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt)<-c("NT1","NT2","test","P.value")

for (i in 1:nrow(Ncomb)) {
        vec1<-sc$EstSC[sc$ref==Ncomb[i,1]]
        vec2<-sc$EstSC[sc$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 

        WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt$test[i]<-"greater"
        WilcoxTest.nt$P.value[i]<-result[[3]]
}   

WilcoxTest.nt2<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","P.value")

for (i in 1:nrow(Ncomb)) {
        vec1<-sc$EstSC[sc$ref==Ncomb[i,1]]
        vec2<-sc$EstSC[sc$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt2$test[i]<-"less"
        WilcoxTest.nt2$P.value[i]<-result[[3]]
}   

WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
write.csv(WilcoxTest.nt,"Output1A/SelCoeff/SC_WilcoxTestResults_byNT.csv")


####

sc3<-sc[sc$makesCpG==0,]

## wilcox test by NT -with or without CpG sites
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))

for (f in 1:2){
        if (f==1) dat<-sc; fname<-"" 
        if (f==2) dat<-sc3; fname<-"_noCpG"
        
        WilcoxTest.nt<-data.frame(matrix(ncol=6,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt)<-c("NT1","NT2","test","P.value","mean.SC.nt1","mean.SC.nt2")
        #SCs.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$EstSC[dat$ref==Ncomb[i,1]]
                vec2<-dat$EstSC[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
                
                WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt$test[i]<-"less"
                WilcoxTest.nt$P.value[i]<-result[[3]]
                WilcoxTest.nt$mean.SC.nt1[i]<-mean(vec1)
                WilcoxTest.nt$mean.SC.nt2[i]<-mean(vec2)
        } 
        
        WilcoxTest.nt2<-data.frame(matrix(ncol=6,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","P.value","mean.SC.nt1","mean.SC.nt2")
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$EstSC[dat$ref==Ncomb[i,1]]
                vec2<-dat$EstSC[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
                
                WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt2$test[i]<-"greater"
                WilcoxTest.nt2$P.value[i]<-result[[3]]
                WilcoxTest.nt2$mean.SC.nt1[i]<-mean(vec1)
                WilcoxTest.nt2$mean.SC.nt2[i]<-mean(vec2)
        }   
        
        WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
        write.csv(WilcoxTest.nt,paste0("Output1A/SelCoeff/SC_WilcoxTestResults_byNT", fname,".csv"))
}




