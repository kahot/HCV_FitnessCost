# Calculate selection coefficients for each site and run statistical tests
source("Rscripts/Pcorrection.R")
source("Rscripts/baseRscript.R")

#Read the estiamted mutatio rates of Geller's
mutrates<-read.csv("Output/Geller/Geller.MutRates.Summary_updated.csv", row.names = 1, stringsAsFactors = F)
#Read the transition mutation frequency file
df<-read.csv("Output/Mutfreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names = 1 )
df<-df[,196:203]
        
#Estiamte SC
df$TSmutrate[df$ref=="a"]<-mutrates$MutRate[mutrates$Mutation=="AG"]
df$TSmutrate[df$ref=="c"]<-mutrates$MutRate[mutrates$Mutation=="CU"]
df$TSmutrate[df$ref=="g"]<-mutrates$MutRate[mutrates$Mutation=="GA"]
df$TSmutrate[df$ref=="t"]<-mutrates$MutRate[mutrates$Mutation=="UC"]
df$EstSC<-""
for (j in 1:nrow(df)){
        df$EstSC[j] <- EstimatedS(df$TSmutrate[j],df$mean[j])
}
df$EstSC<-as.numeric(df$EstSC)

df<-df[df$pos>341&df$pos<8575,]
write.csv(df,"Output/SelCoeff/SC.csv")


## Test of SCs by type and by CpG
#df<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)

# Use A & T only for CpG Analysis
ty<-which(colnames(df)=="Type")
dat<-df

k=1
TypeList<-list()
scsum<-data.frame(Type=c("syn", "nonsyn","stop"))
for (typeofsite in c("syn", "nonsyn","stop")){
        all<-dat$EstSC[dat[,ty]==typeofsite]
        allcpg<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==1]
        allnoncpg1<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==0 & dat$ref=="t"]
        allnoncpg2<-dat$EstSC[dat[,ty]==typeofsite & dat$makesCpG==0 & dat$ref=="a"]

        typev<-paste0(typeofsite,"_all")
        assign(typev, all)
        scsum$mean[k]<-mean(all, na.rm=T)
        scsum$se[k]<-std.error(all,na.rm=T)
        scsum$CI[k]<-std.error(all)*1.96
        
        cpgv1<-paste0(typeofsite,"_allCpG")
        assign(cpgv1, allcpg)
        
        cpgv0<-paste0(typeofsite,"_allnonCpG")
        allnoncpg<-c(allnoncpg1,allnoncpg2)
        assign(cpgv0, allnoncpg)
        k=k+1
}

scsum[4,]<-"all"
scsum$mean[4]<-mean(dat$EstSC, na.rm=T)
scsum$se[4]<-std.error(dat$EstSC,na.rm=T)
scsum$CI[4]<-std.error(dat$EstSC,na.rm=T)*1.96
write.csv(scsum, "Output/SelCoeff/SC_summary_type.csv")

#wilcox.test 
wilcoxtest1<-data.frame("test"=matrix(nrow=4))
        
re1<-wilcox.test(syn_all,nonsyn_all, alternative = "less", paired = FALSE) 
wilcoxtest1$test[1]<-re1[[7]]
wilcoxtest1$rawP[1]<-re1[[3]]
        
re2<-wilcox.test(nonsyn_all,stop_all, alternative = "less", paired = FALSE) 
wilcoxtest1$test[2]<-re2[[7]]
wilcoxtest1$rawP[2]<-re2[[3]]
        
re3<-wilcox.test(syn_allCpG,syn_allnonCpG, alternative = "greater", paired = FALSE) 
wilcoxtest1$test[3]<-re3[[7]]
wilcoxtest1$rawP[3]<-re3[[3]]
        
re4<-wilcox.test(nonsyn_allCpG,nonsyn_allnonCpG, alternative = "greater", paired = FALSE) 
wilcoxtest1$test[4]<-re4[[7]]
wilcoxtest1$rawP[4]<-re4[[3]]

wilcoxtest1<-Pcorrection(wilcoxtest1)
write.csv(wilcoxtest1,"Output/SelCoeff/WilcoxonResults_by_Type.csv")

  
## SC by NT and CpG
S<-data.frame()
se<-data.frame()
S_CpG<-data.frame()
se_CpG<-data.frame()
S_nonCpG<-data.frame()
se_nonCpG<-data.frame()

typeofsite<-c("syn", "nonsyn","stop")
wtnt<-c("a", "t", "c", "g")
k=1
SC<-list()
for (i in 1:3){
     for (j in 1:4 ){
        sc<- dat[dat[,ty]==typeofsite[i] & dat$ref==wtnt[j], c("EstSC","se")]
        S[typeofsite[i],wtnt[j]]<-mean(sc$EstSC, na.rm=T)
        se[typeofsite[i],wtnt[j]]<-std.error(mean(sc$EstSC, na.rm=T)) 
        
        S_NoCpG<-dat[dat$Type==typeofsite[i] & dat$ref==wtnt[j] & dat$makesCpG==0,c("EstSC","se")]
        S_nonCpG[typeofsite[i],wtnt[j]]<-mean(S_NoCpG$EstSC, na.rm=T)
        se_nonCpG[typeofsite[i],wtnt[j]]<-std.error(mean(S_NoCpG$EstSC, na.rm=T))
        
        Sc_CpG<-dat[dat$Type==typeofsite[i] & dat$ref==wtnt[j] & dat$makesCpG==1,c("EstSC","se") ]
        S_CpG[typeofsite[i],wtnt[j]]<-mean(Sc_CpG$EstSC, na.rm=T)
        se_CpG[typeofsite[i],wtnt[j]]<-std.error(mean(Sc_CpG$EstSC, na.rm=T))
        
        SC[[k]]<-sc
        names(SC)[k]<-paste0(typeofsite[i],"_",wtnt[j])
        SC[[k+1]]<-S_NoCpG
        names(SC)[k+1]<-paste0(typeofsite[i],"_",wtnt[j],"_noncpg")
        SC[[k+2]]<-Sc_CpG
        names(SC)[k+2]<-paste0(typeofsite[i],"_",wtnt[j],"_cpg")
        k=k+3
     }
}
rownames(S_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(S_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")

rownames(se_nonCpG)<-c("syn_noncpg","nonsyn_noncpg","stop_noncpg")
rownames(se_CpG)<-c("syn_cpg","nonsyn_cpg","stop_cpg")

SCbyType<-rbind(S,S_nonCpG,S_CpG)
seByType<-rbind(se, se_nonCpG,se_CpG)

write.csv(SCbyType,"Output/SelCoeff/SC.mean_byNt_byTypeCpG.csv")
write.csv(seByType,"Output/SelCoeff/SC.se_byNt_byTypeCpG.csv")


# Wilcoxin Test of SC by nucleotide by CpG 
nuc.sc<-data.frame("syn.cpg"=matrix(nrow=4))
rownames(nuc.sc)<-c("a","t","c","g")

WilcoxTest.nt.sc<-data.frame(matrix(ncol=3,nrow=6))
colnames(WilcoxTest.nt.sc)<-c("nt","test","rawP")
k=1
for (i in c("a","t","c","g")) {
    if (i=="a"|i=="t"){
        syncpg<-SC[[paste0("syn_",i,"_cpg")]]
        syncpg<-syncpg$EstSC
        synnoncpg<-SC[[paste0("syn_",i,"_noncpg")]]
        synnoncpg<-synnoncpg$EstSC
        nonsyncpg<-SC[[paste0("nonsyn_",i,"_cpg")]]
        nonsyncpg<-nonsyncpg$EstSC
        nonsynnoncpg<-SC[[paste0("nonsyn_",i,"_noncpg")]]
        nonsynnoncpg<-nonsynnoncpg$EstSC
    
        if (i=="a"){
                result1<-wilcox.test(syncpg, synnoncpg, alternative = "greater", paired = FALSE) 
                result2<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "greater", paired = FALSE)
                for (r in 1:2){
                        result<-get(paste0('result',r))
                        WilcoxTest.nt.sc$nt[r]<-i
                        WilcoxTest.nt.sc$test[r]<-result[[7]]
                        WilcoxTest.nt.sc$rawP[r]<-result[[3]]}
        }
        if (i=="t"){
                result3<-wilcox.test(syncpg, synnoncpg, alternative = "greater", paired = FALSE) 
                result4<-wilcox.test(nonsyncpg,nonsynnoncpg,alternative = "greater", paired = FALSE)
                for (r in 3:4){
                        result<-get(paste0('result',r))
                        WilcoxTest.nt.sc$nt[r]<-i
                        WilcoxTest.nt.sc$test[r]<-result[[7]]
                        WilcoxTest.nt.sc$rawP[r]<-result[[3]]}
        }
    }
    else {  synnoncpg<-SC[[paste0("syn_",i,"_noncpg")]]
            synnoncpg<-synnoncpg$EstSC
            nonsynnoncpg<-SC[[paste0("nonsyn_",i,"_noncpg")]]
            nonsynnoncpg<-nonsynnoncpg$EstSC
            stopnoncpg<-SC[[paste0("stop_",i,"_noncpg")]]
            stopnoncpg<-stopnoncpg$EstS
        
        if (i =="c") {
                result5<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "less", paired = FALSE) 
                WilcoxTest.nt.sc$nt[5]<-i
                WilcoxTest.nt.sc$test[5]<-result5[[7]]
                WilcoxTest.nt.sc$rawP[5]<-result5[[3]]     }           
        if (i =="g") {
                result6<-wilcox.test(synnoncpg,nonsynnoncpg, alternative = "less", paired = FALSE) 
                WilcoxTest.nt.sc$nt[6]<-i
                WilcoxTest.nt.sc$test[6]<-result6[[7]]
                WilcoxTest.nt.sc$rawP[6]<-result6[[3]]    }
        }
     k=k+1
} 

WilcoxTest.nt.sc<-Pcorrection(WilcoxTest.nt.sc)

write.csv(WilcoxTest.nt.sc,"Output/SelCoeff/WilcoxTestResults_byNTbyCpgType.csv")


## Wilcoxon test on SCs by gene 
df<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
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

#sc by gene
gene.summary<-data.frame(gene=genenames[2:12])
for (i in 1:nrow(gene.summary)){
    sc2<-sc[sc$gene==gene.summary$gene[i],]
    gene.summary$mean[i]<-mean(sc2$EstSC, na.rm = T)
    gene.summary$SE[i]<-std.error(sc2$EstSC, na.rm=T)
    gene.summary$CI[i]<-std.error(sc2$EstSC, na.rm=T)*1.96
    
}
write.csv(gene.summary, "Output/SelCoeff/SC_byGene.csv")

sc2$EstSC


Gcomb<-t(combn(genenames[2:12],2))

WilcoxTest2.gene<-data.frame(matrix(ncol=4,nrow=nrow(Gcomb)))
colnames(WilcoxTest2.gene)<-c("gene1","gene2","test","rawP")

for (i in 1:nrow(Gcomb)) {
        vec1<-sc$EstSC[sc$gene==Gcomb[i,1]]
        vec2<-sc$EstSC[sc$gene==Gcomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest2.gene$gene1[i]<-Gcomb[i,1]
        WilcoxTest2.gene$gene2[i]<-Gcomb[i,2]
        WilcoxTest2.gene$test[i]<-"less"
        WilcoxTest2.gene$rawP[i]<-result[[3]]
}   

WilcoxTest2.gene2<-data.frame(matrix(ncol=4,nrow=nrow(Gcomb)))
colnames(WilcoxTest2.gene2)<-c("gene1","gene2","test","rawP")

for (i in 1:nrow(Gcomb)) {
        vec1<-sc$EstSC[sc$gene==Gcomb[i,1]]
        vec2<-sc$EstSC[sc$gene==Gcomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
        
        WilcoxTest2.gene2$gene1[i]<-Gcomb[i,1]
        WilcoxTest2.gene2$gene2[i]<-Gcomb[i,2]
        WilcoxTest2.gene2$test[i]<-"greater"
        WilcoxTest2.gene2$rawP[i]<-result[[3]]
}        

Wilcox.gene<-rbind(WilcoxTest2.gene,WilcoxTest2.gene2)
Wilcox.gene<-Pcorrection(Wilcox.gene)

write.csv(Wilcox.gene,"Output/SelCoeff/SC_WilcoxTestResults_byGene.csv")


#####
## wilcox test of SC by NT
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))
WilcoxTest.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt)<-c("NT1","NT2","test","rawP")

for (i in 1:nrow(Ncomb)) {
        vec1<-df$EstSC[df$ref==Ncomb[i,1]]
        vec2<-df$EstSC[df$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 

        WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt$test[i]<-"greater"
        WilcoxTest.nt$rawP[i]<-result[[3]]
}   

WilcoxTest.nt2<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","rawP")

for (i in 1:nrow(Ncomb)) {
        vec1<-df$EstSC[df$ref==Ncomb[i,1]]
        vec2<-df$EstSC[df$ref==Ncomb[i,2]]
        result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
        
        WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
        WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
        WilcoxTest.nt2$test[i]<-"less"
        WilcoxTest.nt2$rawP[i]<-result[[3]]
}   

WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
WilcoxTest.nt<-Pcorrection(WilcoxTest.nt)
write.csv(WilcoxTest.nt,"Output/SelCoeff/SC_WilcoxTestResults_byNT.csv")


########

sc3<-df[df$makesCpG==0,]

## wilcox test by NT -with or without CpG sites
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))
for (f in 1:2){
        if (f==1) dat<-df; fname<-"" 
        if (f==2) dat<-sc3; fname<-"_noCpG"
        
        WilcoxTest.nt<-data.frame(matrix(ncol=6,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt)<-c("NT1","NT2","test","rawP","mean.SC.nt1","mean.SC.nt2")
        #SCs.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$EstSC[dat$ref==Ncomb[i,1]]
                vec2<-dat$EstSC[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
                
                WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt$test[i]<-"less"
                WilcoxTest.nt$rawP[i]<-result[[3]]
                WilcoxTest.nt$mean.SC.nt1[i]<-mean(vec1)
                WilcoxTest.nt$mean.SC.nt2[i]<-mean(vec2)
        } 
        
        WilcoxTest.nt2<-data.frame(matrix(ncol=6,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","rawP","mean.SC.nt1","mean.SC.nt2")
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$EstSC[dat$ref==Ncomb[i,1]]
                vec2<-dat$EstSC[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
                
                WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt2$test[i]<-"greater"
                WilcoxTest.nt2$rawP[i]<-result[[3]]
                WilcoxTest.nt2$mean.SC.nt1[i]<-mean(vec1)
                WilcoxTest.nt2$mean.SC.nt2[i]<-mean(vec2)
        }   
        
        
        WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
        WilcoxTest.nt<-Pcorrection(WilcoxTest.nt)
        write.csv(WilcoxTest.nt,paste0("Output/SelCoeff/SC_WilcoxTestResults_byNT", fname,".csv"))
}




