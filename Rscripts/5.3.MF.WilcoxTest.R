# Run the statistical tests on mutation frequency data
source("Rscripts/Pcorrection.R")

### Test of mutation frequency among different mutation types

#test types
types<-c("syn","nonsyn","stop")
comb<-combn(types,2)
comb<-data.frame(t(comb))
# *no need to test syn vs. stop
comb<-comb[-2,]
co_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

#Create a combined transversion syn/nonsyn data frame
## Load the data
Tv1<-read.csv("Output/MutFreq/Filtered.Tv1.Q35.csv", row.names = 1, stringsAsFactors = F)
Tv1<-Tv1[Tv1$pos>=342,1:199 ]
colnames(Tv1)[199]<-"Type"
Tv2<-read.csv("Output/MutFreq/Filtered.Tv2.Q35.csv", row.names = 1, stringsAsFactors = F)
Tv2<-Tv2[Tv2$pos>=342, 1:199]
colnames(Tv2)[199]<-"Type"
Tv<- rbind(Tv1, Tv2)
Ts<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv", row.names = 1, stringsAsFactors = F)
Ts<-Ts[Ts$pos>=342, ]

#Test syn vs. nonsyn mutations for transitions and transversion mutations
testResults<-data.frame()
for (f in 1:2){
    if (f==1) {DF<-Ts; muttype<-"Ts"}
    if (f==2) {DF<-Tv; muttype<-"Tvs"}
    
    testRe<-data.frame(test=co_pairs)
    for (i in 1:length(co_pairs)){
        re<-wilcox.test(DF$mean[DF$Type==comb[i,1]], DF$mean[DF$Type==comb[i,2]], alternative = "greater", paired = FALSE) 
        testRe$rawP[i]<-re[[3]]
    }
    testRe$Type<-muttype
    
    testResults<-rbind(testResults,testRe)
    
}

# Add results for all transition vs. transversion 
tvs<-c(Tv1$mean, Tv2$mean)
r1<-wilcox.test(Ts$mean, tvs, alternative = "greater", paired = FALSE) 

testResults[5,]<-c("Ts.vs.Tvs", r1[[3]], "all" )
testResults$rawP<-as.numeric(testResults$rawP)

#Run the correction
testResults<-Pcorrection(testResults)
#Save teh results
write.csv(testResults,"Output/MutFreq/WilcoxonResults_Ts.vs.Tvs.csv")

##############
# Transition mutations: test by nucleotide and by gene
# Remove CpG creating mutations to see how it differs
Ts2<-Ts[Ts$makesCpG==0,]

# 1. Wilcox test by NT
NT<-c("a","c","t","g")
Ncomb<-t(combn(NT,2))

for (f in 1:2){
        if (f==1) {dat<-Ts;  fname <- "" }
        if (f==2) {dat<-Ts2; fname<-"_noCpG"}
        
        WilcoxTest.nt<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt)<-c("NT1","NT2","test","rawP")
        
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$mean[dat$ref==Ncomb[i,1]]
                vec2<-dat$mean[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
                
                WilcoxTest.nt$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt$test[i]<-"less"
                WilcoxTest.nt$rawP[i]<-result[[3]]
        }   
        
        WilcoxTest.nt2<-data.frame(matrix(ncol=4,nrow=nrow(Ncomb)))
        colnames(WilcoxTest.nt2)<-c("NT1","NT2","test","rawP")
        
        for (i in 1:nrow(Ncomb)) {
                vec1<-dat$mean[dat$ref==Ncomb[i,1]]
                vec2<-dat$mean[dat$ref==Ncomb[i,2]]
                result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
                
                WilcoxTest.nt2$NT1[i]<-Ncomb[i,1]
                WilcoxTest.nt2$NT2[i]<-Ncomb[i,2]
                WilcoxTest.nt2$test[i]<-"greater"
                WilcoxTest.nt2$rawP[i]<-result[[3]]
        }   
        
        WilcoxTest.nt<-rbind(WilcoxTest.nt,WilcoxTest.nt2)
        WilcoxTest.nt<-Pcorrection(WilcoxTest.nt)
        write.csv(WilcoxTest.nt, paste0("Output/MutFreq/WilcoxonResults_byNT", fname,".csv"))
}

# *Removing CpG creating mutations will not change the overall results

## Test by gene
genes<-read.csv("Data/HCV_annotations2.csv", stringsAsFactors = F)
genes$Gene<-as.character(genes$Gene)
#genes$Gene[genes$Gene=="NS1(P7)"]<-"NS1"
genenames<-genes$Gene[2:12]
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
    gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector

#Add gene info to Ts
Ts<-merge(Ts, genetable, by="pos", all.x=T)

## Test MF by gene
Gcomb<-t(combn(genenames,2))
WilcoxTest.gene<-data.frame(matrix(ncol=4,nrow=nrow(Gcomb)))
colnames(WilcoxTest.gene)<-c("gene1","gene2","test","rawP")
for (i in 1:nrow(Gcomb)) {
    vec1<-Ts$mean[Ts$gene==Gcomb[i,1]]
    vec2<-Ts$mean[Ts$gene==Gcomb[i,2]]
    result<-wilcox.test(vec1, vec2, alternative = "less", paired = FALSE) 
    
    WilcoxTest.gene$gene1[i]<-Gcomb[i,1]
    WilcoxTest.gene$gene2[i]<-Gcomb[i,2]
    WilcoxTest.gene$test[i]<-"less"
    WilcoxTest.gene$rawP[i]<-result[[3]]
}   

#Test for another direction 
WilcoxTest.gene2<-data.frame(matrix(ncol=4,nrow=nrow(Gcomb)))
colnames(WilcoxTest.gene2)<-c("gene1","gene2","test","rawP")

for (i in 1:nrow(Gcomb)) {
    vec1<-Ts$mean[Ts$gene==Gcomb[i,1]]
    vec2<-Ts$mean[Ts$gene==Gcomb[i,2]]
    result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
    
    WilcoxTest.gene2$gene1[i]<-Gcomb[i,1]
    WilcoxTest.gene2$gene2[i]<-Gcomb[i,2]
    WilcoxTest.gene2$test[i]<-"greater"
    WilcoxTest.gene2$rawP[i]<-result[[3]]
}        

WilcoxTest.gene<-rbind(WilcoxTest.gene,WilcoxTest.gene2)
WilcoxTest.gene<-Pcorrection(WilcoxTest.gene)
write.csv(WilcoxTest.gene,"Output/MutFreq/WilcoxonResults_byGene.csv")

#hvr<-WilcoxTest.gene[WilcoxTest.gene$gene1=="HVR1"|WilcoxTest.gene$gene2=="HVR1",]
#core<-WilcoxTest.gene[WilcoxTest.gene$gene1=="Core"|WilcoxTest.gene$gene2=="Core",]

# Create a mean mf by gene table
gene.sum<-data.frame(aggregate(Ts[,"mean"], by=list(Ts$gene), mean, na.rm=T))



# Test syn vs. nonsyn mutations in each gene  
Gcomb<-t(combn(genenames,2))
WilcoxTest.gene<-data.frame(matrix(ncol=3,nrow=11))
colnames(WilcoxTest.gene)<-c("gene","test","rawP")

for (i in 1:11) {
    vec1<-Ts$mean[Ts$gene==genenames[i] & Ts$Type=="syn"]
    vec2<-Ts$mean[Ts$gene==genenames[i] & Ts$Type=="nonsyn"]
    result<-wilcox.test(vec1, vec2, alternative = "greater", paired = FALSE) 
    
    WilcoxTest.gene$gene[i]<-genenames[(i)]
    WilcoxTest.gene$test[i]<-"greater"
    WilcoxTest.gene$rawP[i]<-result[[3]]
}   

WilcoxTest.gene<-Pcorrection(WilcoxTest.gene)
write.csv(WilcoxTest.gene,"Output/MutFreq/WilcoxonResults_byGene.synvsNonsyn.csv")
