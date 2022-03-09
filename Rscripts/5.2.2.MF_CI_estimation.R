#estimate confidence intervals for different type of mutation frequencies
library(binom)
source("Rscripts/baseRscript.R")

Reads<-read.csv("Output/ReadDepth_All.csv", row.names=1)
reads<-Reads[Reads$pos>341&Reads$pos<8575,]

#Calculate CI for all positions using "combining stratum specific F-distribution confidence intervals" from Waller et al. 1994
#Treating each population/patient as stratum -so that avoiding overconservative CI estiamtes

#first calculate CI for each position within a population based on CI.lower =(wx/(x+(n-x+1)*F(alpha/2, 2(n-x+1),2x)))
# CI.upper=w*(x+1)*F(alpha/2, 2(x+1), 2(n-x)/((n-x)+(x+1)*F(alpha/2, 2(x+1), 2(n-x))
# x= number of mutation, n= read depth, w=weight (1/sample size)
#Sum all CIs for each position and apply rescaling factor

#Read SeqData files
seqfiles<-list.files("Output/SeqData/",pattern="SeqData_")
seq<-list()
for (i in 1:length(seqfiles)){ 
    df<-read.csv(paste0("Output/SeqData/",seqfiles[i]),stringsAsFactors=FALSE, row.names=1)
    seq[[i]]<-df
    names(seq)[i]<-substr(paste(seqfiles[i]),start=9,stop=15)
}


#get MF Summary files 
files<-c("Ts","Tv1","Tv2","Tvs","All")
mf.files<-list()
summary.mf<-data.frame(type=files)
for (i in 1:length(files)){
    df<-read.csv(paste0("Output/MutFreq/Filtered.",files[i],".Q35.csv"), row.names = 1)
    mf.files[[i]]<-df
    names(mf.files)[i]<-files[i]
    summary.mf$All[i]<-mean(df$mean, na.rm=T)
}


#Add gene information for Ts mf file
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
ts<-mf.files[[1]]
mf.files[[1]]<-merge(ts, genetable, by="pos", all.x=T )

genes<-genes[2:12,]

#remove filtered sites from reads
ts<-ts[ts$pos>341&ts$pos<8575,]
reads3<-reads[,2:196]
reads3[is.na(ts[,1:195])]<-NA
reads3$pos<-reads$pos

Posi.used<-ts$pos
#set alpha to 0.05
p=0.05

#calculate CIs for each position minus the weight
for (f in 1:5){
    mf<-mf.files[[f]]
    mf<-mf[mf$pos>341&mf$pos<8575,]
    print(names(mf.files)[f])
    
    #non-weighted persite CI
    df_low<-data.frame(pos=mf$pos)
    df_up<-data.frame(pos=mf$pos)
    #agresti-coull method
    df_low2<-data.frame(pos=mf$pos)
    df_up2<-data.frame(pos=mf$pos)
    
    for (i in 1:195){
        na.positions<-which(is.na(mf[,i]))
        seq.df<-seq[[i]]
        seq.df<-seq.df[seq.df$pos %in% Posi.used,]
        
        #remove the filtered out sites from seq.df
        seq.df[na.positions, c(1:4, 8)]<-NA
        
        # calculate x (# of mutation of interest)
        if (f==1) seq.df$x<-as.integer(apply(seq.df, 1, function(x) x[paste0(x["transition.ref"])] ))
        if (f==2) seq.df$x<-as.integer(apply(seq.df, 1, function(x) x[transv1(x["ref"])]))
        if (f==3) seq.df$x<-as.integer(apply(seq.df, 1, function(x) x[transv2(x["ref"])]))
        if (f==4) seq.df$x<-as.integer(apply(seq.df, 1, function(x) as.integer(x[transv1(x["ref"])])+as.integer(x[transv2(x["ref"])])))
        if (f==5) seq.df$x<-as.integer(apply(seq.df, 1, function(x) as.integer(x[paste0(x["transition.ref"])])+as.integer(x[transv1(x["ref"])])+as.integer(x[transv2(x["ref"])])))

        
        # calculate n1=n-x+1, here n=TotalReads
        seq.df$n1<-seq.df$TotalReads-seq.df$x+1
        #calculate F values
        seq.df$f1<-apply(seq.df[,c("x","n1")], 1, function(x) ifelse(x["x"]==0, 0, qf(p/2, x["n1"]*2, x["x"]*2, lower.tail=FALSE)))
        seq.df$f2<-qf(p/2, (seq.df$x+1)*2, (seq.df$TotalReads-seq.df$x)*2, lower.tail=FALSE)
        df_low[,names(seq)[i]]<-apply(seq.df[,c("x","n1","f1")], 1, function(x) ifelse(x["x"]==0, 0,x["x"]/(x["x"]+x["n1"]*x["f1"]))) 
        df_up[,names(seq)[i]]<-((seq.df$x+1)*seq.df$f2)/((seq.df$TotalReads-seq.df$x)+(seq.df$x+1)*seq.df$f2)
        
        #Agresti-Coull CIs
        #aggregate across the population
        for (n in 1:nrow(seq.df)){
            if (!is.na(seq.df$x[n])){
                ci.results<- binom.confint(seq.df$x[n], seq.df$TotalReads[n], conf.level=0.95, method="agresti-coull")
                seq.df$ci_low[n]<-ci.results$lower
                seq.df$ci_up[n]<-ci.results$upper
            }
            if (is.na(seq.df$x[n])) {
                seq.df$ci_low[n]<-NA
                seq.df$ci_up[n]<-NA
            }
        }
        seq.df$ci_low[seq.df$ci_low<0]<-0
        
        df_low2[,names(seq)[i]]<-seq.df$ci_low
        df_up2[,names(seq)[i]]<-seq.df$ci_up
    }
    
    
    #calculate the total available to obtain w
    df_low$w<-apply(df_low[2:196], 1, function(x) length(x[!is.na(x)]))
    
    #write.csv(df_low, paste0("Output/MutFreq/ci_estimates_lower_",files[f],".csv"))
    #write.csv(df_up, paste0("Output/MutFreq/ci_estimates_upper_",files[f],".csv"))
    
    #CI esitmation for each position across populations
    df.CI<-data.frame(pos=mf$pos)
    
    #rescaling factor for each site. Assing an equal weight for all samples at each position
    reads3$wi2<-(1/df_low$w)^2
    reads3$w<-1/df_low$w
    r1=apply(reads3, 1, function(x) sum(1/x[1:195], na.rm=T)*x["wi2"]^0.5)
    r2=apply(reads3, 1, function(x) sum(1/sqrt(x[1:195]), na.rm=T)*x["w"])
    
    R=r1/r2
    df.CI[,"lower"]<-rowSums(df_low[2:196], na.rm=T)/df_low$w
    df.CI[,"upper"]<-rowSums(df_up[2:196], na.rm=T)/df_low$w
    #CI adjusted by the rescaling factor
    df.CI[,"lower.R"]<-(mf$mean-df.CI$lower)*R
    df.CI[,"upper.R"]<-(df.CI$upper-mf$mean)*R
    
    write.csv(df.CI, paste0("Output/MutFreq/CIs.",files[f],".csv"))
    
    
    #calculate Agresti-coull CIs for each site
    df.CIac<-data.frame(pos=mf$pos)
    df.CIac$lower=rowMeans(df_low2[2:196], na.rm=T)
    df.CIac$upper=rowMeans(df_up2[2:196], na.rm=T)
    
    write.csv(df.CIac, paste0("Output/MutFreq/CIs_Agresti–Coull.",files[f],".csv"))
}
    
    
  
#create a summary table for different types of mutations
#1. transition mutations
type<-c("syn","nonsyn","stop")

CI.sum<-data.frame(Type=c("All", type)) # for F estimates
CI.sum2<-data.frame(Type=c("All", type)) #for AC estimates
CI.sum3<-data.frame(Type=c("All", type)) #for F without scaling

#summary for each gene for  transition mutations
gene.summary.syn<-data.frame(gene=genes$Gene)
gene.summary.nonsyn<-data.frame(gene=genes$Gene)
gene.summary.syn2<-data.frame(gene=genes$Gene)  #for AC estimates
gene.summary.nonsyn2<-data.frame(gene=genes$Gene)  #for AC estimates
gene.summary.syn3<-data.frame(gene=genes$Gene)  #for F without scaling
gene.summary.nonsyn3<-data.frame(gene=genes$Gene)#for F without scaling

#1. for transition mutations
ci.ts<-read.csv("Output/MutFreq/CIs.Ts.csv", row.names = 1)
ci.ts2<-read.csv("Output/MutFreq/CI_Agresti–Coull_Ts_mean.csv", row.names = 1)
mf<-mf.files[[1]]
mf<-mf[mf$pos>341&mf$pos<8575,]

CI.sum[1,"low"]<-mean(ci.ts$lower.R, na.rm=T)
CI.sum[1,"up"]<-mean(ci.ts$upper.R, na.rm=T)
CI.sum2[1,"low"]<-mean(ci.ts2$lower, na.rm=T)
CI.sum2[1,"up"]<-mean(ci.ts2$upper, na.rm=T)
CI.sum3[1,"low"]<-mean(ci.ts$lower, na.rm=T)
CI.sum3[1,"up"]<-mean(ci.ts$upper, na.rm=T)

for (i in 1:3){
    positions<-mf$pos[mf$Type==type[i]]
    df<-ci.ts[ci.ts$pos %in% positions,]
    df2<-ci.ts2[ci.ts2$pos %in% positions,]
    CI.sum[(i+1),"low"]<-mean(df$lower.R, na.rm=T)
    CI.sum[(i+1),"up"]<- mean(df$upper.R, na.rm=T)
    CI.sum2[(i+1),"low"]<-mean(df2$lower, na.rm=T)
    CI.sum2[(i+1),"up"]<- mean(df2$upper, na.rm=T)
    CI.sum3[(i+1),"low"]<-mean(df$lower, na.rm=T)
    CI.sum3[(i+1),"up"]<- mean(df$upper, na.rm=T)
    
    #mut freq means
    summary.mf[1,type[i]]<-mean(mf$mean[mf$Type==type[i]], na.rm=T)
    
    #For each gene
    for (j in 1:nrow(genes)){
        pos2<-mf$pos[mf$gene==genes$Gene[j]&mf$Type==type[i]]
        df.g1<-ci.ts[ci.ts$pos %in% pos2,]
        df.g2<-ci.ts2[ci.ts2$pos %in% pos2,]
        if (type[i]=="syn"){
            gene.summary.syn$low[j]<- mean(df.g1$lower.R, na.rm=T)
            gene.summary.syn$up[j] <- mean(df.g1$upper.R,na.rm=T)
            gene.summary.syn2$low[j]<-mean(df.g2$lower, na.rm=T)
            gene.summary.syn2$up[j] <-mean(df.g2$upper,na.rm=T)
            gene.summary.syn3$low[j]<-mean(df.g1$lower, na.rm=T)
            gene.summary.syn3$up[j] <-mean(df.g1$upper,na.rm=T)
        }
        if (type[i]=="nonsyn"){
            gene.summary.nonsyn$low[j]<- mean(df.g1$lower.R, na.rm=T)
            gene.summary.nonsyn$up[j] <- mean(df.g1$upper.R,na.rm=T)
            gene.summary.nonsyn2$low[j]<-mean(df.g2$lower,na.rm=T)
            gene.summary.nonsyn2$up[j] <-mean(df.g2$upper,na.rm=T)
            gene.summary.nonsyn3$low[j]<-mean(df.g1$lower, na.rm=T)
            gene.summary.nonsyn3$up[j] <-mean(df.g1$upper,na.rm=T)
        }
    }
}

#add Ts to CI.sum
CI.sum$Mutation<-"Ts"
CI.sum2$Mutation<-"Ts"
CI.sum3$Mutation<-"Ts"
  
#combine the syn and nonsyn tables

gene.summary.syn$Type<-"Syn"
gene.summary.syn2$Type<-"Syn"
gene.summary.syn3$Type<-"Syn"
gene.summary.nonsyn$Type<-"Nonsyn"
gene.summary.nonsyn2$Type<-"Nonsyn"
gene.summary.nonsyn3$Type<-"Nonsyn"

gene.summary<-rbind(gene.summary.syn, gene.summary.nonsyn)
gene.summary2<-rbind(gene.summary.syn2, gene.summary.nonsyn2)
gene.summary3<-rbind(gene.summary.syn3, gene.summary.nonsyn3)

write.csv(gene.summary, "Output/MutFreq/CI.Ts.byGene.csv")
write.csv(gene.summary2, "Output/MutFreq/CI.Ts.byGene_Agresti-coull.csv")
write.csv(gene.summary3, "Output/MutFreq/CI.Ts.byGene_noRescaling.csv")


# trasnversion CI for different types of mutations

mftv<-mf.files[[4]] #tvs
mf1<-mf.files[[2]]
mf2<-mf.files[[3]]

mftv<-mftv[mftv$pos>341&mftv$pos<8575,]
mf1<-mf1[mf1$pos>341&mf1$pos<8575,]
mf2<-mf2[mf2$pos>341&mf2$pos<8575,]

CI.sum.tv<-data.frame(Type=c("All", type)) # for F estimates
CI.sum.tv2<-data.frame(Type=c("All", type)) #for AC estimates
CI.sum.tv3<-data.frame(Type=c("All", type)) #for F without rescaling

tvs<-read.csv("Output/MutFreq/CIs.Tvs.csv", row.names = 1)
tv1<-read.csv("Output/MutFreq/CIs.Tv1.csv", row.names = 1)
tv2<-read.csv("Output/MutFreq/CIs.Tv2.csv", row.names = 1)

tvs.ac<-read.csv("Output/MutFreq/CIs_Agresti–Coull.Tvs.csv", row.names = 1)
tv1.ac<-read.csv("Output/MutFreq/CIs_Agresti–Coull.Tv1.csv", row.names = 1)
tv2.ac<-read.csv("Output/MutFreq/CIs_Agresti–Coull.Tv2.csv", row.names = 1)


CI.sum.tv[1,"low"]<-mean(tvs$lower.R, na.rm=T)
CI.sum.tv[1,"up"]<-mean(tvs$upper.R, na.rm=T)

CI.sum.tv2[1,"low"]<-mean(tvs.ac$lower, na.rm=T) #for agresti-coull method
CI.sum.tv2[1,"up"]<-mean(tvs.ac$upper, na.rm=T)

CI.sum.tv3[1,"low"]<-mean(tvs$lower, na.rm=T) 
CI.sum.tv3[1,"up"] <-mean(tvs$upper, na.rm=T)


type<-c("syn","nonsyn","stop")
#reads3$pos<-reads$pos
for (i in 1: length(type)){
    pos.b<-mftv$pos[mftv$Type.tv1==type[i]&mftv$Type.tv2==type[i]]
    pos.1<-mftv$pos[mftv$Type.tv1==type[i]&mftv$Type.tv2!=type[i]]
    pos.2<-mftv$pos[mftv$Type.tv1!=type[i]&mftv$Type.tv2==type[i]] 
    
    df1<-tv1[tv1$pos %in% pos.1,]
    df2<-tv2[tv2$pos %in% pos.2,]
    df3<-tvs[tvs$pos %in% pos.b,]
    
    df<-rbind(df1, df2,df3)
    
    dfac1<-tv1.ac[tv1.ac$pos %in% pos.1,]
    dfac2<-tv2.ac[tv2.ac$pos %in% pos.2,]
    dfac3<-tvs.ac[tvs.ac$pos %in% pos.b,]
    
    dfac<-rbind(dfac1, dfac2,dfac3)
    
    #average mut freq
    summary.mf[4,type[i]]<-mean(c(mftv$mean[mftv$pos %in% pos.b], mf1$mean[mf1$pos %in% pos.1],  mf2$mean[mf2$pos %in% pos.2]),na.rm=T)
    
    
    #get average CIs
    CI.sum.tv[(i+1),"low"]<-mean(df$lower.R, na.rm=T)
    CI.sum.tv[(i+1),"up"]<- mean(df$upper.R, na.rm=T)
    
    CI.sum.tv2[(i+1),"low"]<-mean(dfac$lower, na.rm=T)
    CI.sum.tv2[(i+1),"up"]<- mean(dfac$upper, naac.rm=T)
    
    CI.sum.tv3[(i+1),"low"]<-mean(df$lower, na.rm=T)
    CI.sum.tv3[(i+1),"up"]<- mean(df$upper, na.rm=T)
}

write.csv(summary.mf, "Output/MutFreq/MF_summary.csv")


#Combine the tv and tvs and two methods
CI.sum.tv$Mutation<-"Tvs"
CI.sum.tv2$Mutation<-"Tvs"
CI.sum.tv3$Mutation<-"Tvs"

CI.sum.tv$Method<-"F"
CI.sum.tv2$Method<-"AC"
CI.sum.tv3$Method<-"no.rescaling"

CI.sum$Method<-"F"
CI.sum2$Method<-"AC"
CI.sum3$Method<-"no.rescaling"


CI.summary<-rbind(CI.sum, CI.sum2,CI.sum3,CI.sum.tv, CI.sum.tv2,CI.sum.tv3)

write.csv(CI.summary,"Output/MutFreq/CI_TypeofMutaions_summary.csv")

