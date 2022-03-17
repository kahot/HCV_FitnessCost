# Create mutation frequency summary data frames and estimate CIs

#Read the filtered overview files
Files<-list.files("Output/OverviewF/",pattern="overviewF.csv")
Overviews<-list()
for (i in 1:length(Files)){ 
        overviews2<-read.csv(paste0("Output/OverviewF/",Files[i]),stringsAsFactors=FALSE, row.names=1)
        Overviews[[i]]<-overviews2
        names(Overviews)[i]<-substr(paste(Files[i]),start=1,stop=7)
}

#Create summary data frames
Ts<-data.frame(sapply(Overviews,"[[","freq.Ts.ref"))
Tv1<-data.frame(sapply(Overviews,"[[","freq.transv1.ref") )
Tv2<-data.frame(sapply(Overviews,"[[","freq.transv2.ref") )
Tvs<-data.frame(sapply(Overviews,"[[","freq.transv.ref") )
All<-data.frame(sapply(Overviews,"[[","freq.mutations.ref"))
MVF<-data.frame(sapply(Overviews,"[[","freq.mutations"))


#Attach the metadata and save them into a list
files<-c("Ts","Tv1","Tv2","Tvs","All","MVF")
cnames<-c("",".tv1",".tv2",".tvs",".all",",mvf")
mf.files<-list()
s<-length(Overviews)

#Extract metadata from one file and attached them to summary data frames
meta<-Overviews[[3]]
colnames(meta)[33:35]<-c("MutAA","MutAA.tv1","MutAA.tv2")
meta<-meta[,c(1,5:42)]
for (i in 1:length(files)) {
        dat<-get(files[i])
        dat$mean<-rowMeans(dat[1:s],na.rm=T)
        if (i==1|i==5|i==6) muttypes<-meta[,c("pos","ref", "Type","WTAA","MutAA","makesCpG","bigAAChange")]
        if(i==2|i==3)       muttypes<-meta[,c("pos","ref", paste0("Type",cnames[i]),"WTAA",paste0("MutAA",cnames[i]), paste0("makesCpG",cnames[i]), paste0("bigAAChange",cnames[i]))]
        if(i==4) muttypes<-meta[,c("pos","ref","Type.tv1","Type.tv2","WTAA","MutAA.tv1","MutAA.tv2","makesCpG.tv1","makesCpG.tv1","bigAAChange.tv1","bigAAChange.tv2")]
        
        dat2<-cbind(dat,muttypes)
        mf.files[[i]]<-dat2
        names(mf.files)[i]<-files[i]
        write.csv(dat2,paste0("Output/MutFreq/Filtered.",files[i],".Q35.csv"))
}


##############################

## Create the read depth summary for all positions for all samples
Reads<-data.frame(sapply(Overviews, "[[", "TotalReads"))
df<-Overviews[[1]]
Reads<-cbind(df[,"pos"],Reads)
colnames(Reads)[1]<-"pos"
write.csv(Reads, "Output/ReadDepth_All.csv")



##############################

# Calculate average and median read depth for all sites
Reads2<-data.frame(Seq=names(Overviews))
for (i in 1:length(Overviews)){
    dat<-Overviews[[i]]
    dat<-dat[!is.na(dat$freq.Ts),]
    Reads2$ave[i]<-mean(dat$TotalReads, na.rm=T) 
    Reads2$median[i]<-median(dat$TotalReads, na.rm = T)
    Reads2$max[i]<-max(dat$TotalReads, na.rm=T)
    Reads2$min[i]<-min(dat$TotalReads, na.rm=T)
    
}

mean(Reads2$ave, na.rm=T) #6218.767
mean(Reads2$median, na.rm=T) #5282.449


### Calculate SE and CI
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

#Reads<-read.csv("Output/ReadDepth_All.csv", row.names=1)
#reads<-Reads[Reads$pos>341&Reads$pos<8575,]


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
ts<-merge(ts, genetable, by="pos", all.x=T )
ts<-ts[,c(2:196, 1,197:204)]
mf.files[[1]]<-ts
genes<-genes[2:12,]

#create a summary table for different types of mutations
type<-c("syn","nonsyn","stop")
CI.sum<-data.frame(Type=c("all", type)) # for F estimates

#1. transition mutations

mf<-mf.files[[1]]
mf<-mf[mf$pos>341&mf$pos<8575,]
    
CI.sum$SE[1]<-std.error(mf$mean, na.rm=T)
CI.sum$CI[1]<-std.error(mf$mean, na.rm=T)*1.96

#summary for each gene for  transition mutations
gene.summary.syn<-data.frame(gene=genes$Gene)
gene.summary.nonsyn<-data.frame(gene=genes$Gene)

#1. for transition mutations

for (i in 1:3){
    df<-mf[mf$Type==type[i],]
    CI.sum[(i+1),"SE"]<-std.error(df$mean, na.rm=T)
    CI.sum[(i+1),"CI"]<- 1.96*std.error(df$mean, na.rm=T)
    
    #mut freq means
    summary.mf[1,type[i]]<-mean(mf$mean[mf$Type==type[i]], na.rm=T)
    
    #For each gene
    for (j in 1:nrow(genes)){
        df<-mf[mf$gene==genes$Gene[j]&mf$Type==type[i],]
        if (type[i]=="syn"){
            gene.summary.syn$mean[j]<-mean(df$mean, na.rm=T)
            gene.summary.syn$se[j]<- std.error(df$mean, na.rm=T)
            gene.summary.syn$ci[j] <- 1.96*std.error(df$mean, na.rm=T)
        }
        if (type[i]=="nonsyn"){
            gene.summary.nonsyn$mean[j]<-mean(df$mean, na.rm=T)
            gene.summary.nonsyn$se[j]<- std.error(df$mean, na.rm=T)
            gene.summary.nonsyn$ci[j] <- 1.96*std.error(df$mean, na.rm=T)
        }
    }
}

#add Ts to CI.sum
CI.sum$Mutation<-"Ts"

#combine the syn and nonsyn tables
gene.summary.syn$Type<-"syn"
gene.summary.nonsyn$Type<-"nonsyn"

gene.summary<-rbind(gene.summary.syn, gene.summary.nonsyn)
write.csv(gene.summary, "Output/MutFreq/SE.Ts.byGene.csv")


# trasnversion SE/CI for different types of mutations

mftv<-mf.files[[4]] #tvs
mf1<-mf.files[[2]]
mf2<-mf.files[[3]]

mftv<-mftv[mftv$pos>341&mftv$pos<8575,]
mf1<-mf1[mf1$pos>341&mf1$pos<8575,]
mf2<-mf2[mf2$pos>341&mf2$pos<8575,]

#summary table
CI.sum.tv<-data.frame(Type=c("all", type)) # for F estimates
CI.sum.tv$SE[1]<-std.error(mftv$mean, na.rm=T)
CI.sum.tv$CI[1]<-std.error(mftv$mean, na.rm=T)*1.96

type<-c("syn","nonsyn","stop")

for (i in 1: length(type)){
    pos.b<-mftv$pos[mftv$Type.tv1==type[i]&mftv$Type.tv2==type[i]]
    pos.1<-mftv$pos[mftv$Type.tv1==type[i]&mftv$Type.tv2!=type[i]]
    pos.2<-mftv$pos[mftv$Type.tv1!=type[i]&mftv$Type.tv2==type[i]] 
    
    df.b<-mftv[mftv$pos %in% pos.b,c("pos","mean")]
    df.1<-mf1[mf1$pos %in% pos.1,c("pos","mean")]
    df.2<-mf2[mf2$pos %in% pos.2,c("pos","mean")]
    
    df<-rbind(df.b, df.1, df.2)
    
    mf.vec<-c()
    mf.vec<-c(mf.vec, mftv$mean[mftv$pos %in% pos.b])
    mf.vec<-c(mf.vec, mf1$mean[mf1$pos %in% pos.1])
    mf.vec<-c(mf.vec, mf2$mean[mf2$pos %in% pos.2])
    
    #write.csv(mf.vec, paste0("Output/MutFreq/MF_tvs.",type[i],".csv"))
    
    assign(paste0("Tvs.", type[i]), mf.vec)
    
    #average mut freq
    summary.mf[4,type[i]]<-mean(df$mean,na.rm=T)
    
    #get average CIs
    CI.sum.tv[(i+1),"SE"]<-std.error(df$mean, na.rm=T)
    CI.sum.tv[(i+1),"CI"]<- std.error(df$mean, na.rm=T)*1.96
}

write.csv(summary.mf, "Output/MutFreq/MF_summary.csv")

tvs.mf<-data.frame(freq=Tvs.syn)
tvs.mf$Type<-"Syn"
tvs.mf<-rbind(tvs.mf, data.frame(freq=Tvs.nonsyn, Type=rep("Nonsyn", times=length(Tvs.nonsyn))))
tvs.mf<-rbind(tvs.mf, data.frame(freq=Tvs.stop, Type=rep("Nonsense", times=length(Tvs.stop))))
tvs.mf$Mutation<-"Tvs"
write.csv(tvs.mf, "Output/MutFreq/MF_tvs_all.csv")


#Combine the tv and tvs and two methods
CI.sum.tv$Mutation<-"Tvs"

CI.summary<-rbind(CI.sum, CI.sum.tv)

write.csv(CI.summary,"Output/MutFreq/SE_TypeofMutaions_summary.csv")



### Structural vs. non-structural genes mut. freq
TS<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)

st<-TS[TS$pos>341& TS$pos<=2579,]
nonst<-TS[TS$pos>2579,]

mean(st$mean, na.rm=T) #0.004964672
mean(nonst$mean, na.rm = T)  #0.004777546

r1<-wilcox.test(st$mean,nonst$mean, alternative = "greater", paired = FALSE) 
r1[[3]]  #P=0.02269177




