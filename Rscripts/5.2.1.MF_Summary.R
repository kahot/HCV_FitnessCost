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
reads<-Reads[Reads$pos>341&Reads$pos<8575,]



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



### Structural vs. non-structural genes mut. freq
TS<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)

st<-TS[TS$pos>341& TS$pos<=2579,]
nonst<-TS[TS$pos>2579,]

mean(st$mean, na.rm=T) #0.004964672
mean(nonst$mean, na.rm = T)  #0.004777546

r1<-wilcox.test(st$mean,nonst$mean, alternative = "greater", paired = FALSE) 
r1[[3]]  #P=0.02269177

