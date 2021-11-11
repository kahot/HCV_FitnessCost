# Create mutation frequency summary datatbles and estimate SE/CI

#Read the filtered overview files
Files<-list.files("Output/Overview3/",pattern="overview3.csv")
Overviews<-list()
for (i in 1:length(Files)){ 
        overviews2<-read.csv(paste0("Output/Overview3/",Files[i]),stringsAsFactors=FALSE, row.names=1)
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

#Extract metadata from one file
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
#Summarize the mean and SE for all types of mutations

## First get the read depths for all positions
Reads<-data.frame(sapply(Overviews, "[[", "TotalReads"))
df<-Overviews[[1]]
Reads<-cbind(df[,"pos"],Reads)
colnames(Reads)[1]<-"pos"
write.csv(Reads, "Output/ReadDepth_All.csv")

#Mean read depth per site
Reads$Depth<-rowMeans(Reads[2:196], na.rm=T)
depth<-Reads[,c("pos","Depth")]
write.csv(depth,"Output/ReadDepth_mean.csv")

#Calculate SE for all positions using sqrt(p(1-p)/n)
reads<-Reads[Reads$pos>341&Reads$pos<8575,]
SE<-list()
for (f in 1:6){
        mf<-mf.files[[f]]
        mf<-mf[mf$pos>341&mf$pos<8575,]
        reads2<-reads[reads$pos %in% mf$pos,]
        
        df<-data.frame(pos=mf$pos)
        for (i in 1:195){
                mf2<-mf[,i]
                re<-reads2[,(i+1)]
                df[,colnames(mf)[i]]<-sqrt(mf2*(1-mf2)/re)
        }
        write.csv(df,paste0("Output/MutFreq/SEmatrix_",names(mf.files)[f],".csv"))
        SE[[f]]<-df
        names(SE)[f]<-names(mf.files)[f]
        
}

#save(SE, file="Output/SE.Rdata")
#load("Output/SE.Rdata")

### Calculate the average statistics for all types of mutations
tb3<-data.frame(type=c("Ts","Ts.syn","Ts.ns","Ts.stop", "Tv1","Tv1.syn","Tv1.ns","Tv1.stop","Tv2","Tv2.syn","Tv2.ns","Tv2.stop","Tvs","All" ))
for (i in 1:5){
        se<-SE[[i]]
        se$mean<-rowMeans(se[,2:196],na.rm = T)
        dt<-mf.files[[i]]
        #coding region only
        dt<-dt[dt$pos>341&dt$pos<8575,]
        
        #For transition mutations
        if (i<=3){
                colnames(dt)[199]<-"Type"
                k<-(i-1)*4+1  #k=>specify the row for the summary output
                tb3$mean[k]<-mean(dt$mean,na.rm=T)
                tb3$se[k]<-mean(se$mean, na.rm = T)
                k<-k+1
                tb3$mean[k]<-mean(dt$mean[dt$Type=="syn"],na.rm=T)
                tb3$se[k]<-mean(se$mean[dt$Type=="syn"], na.rm = T)
                k<-k+1
                tb3$mean[k]<-mean(dt$mean[dt$Type=="nonsyn"],na.rm=T)
                tb3$se[k]<-mean(se$mean[dt$Type=="nonsyn"], na.rm = T)
                k<-k+1
                tb3$mean[k]<-mean(dt$mean[dt$Type=="stop"],na.rm=T)
                tb3$se[k]<-mean(se$mean[dt$Type=="stop"], na.rm = T)
        }
        #Tvs
        if (i==4) {k=13
        tb3$mean[k]<-mean(dt$mean,na.rm=T)
        tb3$se[k]<-mean(se$mean, na.rm = T) }
        if (i==5) {k=14
        tb3$mean[k]<-mean(dt$mean,na.rm=T)
        tb3$se[k]<-mean(se$mean, na.rm = T) }
}


addtb3<-data.frame(type=c("tvs.syn","tvs.ns","tvs.stop", "all.syn","all.ns","all.stop"),
                   mean=c(mean(tb3$mean[6],tb3$mean[10]),mean(tb3$mean[7],tb3$mean[11]),mean(tb3$mean[8],tb3$mean[12]),
                          mean(tb3$mean[2],tb3$mean[6],tb3$mean[10]),mean(tb3$mean[3],tb3$mean[7],tb3$mean[11]),mean(tb3$mean[4],tb3$mean[8],tb3$mean[12])),
                   se  =c(mean(tb3$se[6],tb3$se[10]),mean(tb3$se[7],tb3$se[11]),mean(tb3$se[8],tb3$se[12]),
                          mean(tb3$se[2],tb3$se[6],tb3$se[10]),mean(tb3$se[3],tb3$se[7],tb3$se[11]),mean(tb3$se[4],tb3$se[8],tb3$se[12])))

TB3<-rbind(tb3, addtb3)

write.csv(TB3, "Output/MutFreq/MF.Mean.SE.summary.csv")


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

######
######
#Calculate CI
#Using SE and z-score (alpha=0.95)
TB3$CI_up  <-TB3[,"mean"]+TB3[,"se"]*1.96
TB3$CI_down<-TB3[,"mean"]-TB3[,"se"]*1.96
TB3$CI<-TB3[,"se"]*1.96

write.csv(TB3,"Output/MutFreq/MF.Mean.SE.summary.csv")
