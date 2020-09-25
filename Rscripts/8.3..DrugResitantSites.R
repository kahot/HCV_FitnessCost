library(tidyverse)
source("Rscripts/baseRscript.R")
colors=c("#44AA99","#0077BB","#CC6677" )

####### START ####
#read the modified table of RAV info
DR<-read.csv("Data/RAV_Table_updated.csv",stringsAsFactors = F)
#create an id column
for (i in 1:nrow(DR)){
        if (DR$Need_both[i]=="y") { 
                if (DR$extra[i]=="n") DR$ID[i]<- paste0(DR$Name[i],'.',DR$merged.pos[i])
                if (DR$extra[i]=="y") DR$ID[i]<- paste (DR$Name[i],DR$merged.pos[i],DR$Type[i], sep = ".")}
        if (DR$Need_both[i]=="n") { 
                if (DR$extra[i]=="n") DR$ID[i]<-DR$Name[i]
                if (DR$extra[i]=="y") DR$ID[i]<- paste (DR$Name[i],DR$merged.pos[i],DR$Type[i], sep = ".")}
}

HCVFiles<-list.files("Output1A/Overview2/", pattern="overview2.csv")

DR_mutfreq<-data.frame(ID=factor(DR$ID, levels=c(DR$ID)))
Diff<-data.frame(ID=factor(DR$ID,levels=c(DR$ID)))
diff.count<-list()

for (i in 1:length(HCVFiles)){ 
        df<-read.csv(paste0("Output1A/Overview2/",HCVFiles[i]),stringsAsFactors=FALSE, row.names = 1)
        dname<-substr(paste(HCVFiles[i]),start=1,stop=7)
        dr<-DR
        cname<-"pos.1A"
        DRsites<-df[df$pos %in% dr[,cname],]
        
        #count the number of samples fixed with the RAVs
        dr$obs<-0
        for (k in 1:nrow(dr)){
                pos<-dr[k,cname]
                if (is.na(DRsites$MajNt[DRsites$pos==pos])|is.na(DRsites$ref[DRsites$pos==pos])) next
                if (DRsites$MajNt[DRsites$pos==pos]!=DRsites$ref[DRsites$pos==pos]){
                        if (dr$Type[k]=="Ts") {mutnt<-transition(DRsites$ref[DRsites$pos==pos]) 
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1} 
                        if (dr$Type[k]=='Tv1') {mutnt<-transv1(DRsites$ref[DRsites$pos==pos])
                                        if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1}
                        if (dr$Type[k]=='Tv2') {mutnt<-transv2(DRsites$ref[DRsites$pos==pos])
                                 if (DRsites$MajNt[DRsites$pos==pos]==mutnt) dr$obs[k]<-1} 
                }
                        
        }
        
        #store the observation number in a list
        diff.count[i]<-sum(dr$obs==1)
        names(diff.count)[i]<-dname
        
        #obtain the mutation freq.
        dr$freq<-""
        for (k in 1:nrow(dr)){
                if (dr$Type[k]=='Tv1')    dr$freq[k]<-DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]]
                if (dr$Type[k]=='Tv2')    dr$freq[k]<-DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]]
                if (dr$Type[k]=='Ts')     dr$freq[k]<-DRsites$freq.Ts.ref     [DRsites$pos==dr[k,cname]]
        }
        
        #write.csv(dr, paste0("Output_all/DR/",geno[j],"/DRsites.",dname,".csv"))
        dr$freq<-as.numeric(dr$freq)
        
        #mut.freq
        dr1<-dr[,c("ID","freq")]
        colnames(dr1)[2]<-dname
        DR_mutfreq<-merge(DR_mutfreq,dr1, by="ID")   
        
        # no of samples fixed to RAV 
        dr2<-dr[,c("ID","obs")]
        colnames(dr2)[2]<-dname
        Diff<-merge(Diff, dr2, by="ID")
}

DR.mutated.counts<-do.call(rbind,diff.count)

DR_mutfreq<-DR_mutfreq[order(DR_mutfreq$ID),]
DR_diff<-Diff[order(Diff$ID),]

#count the number of patients fixed with RAV (%)
#count the number of non-NA per RAV
DR_diff$NonNA_count<-apply(DR_mutfreq[,2:(s+1)],1, function(x) sum(!is.na(x)))
DR_diff$total<-apply(DR_diff[2:(s+1)],1,sum, na.rm=T)
DR_diff$Percent<-format(round(DR_diff$total/DR_diff$NonNA_count*100, 1), nsmall=1)
write.csv(DR_diff, "Output1A/DrugRes/RAV.counts.MutFreq_summary.1A.csv")
write.csv(DR_mutfreq, "Output1A/DrugRes/RAV.MutationFreq_summary.1A.csv")


####  START here for plotting only ####
DR_mutfreq<-read.csv("Output1A/DrugRes/RAV.MutationFreq_summary.1A.csv", stringsAsFactors = F, row.names = 1)
DR_diff<-read.csv("Output1A/DrugRes/RAV.MutationFreq_summary_Diff.1A.csv", stringsAsFactors = F, row.names = 1)

######
s<-length(HCVFiles)
ns3n <-nrow(DR[DR$Gene=="NS3",])
ns5an<-nrow(DR[DR$Gene=="NS5A",])
ns5bn<-nrow(DR[DR$Gene=="NS5B",])

#####
#count the # samples that showed drug resistant mutations
DR_mutfreq$Count<-apply(DR_mutfreq[,2:196],1, function(x) sum(x!=0, na.rm=T))
DR_mutfreq$NonNA_count<-apply(DR_mutfreq[,2:(s+1)],1, function(x) sum(!is.na(x)))
DR_mutfreq$Percent<-format(round(DR_mutfreq$Count/DR_mutfreq$NonNA_count*100, 1), nsmall=1)
NatPrev<-DR_mutfreq[,c('ID','Count','Percent')]

write.csv(NatPrev,"Output1A/DrugRes/RAV.NatPrevelence.count.1A.csv")


#####
## Calculate selection coefficients for RAV sites
## this inlucdes transversion so can't use the SCs summary

HCVFiles3<-list.files("Output1A/Overview3/", pattern="overview3.csv")

#read the modified table of RAV info
dr1a<-read.csv("Output1A/DrugRes/RAV_Table_updated.csv",stringsAsFactors = F)

#create an id column
for (i in 1:nrow(dr1a)){
        if (dr1a$Need_both[i]=="y") { 
                if (dr1a$extra[i]=="n") dr1a$ID[i]<- paste0(dr1a$Name[i],'.',dr1a$merged.pos[i])
                if (dr1a$extra[i]=="y") dr1a$ID[i]<- paste (dr1a$Name[i],dr1a$merged.pos[i],dr1a$Type[i], sep = ".")}
        if (dr1a$Need_both[i]=="n") { 
                if (dr1a$extra[i]=="n") dr1a$ID[i]<-dr1a$Name[i]
                if (dr1a$extra[i]=="y") dr1a$ID[i]<- paste (dr1a$Name[i],dr1a$merged.pos[i],dr1a$Type[i], sep = ".")}
}


dr.mf2<-data.frame(ID=dr1a$ID)
for (i in 1:length(HCVFiles3)){ 
        df<-read.csv(paste0("Output1A/Overview3/",HCVFiles3[i]),stringsAsFactors=FALSE, row.names = 1)
        dname<-substr(paste(HCVFiles3[i]),start=1,stop=7)
        dr<-dr1a
        cname<-"pos.1A"
        DRsites<-df[df$pos %in% dr[,cname],]
        
        #obtain the mutation freq.
        dr$freq<-""
        for (k in 1:nrow(dr)){
                if (dr$Type[k]=='Tv1')    dr$freq[k]<-DRsites$freq.transv1.ref[DRsites$pos==dr[k,cname]]
                if (dr$Type[k]=='Tv2')    dr$freq[k]<-DRsites$freq.transv2.ref[DRsites$pos==dr[k,cname]]
                if (dr$Type[k]=='Ts')     dr$freq[k]<-DRsites$freq.Ts.ref     [DRsites$pos==dr[k,cname]]
        }
        
        dr$freq<-as.numeric(dr$freq)
        dr<-dr[,c("ID","freq")]
        colnames(dr)[2]<-dname
        dr.mf2<-merge(dr.mf2,dr, by="ID")
}

write.csv(dr.mf2, "Output1A/DrugRes/RAV.MF_summary.1A.Filtered.csv")


dr.mf2$mean<-rowMeans(dr.mf2[,2:196],na.rm=T)
dr_sc<-dr.mf2[,c("ID","mean")]
dr_sc<-merge(dr_sc, dr1a, by="ID")

# attach the mut rates info:
mutrates<-read.csv("Output1A/Geller/Geller.MutRates.Summary_updated.csv", row.names = 1, stringsAsFactors = F)

dr_sc$MR[dr_sc$ref=="c"&dr_sc$Type=="Ts"] <-mutrates$Mean[mutrates$Mutation=="CU"]
dr_sc$MR[dr_sc$ref=="a"&dr_sc$Type=="Ts"] <-mutrates$Mean[mutrates$Mutation=="AG"]
dr_sc$MR[dr_sc$ref=="g"&dr_sc$Type=="Ts"] <-mutrates$Mean[mutrates$Mutation=="GA"]
dr_sc$MR[dr_sc$ref=="t"&dr_sc$Type=="Ts"] <-mutrates$Mean[mutrates$Mutation=="UC"]
dr_sc$MR[dr_sc$ref=="a"&dr_sc$Type=="Tv1"]<-mutrates$Mean[mutrates$Mutation=="AC"]
dr_sc$MR[dr_sc$ref=="c"&dr_sc$Type=="Tv1"]<-mutrates$Mean[mutrates$Mutation=="CA"]
dr_sc$MR[dr_sc$ref=="g"&dr_sc$Type=="Tv1"]<-mutrates$Mean[mutrates$Mutation=="GC"]
dr_sc$MR[dr_sc$ref=="t"&dr_sc$Type=="Tv1"]<-mutrates$Mean[mutrates$Mutation=="UA"]
dr_sc$MR[dr_sc$ref=="a"&dr_sc$Type=="Tv2"]<-mutrates$Mean[mutrates$Mutation=="AU"]
dr_sc$MR[dr_sc$ref=="c"&dr_sc$Type=="Tv2"]<-mutrates$Mean[mutrates$Mutation=="CG"]
dr_sc$MR[dr_sc$ref=="g"&dr_sc$Type=="Tv2"]<-mutrates$Mean[mutrates$Mutation=="GU"]
dr_sc$MR[dr_sc$ref=="t"&dr_sc$Type=="Tv2"]<-mutrates$Mean[mutrates$Mutation=="UG"]


dr_sc$EstSC<-""
for (j in 1:nrow(dr_sc)){
        dr_sc$EstSC[j] <- EstimatedS(dr_sc$MR[j],dr_sc$mean[j])
}
dr_sc$EstSC<-as.numeric(dr_sc$EstSC)
dr_sc$ID<-factor(dr_sc$ID, levels = c(DR$ID))
dr_sc<-dr_sc[order(dr_sc$ID),]
write.csv(dr_sc, "Output1A/DrugRes/RAV.SC_summary.1A.csv")


### Create a figure ###

dr_sc<-read.csv("Output1A/DrugRes/RAV.SC_summary.1A.csv")

pdf(paste0("Output1A/DrugRes/RAVs.sc.plots.pdf"), height = 11, width = 15.5)
layout(matrix(c(1,2), byrow = TRUE), heights=c(1,3))
par(mar=c(0, 4, 4, .5))

genesDF<-data.frame("name"= c("NS3","NS5A","NS5B"), "Begin"= c(1,ns3n+1,(ns3n+ns5an+1)),"End"= c(ns3n,(ns3n+ns5an),nrow(DR_mutfreq)))

ymin1=-4
plot(0, type = "n", xlim = c(1, nrow(dr_sc)), ylim = c(ymin1, -1), axes = FALSE, ylab = "Cost", xlab = "")
axis(side = 2, at = seq(-1, ymin1, by=-1), labels = expression(10^-1, 10^-2, 10^-3, 10^-4), las = 2)
abline(v=c((0:nrow(dr_sc))+0.5), col="gray60", lty=1, lwd=.1)
#abline(h=-4, col="gray50")
segments(i,-3.9,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[1], lend=2)

for (i in 1:nrow(dr_sc)){
        if (i<=genesDF$End[1]) segments(i,ymin1,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[1],lend=2)
        if (i<=genesDF$End[2]& i>=genesDF$Begin[2]) segments(i,ymin1,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[2],lend=2)
        if (i<=genesDF$End[3]& i>=genesDF$Begin[3]) segments(i,ymin1,i,log10(dr_sc$EstSC[i]), lwd=10, col=colors[3],lend=2)
}

abline(v=ns3n+.5,col='gray50', lwd=3)
abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )

par(mar=c(6, 4, 4, .5))
ymin <- -5
plot(0, type = "n", xlim = c(1, nrow(DR_mutfreq)), ylim = c(ymin, 0), axes = FALSE, ylab = "Frequency of resistance-assocaited variants",xlab = "")
axis(side = 2, at = seq(0, ymin, by=-1), labels = expression(10^0, 10^-1, 10^-2, 10^-3, 10^-4,10^-5), las = 2)
n<-seq(1, by = 2, len = (nrow(DR_mutfreq)/2))

for (i in 1:3){
        nvec<-seq(genesDF$Begin[i],genesDF$End[i],2)
        nvec2<-seq(genesDF$Begin[i]+1,genesDF$End[i],2)
        
        if (nvec[1]%%2==0) {v2<-nvec; v1<-nvec2}
        if (nvec[1]%%2==1) {v1<-nvec; v2<-nvec2}
        
        abline(v=v1, col=paste0(colors[i],"66"),lty=1, lwd=16)
        abline(v=v2, col=paste0(colors[i],"1A"),lty=1, lwd=16)
}

for (i in 1:nrow(DR_mutfreq)){
        if (i<=genesDF$End[1]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[1])}
        if (i<=genesDF$End[2]& i>=genesDF$Begin[2]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.4,col=colors[2])}
        if (i<=genesDF$End[3]& i>=genesDF$Begin[3]){
                xjit <- rnorm(s, 0, .1)
                points(rep(i,s)+xjit, log10(DR_mutfreq[i,2:(s+1)]),pch=16,cex=0.3, col=colors[3])}
}

# add % of fixed samples
for (i in 1:nrow(DR)){
        if (DR$Type[i]=='Ts') mtext(side = 3, at = i, text = paste(DR_diff$Percent[i]), padj=0, las=2, cex = .8)
        if (DR$Type[i]=='Tv1'|DR$Type[i]=='Tv2')  mtext(side = 3, at = i, text = paste(DR_diff$Percent[i]), las=2,  padj=0, cex = .8, col="#9437FF")
}


#add IDs
for (i in 1:nrow(DR)){
        if (DR$Type[i]=='Ts') mtext(side = 1, at = i, text = paste(dr$ID[i]), las=2, padj=0, cex = .8)
        if (DR$Type[i]=='Tv1'|DR$Type[i]=='Tv2')  mtext(side = 1, at = i, text = paste(dr$ID[i]), las=2, padj=0, cex = .8, col="#9437FF")
}

#add % of samples with RAV observed 
for (i in 1:nrow(DR)){
        if (DR$Type[i]=='Ts') mtext(side = 3, at = i, line=2, text = paste(NatPrev$Percent[i]), las=2, padj=0, cex = .8)
        if (DR$Type[i]=='Tv1'|DR$Type[i]=='Tv2')  mtext(side = 3, at = i, line=2, text = paste(NatPrev$Percent[i]), las=2, padj=0, cex = .8, col="#9437FF")
}
abline(v=ns3n+.5,col='gray50', lwd=3)
abline(v=ns3n+ns5an+.5,col='gray50', lwd=3 )

rect(0.5,-5,ns3n+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[1])
text(ns3n/2+.5,-5.08,paste0(genesDF$name[1]),col="black", cex=.8)
rect(ns3n+.5,-5, ns3n+ns5an+.5 ,-5.19,density = NULL, angle = 45,col="white",border =colors[2])
text(ns3n+ns5an-ns5an/2+.5,-5.08,paste0(genesDF$name[2]),col="black", cex=.8)
rect(ns3n+ns5an+.5,-5,nrow(dr)+.5,-5.19,density = NULL, angle = 45,col="white",border =colors[3])
text(nrow(dr)-(nrow(dr)-ns3n-ns5an)/2+.5,-5.08,paste0(genesDF$name[3]),col="black", cex=.8)

dev.off()

## Natural Prevelence vs. fixed observations
natprev<-NatPrev
natprev$Percent<-as.numeric(natprev$Percent)
natprev$Percent_fixed<-DR_diff$Percent
natprev$Percent_fixed<-as.numeric(natprev$Percent_fixed)


range(natprev$Percent)
# 3.6 99.5
mean(natprev$Percent, na.rm=T)
#61.55057
mean(natprev$Percent_fixed, na.rm=T)
#1.30%

cor.test(natprev$Percent, natprev$Percent_fixed,method = "spearman")
#p-value = 0.0002183
#rho=0.3864204



##### SC comparison between RAV sites vs. non-RAV sites
#Read the summary RAV Sel coeff table
dr_sc<-read.csv("Output1A/DrugRes/RAV.SC_summary.1A.csv", row.names = 1, stringsAsFactors = F)
drsites<-dr_sc[,c("ID","pos.1A","EstSC","Type","Gene","MR")]
drsitesT<-dr_sc[dr_sc$Type=="Ts",]

#get the non-DRsites sel coeff info
sc<-read.csv("Output1A/SelCoeff/SC.csv", row.names = 1, stringsAsFactors = F)
sc<-sc[sc$pos>341,]

### addd gene annotation info for all genes
genes<-read.csv("Data/HCV_annotations2.csv")
genenames<-genes$Gene
gene<-c()
for (i in 2:12){
    gene<-c(gene, rep(i, times=genes$end[i]-genes$start[i]+1))
}

n<-data.frame(pos=342:(length(gene)+341))
g<-cbind(n,gene)

sc<-merge(sc, g, by ="pos")

drpos<-drsitesT$pos.1A 
NonDRsites<-sc[!(sc$pos %in% drpos),]

#average SC of drug resistance sites 
mean(drsites$EstSC) 
#0.004722373
#average Transition SC sites
mean(drsitesT$EstSC) 
#0.002482708
mean(dr_sc$EstSC[dr_sc$Type=='Tv1'|dr_sc$Type=='Tv2'], na.rm=T)
#0.006812726

#Genome wide average SC
mean(sc$EstSC)
#0.001977764

#genome wide average SC of non-drug resisitant sites
mean(NonDRsites$EstSC)
#0.001975877

# Average of nonDR sites for the 3 genes
mean(NonDRsites$EstSC[NonDRsites$gene==8|NonDRsites$gene==11|NonDRsites$gene==12])
# 0.002052286
std.error(NonDRsites$EstSC[NonDRsites$gene==8|NonDRsites$gene==11|NonDRsites$gene==12])


#Test the difference
wilcox.test(drsitesT$EstSC, NonDRsites$EstSC[NonDRsites$gene==8|NonDRsites$gene==11|NonDRsites$gene==12], alternative = "greater", paired = FALSE)
#p-value = 0.002332
wilcox.test(drsitesT$EstSC[drsitesT$Gene=="NS3"], NonDRsites$EstSC[NonDRsites$gene==8], alternative = "greater", paired = FALSE)
#p-value = 0.077
wilcox.test(drsitesT$EstSC[drsitesT$Gene=="NS5A"], NonDRsites$EstSC[NonDRsites$gene==11], alternative = "greater", paired = FALSE)
#p-value = 0.01897
wilcox.test(drsitesT$EstSC[drsitesT$Gene=="NS5B"], NonDRsites$EstSC[NonDRsites$gene==12], alternative = "greater", paired = FALSE)
#p-value = 0.03766

