# Calculate selection coefficients using mutation rates estimated for each site (i.e. in vitro mutation frequency/23.3)

library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)
library(cowplot)
scCols<-c("#E16A86","#009ADE")
source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")
###########

#### Site by Site mut rate:
mf<-read.csv("Output/Geller/Geller_mf_overview.csv", row.names = 1, stringsAsFactors = F)
min(mf$mean[mf$mean>0], na.rm = T) #2.780489e-06
nrow(mf[mf$mean==0,]) #1044

#replace 0 with 1*10^-6
mf2<-mf
mf2$mean[mf2$mean==0]<-1*10^-6

mf2$mr<-mf2$mean/23.3

df<-read.csv("Output/Mutfreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names = 1 )
df<-df[,c(196:203)]

df.ind<-merge(df, mf2[,c("pos","mr")], by="pos")
df.ind$EstSC<-apply(df.ind[c("mr","mean")],1, function(x) EstimatedS(x[1],x[2]))
df.ind$EstSC<-as.numeric(df.ind$EstSC)
df<-df.ind
mean(df$EstSC, na.rm = T) # 0.002026144 vs. 0.002036205 (single mutrate for each nt)

depth<-read.csv("Output/ReadDepth_sum.csv", stringsAsFactors = F, row.names = 1)
df<-merge(df,depth, by="pos", all.x=T)
#write.csv(df,"Output/SelCoeff/SC.csv") #add the rad depth to SC files
df<-df[!is.na(df$EstSC),]


#### 
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-df$pos[nrow(df)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]

sc<-merge(df, genetable, by="pos")


#### Histogram comparison
SC<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
colnames(SC)[which(colnames(SC)=="EstSC")]<-"SC_original"
sc<-merge(sc, SC[,c("pos","SC_original")], by="pos")


scDF<-sc[,c("SC_original","EstSC")]
range(scDF$SC_original) #0.0001416349 0.0101180076
range(scDF$EstSC, na.rm = T) # 1.63639e-06 2.19968e-01

colnames(scDF)<-c("Uniform", "Site-level")
scD<-melt(scDF)
colnames(scD)[1]<-c("SC")
library(gridExtra)

p1<-ggplot(scD[scD$SC=="Uniform",], aes(x=value))+ggtitle("Single mut rate")+geom_histogram(color="gray60", alpha = 0.5)+theme_bw()+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.000001,.5))+xlab("selection coefficient")
p2<-ggplot(scD[scD$SC=="Site-level",], aes(x=value))+ggtitle("Site-level mut rate")+geom_histogram(color="gray60", alpha = 0.5)+theme_bw()+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.000001,.5))+xlab("selection coefficient")

pdf("Output/SelCoeff_indivMutRates/SC_comparison_all.pdf", width = 4, height = 4)
grid.arrange(p1,p2,nrow=2)
dev.off()

