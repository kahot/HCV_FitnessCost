# Examine 'conserved sites' vs. 'variable sites' across all sampeles in 1A samples

library(plotrix)
library(reshape)
library(tidyverse)
library(zoo)
library(purrr)
library(ape)
library(seqinr)
library(DescTools)
library(ggplot2)
library(colorspace)
library(emojifont)
source("Rscripts/baseRscript.R")
colors2<-qualitative_hcl(6, palette="Dark3")


align<-read.fasta("Data/HCV1A_Consensus_Alignment.fasta",as.string=TRUE)
consensus<-read.fasta("Data/HCVref.fasta",as.string=TRUE)
consen<-paste(consensus)
consen<-unlist(strsplit(consen, "") )
consen<-consen[289:8613]

al.ta<-data.frame(pos=289:8613)
for (i in 1: length(align)){
        seq<-paste(align[i])
        sname<-substr(names(align)[i], start=1,stop =6)
        seq<-unlist(strsplit(seq, "") )
        al.ta[,sname]<-seq
}


dt<-al.ta[,2:196]
dt2<-dt
dt2[,]<-0
for (i in 1:nrow(dt)){
        for (k in 1: ncol(dt)){
                if (dt[i,k]==consen[i]) dt2[i,k]<-1
        }
}

dt2$Sum<-rowSums(dt2[,1:195])      
write.csv(dt2, "Output/Identical_sites_matrix.csv")

###
dt2<-read.csv("Output/Identical_sites_matrix.csv", row.names = 1, stringsAsFactors = F)
same.sites<-which(dt2$Sum==195)

### Use filtered data first  ###

# separate identical sites and non-identical sites
Ts<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",row.names = 1,stringsAsFactors = F)
Ts<-Ts[Ts$pos>=289& Ts$pos<=8613,]
T_same<-Ts[same.sites,]
T_diff<-Ts[-same.sites,]

mean(T_same$mean, na.rm=T)  #0.004690898
mean(T_diff$mean, na.rm=T) #0.004929429

wilcox.test(T_same$mean, T_diff$mean,  alternative = "less", paired = FALSE)
#W = 7665700, p-value = 0.0003269

nt<-data.frame(table(T_same$ref))
nt2<-data.frame(table(T_diff$ref))

table(T_same$makesCpG)
table(T_diff$makesCpG)

nt$Diff<-nt2$Freq
colnames(nt)[1:2]<-c("nt", "Same")
GTest(nt[,2:3])
#G = 1.8794, X-squared df = 3, p-value = 0.5978


##ad SC
SC<-read.csv("Output/SelCoeff/SC.csv", row.names = 1, stringsAsFactors = F)

SC<-SC[SC$pos>=289& SC$pos<=8613,]
Ts_same<-SC[same.sites,]
Ts_diff<-SC[-same.sites,]

mean(Ts_same$EstSC, na.rm = T) #0.002680869
mean(Ts_diff$EstSC, na.rm = T) #0.002578303

wilcox.test(Ts_same$EstSC, Ts_diff$EstSC,  alternative = "greater", paired = FALSE)
#W = 8210700, p-value = 0.03153

## Coding region only
Ts2<-Ts[Ts$pos>=342& Ts$pos<=8613,]
T_same<-Ts2[same.sites,]
T_diff<-Ts2[-same.sites,]

mean(T_same$mean, na.rm=T)  #0.004779994
mean(T_diff$mean, na.rm=T) #0.00487539
wilcox.test(T_same$mean, T_diff$mean,  alternative = "less", paired = FALSE)
#W = 7837100, p-value = 0.2332

SC2<-SC[SC$pos>=342& SC$pos<=8613,]
Ts_same<-SC2[same.sites,]
Ts_diff<-SC2[-same.sites,]
mean(Ts_same$EstSC, na.rm = T) #0.002621477
mean(Ts_diff$EstSC, na.rm = T) #0.002583721

wilcox.test(Ts_same$EstSC, Ts_diff$EstSC,  alternative = "greater", paired = FALSE)
#W = 8210700, p-value = 0.03736

