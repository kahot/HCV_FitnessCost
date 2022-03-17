vl<-read.csv("/Users/kahotisthammer/Projects/HCV_mutation_patterns/Data/Viral_load.csv")
vl<-vl[vl$Subtype=="1A",]

#viral load vs. primer ID generated copies comparisons

siv.vl<-read.csv("/Users/kahotisthammer/Projects/SIV-R01/Data/ViralLoads.csv")
siv.reads<-read.csv("/Users/kahotisthammer/Projects/SIV-R01/Output/ReadDeapth_all.csv")
samples<-read.csv("/Users/kahotisthammer/Projects/SIV-R01/Data/SamplesNoduplicates.csv")
samples<-samples[samples$Tissue=="Plasma",]
siv.vl$id<-paste0(siv.vl$Monkey,".",siv.vl$Week)
samples$id<-paste0(gsub("A",'',samples$Monkey),".",samples$Week)
siv.vl<-merge(siv.vl, samples[,c("File.name", "id")], by="id")

siv<-merge(siv.vl, siv.reads[c("File.name","Average")], by="File.name")

siv$ratio<-siv$Average/siv$VL
mean(siv$ratio)*100
#0.01697388 %