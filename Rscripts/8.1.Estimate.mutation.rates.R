library(tidyverse)
library(zoo)
library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(sfsmisc)
library(colorspace)
source("Rscripts/label_scientific.R")
source("Rscripts/baseRscript.R")

colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])
color.genes<-qualitative_hcl(11, palette="Dark3")


TS<-read.csv("Output/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)
stops<-TS[TS$Type=="stop",]
#transition rates (only c and g are present)
stopC<-stops[stops$ref =="c",]
stopG<-stops[stops$ref =="g",]

#G to A mutation
mean(stopG$mean) #0.001915474
#C to T mutation
mean(stopC$mean) #0.002165682


Tv1<-read.csv("Output/MutFreq.filtered/Filtered.Tv1.MutFreq.Q35.csv",stringsAsFactors = F, row.names = 1)
Tv2<-read.csv("Output/MutFreq.filtered/Filtered.Tv2.MutFreq.Q35.csv",stringsAsFactors = F, row.names = 1)

stopTv1<-Tv1[Tv1$Type.tv1=="stop",]
stopTv2<-Tv2[Tv2$Type.tv2=="stop",]

table(stopTv1$ref)
#c   t 
#175  97

table(stopTv2$ref)
#a   c   g   t 
#107  80 138  36 


#C to A mutation
mean(stopTv1$mean[stopTv1$ref=="c"]) #0.0005057763
#T to A mutation
mean(stopTv1$mean[stopTv1$ref=="t"]) #0.0005727744

#A to T mutation
mean(stopTv2$mean[stopTv2$ref=="a"]) #0.0006335722
#C to G mutation
mean(stopTv2$mean[stopTv2$ref=="c"]) # 3.796976e-05
#G to T mutation
mean(stopTv2$mean[stopTv2$ref=="g"]) # 0.001000036
#T to G mutation
mean(stopTv2$mean[stopTv2$ref=="t"]) # 0.0001766808


## mutation rates from Geller
geller<-read.csv("Data/Geller.mutation.rates.csv", stringsAsFactors = F, row.names = 1)
geller<-geller[!(geller$mutation %in% c("AA","CC","GG","UU")),]

Mutrates<-data.frame(mutation=c("CT","GA","CA","TA","AT","CG","GT","TG"),
           Estimated=c(mean(stopG$mean), 
                      mean(stopC$mean),
                      mean(stopTv1$mean[stopTv1$ref=="c"]), 
                      mean(stopTv1$mean[stopTv1$ref=="t"]), 
                      mean(stopTv2$mean[stopTv2$ref=="a"]), 
                      mean(stopTv2$mean[stopTv2$ref=="c"]), 
                      mean(stopTv2$mean[stopTv2$ref=="g"]), 
                      mean(stopTv2$mean[stopTv2$ref=="t"]) ))

geller$mutation<-gsub("U","T",geller$mutation)

Mutrates<-merge(Mutrates, geller[,1:2], by="mutation", all.y = T )
colnames(Mutrates)[2:3]<-c("In vivo estimation", "In vitro")

#order based on mut rates
Mutrates<-Mutrates[order(Mutrates$`In vitro`, decreasing = T),]
Mutrates$mutation<-factor(Mutrates$mutation, levels=paste0(Mutrates$mutation))                           

MutratesM<-melt(Mutrates)
ggplot(MutratesM, aes(x= mutation, y=value, color=variable ))+
        geom_point(size=2)+
        scale_y_continuous(trans = "log",labels=label_scientific, breaks=c(10^-7, 10^-6,10^-5, 10^-4, 10^-3))+
        scale_color_manual(values=colors2[c(1,4)])+
        theme_bw()+
        ylab("Mutation rate")+xlab('Mutation')+
        theme(legend.title = element_blank(), panel.grid.minor.y  = element_blank())+
        scale_x_discrete(breaks=c("CG","GC","CA","TA","GT","TG","AT","AC","GA","CT","TC","AG"),labels=c(expression(C%->%G),expression(G%->%C),expression(C%->%A),expression("T"%->%A),
                                              expression(G%->%"T"),expression("T"%->%G),expression(A%->%"T"),expression(A%->%C),expression(G%->%A),expression(C%->%"T"), expression("T"%->%C), expression(A%->%G)))+
        theme(axis.text.x = element_text(angle = 90))
ggsave("Output/SummaryFig/Estimated_MutRates_invivo-invitro.pdf", width = 6, height = 4)


cor.test(Mutrates$`In vivo estimation`, Mutrates$`In vitro`, method="spearman")
#S = 18, p-value = 0.02793
#rho 0.7857143

#correlation plot
mrates2<-Mutrates[!is.na(Mutrates$`In vivo estimation`),]
colnames(mrates2)[2:3]<-c("In vivo", "In vitro")

ggplot(mrates2, aes(x=`In vivo`, y=`In vitro`))+
        geom_point(size=2, color=colors2[5])+
        theme_bw()+
        geom_smooth(method='lm', formula= y~x, se=F, lwd=.2, color="gray50")
ggsave("Output/SummaryFig/Invivo.Invtro.Mutation.rates.correlation.pdf", width = 4.5, height = 4)
