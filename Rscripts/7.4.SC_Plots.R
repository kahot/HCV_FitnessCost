# Plot the selection coefficient estimation results
library(ggplot2)
library(reshape2)
library(sfsmisc)
library(colorspace)
library(cowplot)
source("Rscripts/label_scientific.R")


colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")

###########
df<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
#coding regions only
df<-df[df$pos>=342,]

#######
#Plot summary of SC by type by nucleotide (Fig 4B)      
#2 no CPG creating mutations
sc3<-sc[sc$makesCpG==0,]
k=1
transSC<-list()
for (i in c("a","t","c","g")) {
    for (type in c("syn","nonsyn")){
        datavector<-sc3$EstSC[sc3$Type==type & sc3$ref==i]
        nt<-toupper(i)
        vname<-paste0(nt,".",type)
        dat<-data.frame(base=rep(nt,times=length(datavector)),
                        type=rep(type, times=length(datavector)), S.C.=datavector)
        transSC[[k]]<-dat
        names(transSC)[k]<-vname
        k=k+1
    }
}

scdata2<-do.call(rbind, transSC)

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 

scdata2$base<-factor(scdata2$base, levels=c("A","T", "C", "G"))
scdata2$type<-factor(scdata2$type, levels=c("syn","nonsyn"))

p1<-ggplot(scdata2,aes(x=base, y=S.C., fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
    scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
    labs(x="Nucleotide",y="Estimated selection coefficient")+
    scale_fill_manual(values=colors2[c(5,1)]) + theme_bw()+
    theme(legend.title = element_blank())+theme(axis.text.x = element_text(size =10, color=1), axis.title.x = element_blank())+
    theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
    geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
    scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))

#Add the mutation frequency summary plot for contrast (Fig. 4C)
Ts <-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)
#coding region only
Ts<-Ts[Ts$pos>=342,c("pos","ref","mean","Type","makesCpG") ]
Ts<-Ts[Ts$makesCpG==0,]

#Plot summary of mutation frequency by type by nucleotide       
k=1
transMF<-list()
for (i in c("a","t","c","g")) {
    for (type in c("syn","nonsyn")){
        datavector<-Ts$mean[Ts$Type==type & Ts$ref==i]
        nt<-toupper(i)
        vname<-paste0(nt,".",type)
        dat<-data.frame(base=rep(nt,times=length(datavector)),
                        type=rep(type, times=length(datavector)), MF=datavector)
        transMF[[k]]<-dat
        names(transMF)[k]<-vname
        k=k+1
    }
}

mfdata<-do.call(rbind, transMF)
mfdata$base<-factor(mfdata$base, levels=c("A","T", "C", "G"))
mfdata$type<-factor(mfdata$type, levels=c("syn","nonsyn"))

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 
col2_light1<-qualitative_hcl(6, palette="Set2")

p2<-ggplot(mfdata,aes(x=base, y=MF, fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
    scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
    labs(x="",y="Mutation frequency")+
    scale_fill_manual(values=col2_light1[c(5,1)]) + theme_bw()+
    theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10, color=1))+
    theme(axis.text.y = element_text(size =10), axis.title.y= element_text(size =12), panel.grid.major.x=element_blank())+
    geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
    scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))


pdf("Output/SelCoeff/SC-MF-boxplot.pdf",width=12, height=4)
ggdraw()+
    draw_plot(p1,x=0,y=0.05,width=0.45, height=0.9)+
    draw_plot(p2,x=0.5,y=0,width=0.45, height=0.95)+
    draw_plot_label(c("B", "C"), c(0, 0.5), c(1, 1), size = 18,fontface = "plain")
dev.off()    

