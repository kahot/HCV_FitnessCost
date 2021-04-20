library(ggplot2)
library(ggthemes)
library(colorspace)
source("Rscripts/label_scientific.R")
source("Rscripts/baseRscript.R")

colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
div.colors<-c(colors2[1],col2_light[1],colors2[3],col2_light[3],colors2[5],col2_light[5])
color.genes<-qualitative_hcl(11, palette="Dark3")


#Geller's in vitro mut freq
dt2<-read.csv("Output/Geller/Geller_MutFrequency.csv", stringsAsFactors = F, row.names = 1)
dt2<-dt2[!is.na(dt2$mean),]

transmf1<-list()
k<-1
for (i in c("a","t","c","g")) {
        datavector<-dt2$mean[dt2$ref==i]
        nt<-toupper(i)
        dat<-data.frame(NT=rep(nt,times=length(datavector)), mf=datavector)
        transmf1[[k]]<-dat
        names(transmf1)[k]<-nt
        k=k+1
}
mfdata1<-do.call(rbind, transmf1)
mfdata1$Study<-"In vitro"


TS<-read.csv("Output/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)

k=1
transmf2<-list()
for (i in  c("a","t","c","g")) {
        datavector<-TS$mean[TS$ref==i]
        nt<-toupper(i)
        dat<-data.frame(NT=rep(nt,times=length(datavector)), mf=datavector)
        transmf2[[k]]<-dat
        names(transmf2)[k]<-nt
        k=k+1
}
mfdata2<-do.call(rbind, transmf2)
mfdata2$Study<-"In vivo"

mfdata<-rbind(mfdata1, mfdata2)
mfdata$Study<-factor(mfdata$Study, level=c("In vivo", "In vitro"))

#Create the averages and SE for each mut type
nuc<-c("A","T","C","G")
tb<-data.frame(nt=c(paste0(nuc, ".inVitro"),paste0(nuc,".inVivo")))

for (i in 1:nrow(tb)){
        if (i<=4){
        tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        m<-mean(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        l<-length(!is.na(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"]))
        tb$SE[i]<-sqrt(m*(1-m)/l)
        }
        if (i>4){
                
                tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
                m<-mean(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
                l<-length(!is.na(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"]))
                tb$SE[i]<-sqrt(m*(1-m)/l)
        }
}
write.csv(tb, "Output/MutFreq.filtered/Invitro.invivo.mf.summary.csv")



ggplot()+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfdata, aes(x=NT, y=mf, middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,4)]) +
        scale_fill_manual(values=paste0(colors2[c(1,4)],"66" )) +
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        guides(shape = guide_legend(override.aes = list(size = 7)))+
        geom_vline(xintercept = c(1:3)+0.5,  
                   color = "gray60", size=.4)
ggsave(filename="Output/SummaryFig/Invivo.Invitro.comparison.pdf",width=6, height=4, units='in',device='pdf')


################
#in vivo vs. in vitro by gene
dt2<-read.csv("Output/Geller/Geller_MutFrequency.csv", stringsAsFactors = F, row.names = 1)
colnames(dt2)[1]<-"pos"

genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-TS$pos[nrow(TS)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]

TS2<-merge(dt2, genetable, by="pos")

## Boxplot
transmf1<-list()
k<-1
for (i in 1:11) {
        datavector<-TS2$mean[TS2$gene==genes$Gene[(i+1)]]
        dat<-data.frame(Gene=rep(genes$Gene[(i+1)],times=length(datavector)), MF=datavector)
        transmf1[[k]]<-dat
        names(transmf1)[k]<-genes$Gene[(i+1)]
        k=k+1
}
mfGene1<-do.call(rbind, transmf1)
mfGene1$Study<-"In vitro"


k=1
transmf2<-list()
for (i in  1:11) {
        datavector<-TS$mean[TS$gene==genes$Gene[(i+1)]]
        dat<-data.frame(Gene=rep(genes$Gene[(i+1)],times=length(datavector)), MF=datavector)
        transmf2[[k]]<-dat
        names(transmf2)[k]<-genes$Gene[(i+1)]
        k=k+1
}
mfGene2<-do.call(rbind, transmf2)
mfGene2$Study<-"In vivo"

mfGene<-rbind(mfGene1, mfGene2)
mfGene$Study<-factor(mfGene$Study, level=c("In vivo", "In vitro"))
mfGene$Gene<-factor(mfGene$Gene, levels=c(genes$Gene[2:12]))

z=c(0.7,0.3,0.7,0.3)
ggplot()+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfGene, aes(x=Gene, y=MF,middle=mean(MF), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,4)]) +
        scale_fill_manual(values=paste0(colors2[c(1,4)],"66" )) +
        theme_bw()+
        theme(axis.text = element_text(size =11, color="black"), axis.title.y = element_text(size=12))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  color = "gray60", size=.4)
ggsave(filename="Output/SummaryFig/Invivo.Invitro.byGene.comparison.pdf",width=, height=4, units='in',device='pdf')
