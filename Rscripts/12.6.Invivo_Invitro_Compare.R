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
dt2<-read.csv("Output1A/Geller/Geller_MutFrequency.csv", stringsAsFactors = F, row.names = 1)
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


TS<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)

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


z=c(0.7,0.3,0.7,0.3)
ggplot(mfdata,aes(x=NT,y=mf,color=Study, fill=Study))+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(aes(middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,4)]) +
        scale_fill_manual(values=paste0(colors2[c(1,4)],"66" )) +
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.position = "none")+
        theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        geom_vline(xintercept = c(1:3)+0.5, color="gray60")
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.comparison.pdf",width=4, height=4, units='in',device='pdf')




#Add between-host mf 
among<-read.csv("Output1A/HCV1A_Combined_mutfreq.csv", row.names = 1, stringsAsFactors = F)

k=1
transmf3<-list()
for (i in  c("a","t","c","g")) {
        datavector<-among$freq.Ts[among$ref==i]
        nt<-toupper(i)
        dat<-data.frame(NT=rep(nt,times=length(datavector)), mf=datavector)
        transmf3[[k]]<-dat
        names(transmf3)[k]<-nt
        k=k+1
}
mfdata3<-do.call(rbind, transmf3)
mfdata3$Study<-"Among hosts"

mfdata<-rbind(mfdata1, mfdata2)
mfdata<-rbind(mfdata, mfdata3)
mfdata$Study<-factor(mfdata$Study, level=c( "In vitro","In vivo", "Among hosts"))


###
#Create the averages and SE for each mut type
nuc<-c("A","T","C","G")
tb<-data.frame(nt=c(paste0(nuc, ".inVitro"),paste0(nuc,".inVivo"), paste0(nuc,".AmongHosts")))

for (i in 1:nrow(tb)){
        if (i<=4){
        tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        m<-mean(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"], na.rm=T)
        l<-length(!is.na(mfdata$mf[mfdata$NT==nuc[i] & mfdata$Study=="In vitro"]))
        tb$SE[i]<-sqrt(m*(1-m)/l)
        }
        if (i>4& i<=8){
                
                tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
                m<-mean(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"], na.rm=T)
                l<-length(!is.na(mfdata$mf[mfdata$NT==nuc[i-4] & mfdata$Study=="In vivo"]))
                tb$SE[i]<-sqrt(m*(1-m)/l)
        }
        if (i>8){
                tb$Mean[i]<-mean(mfdata$mf[mfdata$NT==nuc[i-8] & mfdata$Study=="Among hosts"], na.rm=T)
                m<-mean(mfdata$mf[mfdata$NT==nuc[i-8] & mfdata$Study=="Among hosts"], na.rm=T)
                l<-length(!is.na(mfdata$mf[mfdata$NT==nuc[i-8] & mfdata$Study=="Among hosts"]))
                tb$SE[i]<-sqrt(m*(1-m)/l)
         }
}
write.csv(tb, "Output1A/MutFreq.filtered/Invitro.invivo.amonghosts.mf.summary.csv")
tb$NT<-rep(nuc, times=3)
tb$Study<-c(rep("In vitro", times=4), rep("In vivo", times=4),rep("Among hosts", times=4))
tb$Study<-factor(tb$Study, level=c("In vitro","In vivo", "Among hosts"))


ggplot()+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfdata, aes(x=NT, y=mf, middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,3,5)]) +
        scale_fill_manual(values=paste0(colors2[c(1,3,5)],"66" )) +
        geom_point(data=tb, aes(x=NT, y=Mean, group=Study), position=position_dodge(.75), size=0.6, color="black") +
        geom_errorbar(data=tb, aes(x=NT,ymin=Mean-SE, ymax=Mean+SE, group=Study), width=.2, size=.2, position=position_dodge(width=0.75))+
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        guides(shape = guide_legend(override.aes = list(size = 7)))+
        geom_vline(xintercept = c(1:3)+0.5,  
                   color = "gray60", size=.4)
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.Amonghosts.comparison.pdf",width=6, height=4, units='in',device='pdf')

## in vivo & in vitro only with mean and error bars
mfDF2<-mfdata[mfdata$Study=="In vitro"|mfdata$Study=="In vivo",]
mfDF2$Study<-factor(mfDF2$Study, levels=c("In vivo", "In vitro"))
tb2<-tb[1:8,]
tb2$Study<-factor(tb2$Study, levels=c("In vivo", "In vitro"))

ggplot()+scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mfDF2, aes(x=NT, y=mf, middle=mean(mf), color=Study, fill=Study),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(1,4)]) +
        scale_fill_manual(values=paste0(colors2[c(1,4)],"66" )) +
        geom_point(data=tb2, aes(x=NT, y=Mean, group=Study), position=position_dodge(.75), size=0.6, color="black") +
        geom_errorbar(data=tb2, aes(x=NT,ymin=Mean-SE, ymax=Mean+SE, group=Study), width=.2, size=.2, position=position_dodge(width=0.75))+
        theme_bw()+
        theme(axis.text = element_text(size =12, color="black"), axis.title.y = element_text(size=13))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        guides(shape = guide_legend(override.aes = list(size = 7)))+
        geom_vline(xintercept = c(1:3)+0.5,  
                   color = "gray60", size=.4)+
        scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.MFcomparison.pdf",width=5, height=4, units='in',device='pdf')




################
####
#in vivo vs. in vitro by gene
dt2<-read.csv("Output1A/Geller/Geller_MutFrequency.csv", stringsAsFactors = F, row.names = 1)
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

TS<-merge(TS, genetable, by="pos")
TS2<-merge(dt2, genetable, by="pos")

#TS<-TS[TS$pos>=TS2$pos[1],]
#TS2<-TS2[TS2$pos<=TS$pos[nrow(TS)],]

summary1.mean<-aggregate(TS$mean,by=list(TS$gene),FUN=mean)
#summary1.se<-aggregate(TS$mean,by=list(TS$gene),FUN=std.error)
summary2.mean<-aggregate(TS2$mean,by=list(TS2$gene),FUN=mean, na.rm=T)
#summary2.se<-aggregate(TS2$mean,by=list(TS2$gene),FUN=std.error)

#summary<-merge(summary1.mean,summary1.se, by="Group.1")
summary1.mean$Study<-"In vivo"

#summary2<-merge(summary2.mean,summary2.se, by="Group.1")
summary2.mean$Study<-"In vitro"

summary<-rbind(summary1.mean,summary2.mean)
colnames(summary)<-c("Gene", "Mean", "Study")
summary$Gene<-factor(summary$Gene, levels=c(genes$Gene[2:12]))
summary$Study<-factor(summary$Study, levels=c("In vivo", "In vitro"))

ggplot(summary,aes(x=Gene,y=Mean,color=Study, fill=Study))+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_point(data=summary, position=position_dodge(width=0.3), size=2,shape = 21)+scale_color_manual(values=colors2[c(1,4)])+
        scale_fill_manual(values=paste0(colors2[c(1,5)],"66"))+
        #geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(width=0.3))+
        theme_bw()+theme(axis.text=element_text(size=11), axis.title=element_text(size=13))+ylab("Average mutation frequency")+
        theme(panel.grid.major.x=element_blank(),axis.title.y = element_text(size=13))+
        geom_vline(xintercept = c(1:10)+0.5, color = "gray70", size=.5)+
        theme(axis.title.x=element_blank())
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.byGene.pdf",width=7, height=4)


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
        #geom_point(data=summary, aes(x=Gene, y=Mean, group=Study), position=position_dodge(.75), size=0.6, color="black") +
        #geom_errorbar(data=summary, aes(x=Gene,ymin=Mean-SE, ymax=Mean+SE, group=Study), width=.2, size=.2, position=position_dodge(width=0.75))+
        
        theme_bw()+
        theme(axis.text = element_text(size =11, color="black"), axis.title.y = element_text(size=12))+
        theme(legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
        geom_vline(xintercept = c(1:10)+0.5,  color = "gray60", size=.4)
ggsave(filename="Output1A/SummaryFig.Filtered/Invivo.Invitro.byGene.comparison.pdf",width=, height=4, units='in',device='pdf')





#### Average number of samples per site in the filtered dataset used to calculate the mean mut freqs
TS<-TS[TS$pos>=342,]
TS$Total<-apply(TS[,2:196], 1, function(x) sum(!is.na(x)) )
mean(TS$Total, na.rm=T) #141.6

6200*mean(TS$Total, na.rm=T) #877903
