library(ggplot2)
library(colorspace)
source("Rscripts/label_scientific.R")

colors2<-qualitative_hcl(6, palette="Dark3")
source("Rscripts/baseRscript.R")

mfs<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv",stringsAsFactors = F, row.names = 1)
mfs<-mfs[mfs$pos>341,]
### Plot mutation freq. across the genome based on the mutation types 

n<-data.frame("pos"=c(1:mfs$pos[nrow(mfs)]))
mfs<-merge(n,mfs,by="pos",all.x=T)
mfs2<-mfs[,c(1,197:203)]
hmlm<-read.csv("Data/HMLMsites.csv")

HL<-list()
for (type in c("LM","HM")){
        hl<-hmlm[hmlm$Type==type,]
        MF<-list()
        SE<-list()
        for (i in 1:nrow(hl)){
                df<-mfs2[mfs2$pos>=hl$Start[i]&mfs2$pos<=hl$End[i],]
                MF[[i]]<-mean(df$mean, na.rm = T)
                m<-mean(df$mean, na.rm=T)
                l<-length(!is.na(df$mean))
                SE[[i]]<-sqrt(m*(1-m)/l)
                names(MF)[i]<-paste0(hl$Type[i],i,"(",hl$Gene[i],")")
                names(SE)[i]<-paste0(hl$Type[i],i,"(",hl$Gene[i],")")
        }
        MFave<-as.data.frame(do.call(rbind,MF))
        colnames(MFave)[1]<-"MF"
        MFave$ID<-rownames(MFave)
        se<-as.data.frame(do.call(rbind,SE))
        MFave$SE<-se$V1
        MFave$Type<-type
        HL[[type]]<-MFave
        
}

MFave<-do.call(rbind, HL)
MFave$ID<-factor(MFave$ID, levels=paste0(MFave$ID))

#######     
#Box plotall mutation freq ()

mf<-mfs2[,c("pos", "mean")]
HMLM<-list()

for (type in c("LM","HM")){
        hl<-hmlm[hmlm$Type==type,]
        MF<-list()
        for (i in 1:nrow(hl)){
                df<-mf[mf$pos>=hl$Start[i]&mf$pos<=hl$End[i],]
                id<-paste0(hl$Type[i],i,"(",hl$Gene[i],")")
                df$ID<-id
                MF[[i]]<-df
                names(MF)[i]<-id
        }
        mufreq<-as.data.frame(do.call(rbind,MF))
        mufreq$Type<-type
        colnames(mufreq)[2]<-"MF"
        HMLM[[type]]<-mufreq
        
}

mufreq<-do.call(rbind, HMLM)
mufreq$ID<-factor(mufreq$ID, levels=paste0(MFave$ID))

#allmean<-mean(mfs$mean,na.rm = T)

ggplot()+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
        geom_boxplot(data=mufreq, aes(x=ID, y=MF, color=Type, fill=Type),outlier.alpha = 0.2)+
        labs(x="",y="Transition mutation frequency")+
        scale_color_manual(values=colors2[c(4,6)])+
        scale_fill_manual(values=paste0(colors2[c(4,6)],"66"))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=90, hjust=0, color=1))+
        geom_hline(yintercept=allmean, color="gray60", linetype=2)+
        geom_text(label="Mean", x= 20.5, y= log10(allmean), size=3, color="gray20")+
        theme(plot.margin = unit(c(.5, 1, .5, .5), "cm"))+
        coord_cartesian(xlim = c(1, 19), clip = 'off') +
        theme(legend.justification = "top")

ggsave("Output1A/SummaryFig/HM.LM.sites.boxplots.pdf", width = 6.5, height = 4)
