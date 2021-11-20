#Create figures of features importance from random forest results

library(ggplot2)
library(gridExtra)
library(reshape2)
library(colorspace)
colors2<-qualitative_hcl(6, palette="Dark3")
cols.60<-paste0(colors2,"66")

feacol<-c("#7071FF",colors2[4])


#1. Fitness Costs
#Results from classification
feSC<-read.csv("Output/ML/FeatureImpotance_regression_all_SC.csv")
feSC<-feSC[order(feSC$Importance, decreasing = T),]
feSC$Feature<-factor(feSC$Feature, levels=paste(feSC$Feature))

ggplot(data=feSC,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=colors2[5], fill=paste0(colors2[5],"66"))+
    #geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_bw()+ylab("Importance %")+
    ggtitle("Regression")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
#ggsave("Output/ML/SC_Regression_all.pdf", width = 10, height = 3)


#top20
feSC20<-feSC[1:20,]
xlab<-c("Nonsyn","AAChange", "Nonsense","G", "Core","Nonpolar AA", "NS5B","HVR1","V-ogAA","W-ogAA","A",
        "NS1","E2","NS2","T","makesGpC","V-MutAA","NS5A", "makesCpA","A-ogAA")
p1<-ggplot(data=feSC20,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_classic()+ylab("Importance %")+
    ylim(0,60)+ggtitle("Selection coeffcients")+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=18, y=57, size=3.5, label=paste("R^2 Score = 0.623"))+
    scale_x_discrete(labels=xlab)


#Top8 of Regression
feSC8<-feSC[1:8,]
feSC8$Feature
xlabels<-c("Nonsyn","AAChange","Nonsense","G","Core","Nonpolar AA","NS5B","HVR1")

ggplot(data=feSC8,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=sccol[2], fill=paste0(sccol[1],"CC"))+
    theme_classic()+ylab("Importance %")+
    ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=7, y=57, size=3.5, label=paste("R^2 Score = 0.623"))+
    scale_x_discrete(labels=xlabels)
ggsave("Output/ML/SC_Regression_Top8.pdf",width = 3.7,height = 3.1)


### MF ###

#Results from classification
feMF<-read.csv("Output/ML/FeatureImpotance_regression_all_MF.csv")

feMF<-feMF[order(feMF$Importance, decreasing = T),]
feMF$Feature<-factor(feMF$Feature, levels=paste(feMF$Feature))

ggplot(data=feMF,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_bw()+ylab("mportance %")+
    ggtitle("Regression")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
#ggsave("Output/ML/Regression_MF_all.pdf", width = 10, height = 3)

#top20
feMF20<-feMF[1:20,]
feMF20$Feature
xlabel<-c("Nonsyn","AAChange", "Nonsense","A", "T","Core","NS5B","HVR1","V-ogAA","Hydrophobic AA",
          "NS2","C","NS1","NS5A","E2","G","W-ogAA","I-ogAA", "makesApA",  "Polar AA ")
p2<-ggplot(data=feMF20,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_classic()+ylab("Importance %")+
    ggtitle("Mutation frequency")+
    ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=18, y=57, size=3.5, label=paste("R^2 = 0.597"))
    scale_x_discrete(labels=xlabel)

#top8
feMF8<-feMF[1:8,]
feMF8$Feature
xlabel<-c("Nonsyn","AAChange", "Nonsense","A", "T","Core","NS5B","HVR1")

ggplot(data=feMF8,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_classic()+ylab("Importance %")+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=7.5, y=55, size=3.5, label=paste("R^2 = 0.597"))+
    scale_x_discrete(labels=xlabel)
ggsave("Output/ML/MF_Top8_features.pdf", width = 4, height = 3)


#Plot Top20 together

pdf("Output/ML/Both_Top20_features.pdf", width = 7, height = 7)
grid.arrange(p2, p1, nrow = 2)
dev.off()



#Overlapping top20 features
sccol<-c("#7071FF", "#3F4CE1")

fe<-feSC20[,1:2]
colnames(fe)[2]<-"SC" 
comm<-merge(fe, feMF20[,1:2], by="Feature")
#comm2<-merge(fe, feSC20[,1:2], by="Feature", all=T)
colnames(comm)[3]<-"MF"

comm<-comm[order(comm$SC, decreasing = T),]
comm$Feature<-factor(comm$Feature, levels=paste(comm$Feature))
commM<-melt(comm)
colnames(commM)[2:3]<-c("Model","Importance")
commM$Feature
labels<-c("Nonsyn","AAChange","Nonsense","G","Core","NS5B","HVR1","V-ogAA","W-ogAA","A","NS1","E2","NS2","T","NS5A")
commM$Model<-factor(commM$Model, levels = c("MF","SC"))
ggplot(data=commM,aes(x=Feature,y=Importance*100, color=Model, fill=Model))+
    geom_bar(stat="identity", position=position_dodge(width = 0.9))+
    theme_classic()+ylab("Feature importance %")+
    scale_fill_manual(values=c(paste0(colors2[4],"99"),paste0(sccol[1],"CC")) )+
    scale_color_manual(values=c(colors2[4],sccol[2]))+
    scale_x_discrete(labels=labels)+
    theme(axis.text.x = element_text(angle=90, hjust=1,size=11),axis.title.x = element_blank())+
    theme(legend.position = "none") 
ggsave("Output/ML/Feature_importance_common15.pdf", width = 6, height = 4)

comm10<-melt(comm[])
