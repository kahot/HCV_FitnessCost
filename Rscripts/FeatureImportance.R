#Create figures of features importance from random forest results

library(ggplot2)
library(gridExtra)
library(reshape2)
library(colorspace)
colors2<-qualitative_hcl(6, palette="Dark3")
cols.60<-paste0(colors2,"66")

feacol<-c(colors2[5],"#d95f02")
#          "#ff7f00")

    
#1. Fitness Costs
#Results from classification
feReg<-read.csv("Output/ML/FeatureImpotance_regression_all_SC.csv")
feCla<-read.csv("Output/ML/FeatureImpotance_classification_all_SC.csv")

feReg<-feReg[order(feReg$Importance, decreasing = T),]
feReg$Feature<-factor(feReg$Feature, levels=paste(feReg$Feature))
feCla<-feCla[order(feCla$Importance, decreasing = T),]
feCla$Feature<-factor(feCla$Feature, levels=paste(feCla$Feature))

ggplot(data=feReg,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=colors2[5], fill=paste0(colors2[5],"66"))+
    #geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_bw()+ylab("Importance %")+
    ggtitle("Regression")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
ggsave("Output/ML/SC_Regression_all.pdf", width = 10, height = 3)

ggplot(data=feCla,aes(x=Feature,y=Importance*100))+
    #geom_bar(stat="identity", color="deeppink", fill="pink")+
    geom_bar(stat="identity", color=feacol[2], fill=paste0(feacol[2],"66"))+
    theme_bw()+ylab("Importance %")+
    ggtitle("Classification")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
ggsave("Output/ML/SC_Classification_all.pdf", width = 10, height = 3)

feReg$Method<-"Regression"
feCla$Method<-"Classification"
fea<-rbind(feReg, feCla)

ggplot(data=fea,aes(x=Feature,y=Importance*100, fill=Method))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))+
    theme_bw()+ylab("Importance %")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())

ggplot(data=fea,aes(x=Feature,y=Importance*100, fill=Method, color=Method))+
    facet_wrap(~Method, ncol=1)+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))+
    theme_bw()+ylab("Importance %")+
    scale_fill_manual(values=paste0(feacol,"66"))+
    scale_color_manual(values=feacol)+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())


#top20
feReg20<-feReg[1:20,]
feCla20<-feCla[1:20,]

ggplot(data=feReg20,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_classic()+ylab("Importance %")+
    ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=18, y=57, size=3.5, label=paste("R^2 Score = 0.623"))
ggsave("Output/ML/SC_Regression_Top20.pdf",width = 6,height = 3)


p1<-ggplot(data=feReg20,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_classic()+ylab("Importance %")+
    ggtitle("Regression")+ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=18, y=57, size=3.5, label=paste("R^2 Score = 0.623"))

p2<-ggplot(data=feCla20,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[2], fill=paste0(feacol[2],"66"))+
    theme_classic()+ylab("Importance %")+
    ggtitle("Classification")+ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=18, y=57, size=3.5, label=paste("Accuracy = 0.611"))
pdf("Output/ML/SC_Both_Top20_features.pdf", width = 7, height = 6)
grid.arrange(p1, p2, nrow = 2)
dev.off()


fe20<-rbind(feReg20,feCla20)
ggplot(data=fe20,aes(x=Feature,y=Importance*100, fill=Method))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))+
    theme_bw()+ylab("Importance %")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())

fe20$Method<-factor(fe20$Method, levels=c("Regression","Classification"))
ggplot(data=fe20,aes(x=Feature,y=Importance*100, fill=Method, color=Method))+
    facet_wrap(~Method, ncol=1)+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))+
    theme_bw()+ylab("Importance %")+
    scale_fill_manual(values=paste0(feacol,"66"))+
    scale_color_manual(values=feacol)+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank(),legend.position = "none") 
ggsave("Output/ML/SC_Both_Top20_facet_plot.pdf", width = 7, height = 6)

#Top8 of Regression
feReg8<-feReg[1:8,]
xlabels<-c("Nonsyn","AAChange","Nonsense","G","Core","Nonpolar AA","NS5B","HVR1")

ggplot(data=feReg8,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_classic()+ylab("Importance %")+
    ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=7, y=57, size=3.5, label=paste("R^2 Score = 0.623"))+
    scale_x_discrete(labels=xlabels)
ggsave("Output/ML/SC_Regression_Top8.pdf",width = 4,height = 3)





#Overlapping top20 features

fe<-feReg20[,1:2]
colnames(fe)[2]<-"Regression" 
comm<-merge(fe, feCla20[,1:2], by="Feature")
#comm2<-merge(fe, feReg20[,1:2], by="Feature", all=T)
colnames(comm)[3]<-"Classification"
commM<-melt(comm)
colnames(commM)[2:3]<-c("Method","Importance")
labels<-c("Nonsyn","G","AAChange","Nonsense","A","Core","NS5B","V-ogAA","CpA","A-ogAA")
ggplot(data=commM,aes(x=Feature,y=Importance*100, color=Method, fill=Method))+
    geom_bar(stat="identity", position=position_dodge(width = 0.9))+
    theme_classic()+ylab("Feature importance %")+
    scale_fill_manual(values=c(paste0(feacol[1],"66"),paste0(feacol[2],"66")) )+
    scale_color_manual(values=feacol)+
    scale_x_discrete(labels=labels)+
    theme(axis.text.x = element_text(angle=90, hjust=1,size=11),axis.title.x = element_blank())
ggsave("Output/ML/Feature_importance_common10.pdf", width = 6, height = 4)

#ggplot(data=commM,aes(x=Feature,y=Importance*100, color=Method, fill=Method))+
#    geom_bar(stat="identity", position=position_dodge(width = 0.8))+
#    theme_bw()+ylab("Feature importance %")+
#    scale_fill_manual(values=c("pink","lightblue"))+
#    scale_color_manual(values=c("deeppink","blue"))+
#    scale_x_discrete(labels=labels)+
#    theme(axis.text.x = element_text(angle=90, hjust=1,size=11),axis.title.x = element_blank())
#ggsave("Output/ML/Feature_importance_common11_v2.pdf", width = 6, height = 4)
#


#Selected features

feSReg<-read.csv("Output/ML/FeatureImpotance_regression_selected_SC.csv")
feSReg<-feSReg[order(feSReg$Importance, decreasing = T),]
feSReg$Feature<-factor(feSReg$Feature, levels=paste(feSReg$Feature))
feSCla<-read.csv("Output/ML/FeatureImpotance_classification_selected_SC.csv")
feSCla<-  feSCla[order(feSCla$Importance, decreasing = T),]
feSCla$Feature<-factor(feSCla$Feature, levels=paste(feSCla$Feature))

p3<-ggplot(data=feSReg,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[1], fill=paste0(feacol[1],"66"))+
    theme_classic()+ylab("Importance %")+
    ggtitle("Regression")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())+
    annotate('text', x=6.8, y=75, size=3, label=paste("R^2 Score = 0.583"))
p4<-ggplot(data=feSCla,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color=feacol[2], fill=paste0(feacol[2],"66"))+
    theme_classic()+ylab("Importance %")+
    ggtitle("Classification")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())+
    annotate('text', x=9.8, y=42, size=3, label=paste("Accuracy = 0.608"))

pdf("Output/ML/SC_Selected_features.pdf", width = 4.5, height = 6)
grid.arrange(p3, p4, nrow = 2)
dev.off()
    

### MF ###

#Results from classification
feReg<-read.csv("Output/ML/FeatureImpotance_regression_all_MF.csv")

feReg<-feReg[order(feReg$Importance, decreasing = T),]
feReg$Feature<-factor(feReg$Feature, levels=paste(feReg$Feature))

ggplot(data=feReg,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_bw()+ylab("Importance %")+
    ggtitle("Regression")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
ggsave("Output/ML/Regression_MF_all.pdf", width = 10, height = 3)


#top20
feReg20<-feReg[1:20,]
xlabel<-c("Nonsyn","AAChange", "Nonsense-MutAA","A", "T","Core","NS5B","HVR1","V-ogAA","Hydrophobic AA",
          "NS2","C","NS1","NS5A","E2","G","W-ogAA","I-ogAA", "makesApA",  "Polar AA ")

ggplot(data=feReg20,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_classic()+ylab("Importance %")+
    #ylim(0,60)+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=18, y=55, size=3.5, label=paste("R^2 = 0.597"))+
    scale_x_discrete(labels=xlabel)
ggsave("Output/ML/MF_Top20_features.pdf", width = 7, height = 3)

#top8
feReg8<-feReg[1:8,]
xlabel<-c("Nonsyn","AAChange", "Nonsense","A", "T","Core","NS5B","HVR1")

ggplot(data=feReg8,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_classic()+ylab("Importance %")+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank(), panel.grid.major.x = element_blank())+
    annotate('text', x=7.5, y=55, size=3.5, label=paste("R^2 = 0.597"))+
    scale_x_discrete(labels=xlabel)
ggsave("Output/ML/MF_Top8_features.pdf", width = 4, height = 3)

feReg20$Feature
#xlabel<-c("Nonsyn","AAChange", "Nonsense-MutAA","A", "T","Core","NS5B","HVR1","V-ogAA","Hydrophobic AA",
          "NS2","C")

#Plot with beta regression results?
### Plot the effect size ###
effects<-read.csv("Output/BetaReg/BetaReg.effects.csv",stringsAsFactors = F)
effects<-effects[order(abs(effects$percent), decreasing = T),]

effects$factor<-factor(effects$factor, levels=effects$factor[1:17])

effects$percent<-effects$percent*100

#Horizontal
ggplot(effects, aes(factor,percent)) +
    geom_bar(stat="identity", color=colors2[4], fill=paste0(colors2[4],"CC"))+
    theme_classic() +
    #theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    #theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")
ggsave("Output/ML/BetaReg_mf_effects_ordered1.pdf", width =5.5 ,height = 3.7)

#order2
effects<-effects[order(abs(effects$effect), decreasing = T),]
effects$factor<-factor(effects$factor, levels=effects$factor[1:17])
ggplot(effects, aes(factor,percent)) +
    geom_bar(stat="identity", color=colors2[4], fill=paste0(colors2[4],"CC"))+
    theme_classic() +
    #theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    #theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")
ggsave("Output/ML/BetaReg_mf_effects_ordered1.pdf", width =5.5 ,height = 3.7)





#Selected features
feSReg<-read.csv("Output/ML/FeatureImpotance_regression_selected_MF.csv")
feSReg<-feSReg[order(feSReg$Importance, decreasing = T),]
feSReg$Feature<-factor(feSReg$Feature, levels=paste(feSReg$Feature))

ggplot(data=feSReg,aes(x=Feature,y=Importance*100))+
    geom_bar(stat="identity", color="darkgreen", fill=cols.60[4])+
    theme_classic()+ylab("Importance %")+
    theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())+
    annotate('text', x=7.5, y=67, size=3, label=paste("R^2 Score = 0.543"))

ggsave("Output/ML/MF_Selected_features.pdf", width = 4, height = 3)




data<-read.csv("~/programs/HCV_ML/Data/allfeatturesv8_dinucs_nucs_MutFreqBelow50%.csv")
hist(data$Avg_Mutation_Freq)
hist(log(data$Avg_Mutation_Freq))
hist(log(data$Avg_Mutation_Freq,10))
hist(log(data$Avg_Mutation_Freq,2))

hist(data$Costs)
hist(log(data$Costs))
hist(log10(data$Costs))
hist(log2(data$Costs))

range(data$Costs)
#0.00000637 0.59788360

ggplot(data, aes(x=pos, y=Costs))+
        geom_bar(stat = 'identity')+
        scale_y_continuous(trans="log")
ggplot(data, aes(x=pos, y=Costs))+
        geom_point(color="blue", alpha=0.5, size=.2)+
        scale_y_continuous(trans="log")
nrow(data[data$Costs>=0.01,]) #51


nrow(data[data$Avg_Mutation_Freq==0,]) #0
table(data$Cost.Label)
# High       Low Very High  Very low 
#2495      1937      2291      1433 

aggregate(data["Avg_Mutation_Freq"], by=list(data$Cost.Label), mean)
#    Group.1 Avg_Mutation_Freq
#1      High       0.003372113
#2       Low       0.039830187
#3 Very High       0.002844994
#4  Very low       0.181862666
aggregate(data["Avg_Mutation_Freq"], by=list(data$Cost.Label), range)
#    Group.1 Avg_Mutation_Freq.1 Avg_Mutation_Freq.2
#1      High         0.001280108         0.014400394
#2       Low         0.003478993         0.175304132
#3 Very High         0.000009380         0.005324279
#4  Very low         0.042419238         0.503618179

aggregate(data["Costs"], by=list(data$Cost.Label), mean)
aggregate(data["Costs"], by=list(data$Cost.Label), range)


features2<-read.csv("~/programs/HCV_ML/FeatureImpotance_all_regression_fq.csv")


features2<-features2[order(features2$Importance, decreasing = T),]
features2$Feature<-factor(features2$Feature, levels=paste(features2$Feature))

ggplot(data=features2,aes(x=Feature,y=Importance*100))+
        geom_bar(stat="identity", color="red", fill="pink")+
        theme_bw()+ylab("Importance %")+
        theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
ggsave("~/programs/HCV_ML/Figures/Feature_importance_mf_regression.pdf", width = 10, height = 3)


fe2<-features2[1:20,]
ggplot(data=fe2,aes(x=Feature,y=Importance*100))+
        geom_bar(stat="identity", color="violetred", fill="lightpink")+
        theme_bw()+ylab("Importance %")+ggtitle("Mutation frequencies")+
        annotate('text', x=16, y=50, size=3, label=expression(paste("Regression: ",r^{2},"=0.59")))+
        theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())
ggsave("~/programs/HCV_ML/Figures/Feature_importance_mf_regression_top20.pdf", width = 4.8, height = 3)

#find common elements from classification and regression
intersect(fe20$Feature, fe2$Feature)
#[1] "Nonsyn"      "g-refN"      "bigAAChange" "c-refN"      "*-MutAA"    
#[6] "a-refN"      "NS5B"        "Core"        "t-refN"      "5' UTR"     
#[11] "E1"   

#11 are overlapping


features3<-read.csv("~/programs/HCV_ML/FeatureImpotance_all_regression_costs.csv")

features3<-features3[order(features3$Importance, decreasing = T),]
features3$Feature<-factor(features3$Feature, levels=paste(features3$Feature))

ggplot(data=features3,aes(x=Feature,y=Importance*100))+
        geom_bar(stat="identity", color="red", fill="pink")+
        theme_bw()+ylab("Importance %")+ggtitle("Fitness costs")+
        theme(axis.text.x = element_text(angle=90, hjust=1),axis.title.x = element_blank())
ggsave("~/programs/HCV_ML/Figures/Feature_importance_costs_regression.pdf", width = 10, height = 3)


fe3<-features3[1:20,]
ggplot(data=fe3,aes(x=Feature,y=Importance*100))+
        geom_bar(stat="identity", color="violetred", fill="lightpink")+
        theme_bw()+ylab("Importance %")+ggtitle("Fitness costs")+
        annotate('text', x=16, y=56, size=3, label=expression(paste("Regression: ",r^{2},"=0.62")))+
        theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())
ggsave("~/programs/HCV_ML/Figures/Feature_importance_costs_regression_top20.pdf", width = 5, height = 3)




#divide the mutation frequencies into 10 categories:

aggregate(data["Avg_Mutation_Freq"], by=list(data$Cost.Label), range)
#    Group.1 Avg_Mutation_Freq.1 Avg_Mutation_Freq.2
#1      High         0.001280108         0.014400394
#2       Low         0.003478993         0.175304132
#3 Very High         0.000009380         0.005324279
#4  Very low         0.042419238         0.503618179

intersect(fe20$Feature, fe3$Feature)
# [1] "Nonsyn"      "g-refN"      "bigAAChange" "c-refN"      "*-MutAA"    
#[6] "a-refN"      "NS5B"        "Core"        "t-refN"      "5' UTR"     

intersect(fe2$Feature, fe3$Feature)
#[1] "Nonsyn"         "bigAAChange"    "*-MutAA"        "a-refN"        
#[5] "5' UTR"         "t-refN"         "Core"           "NS5B"          
#[9] "Polar AA"       "g-refN"         "HVR1"           "c-refN"        
#[13] "Hydrophobic AA" "NS5A"           "V-ogAA"         "W-ogAA"        
#[17] "NS1"#