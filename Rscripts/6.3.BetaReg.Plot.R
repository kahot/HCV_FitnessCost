# run BetaReg first (11.1 & 11.2)
#library(plotrix)
#library(sfsmisc)
library(colorspace)
library(reshape2)
library(ggplot2)
#source("Rscripts/baseRscript.R")
source("Rscripts/BetaRegPlotFunction.R")

colors2<-qualitative_hcl(6, palette="Dark3")

# Read the mut freq. data
d<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,]

# Beta regression results file (model coefficients)
modcoef<-read.csv("Output/BetaReg/BetaReg_MF_BestModel.csv",stringsAsFactors = F)
rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
modcoef<-modcoef[c("(Intercept)","t","c","g","CpG","Nonsyn","bigAAChange","t:Nonsyn","c:Nonsyn","g:Nonsyn"),]

coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

# plot syn vs nonsyn, estimated and observed MF

pdf("Output/BetaReg/Betareg.plot.pdf", width = 4.8,height = 5)
col.par <- "gray95"
makePlotAll(main="")
plotMFs(0, paste0(colors2[5],"66"),  -.1)
plotMFs(1, paste0(colors2[1],"33"),  .1)

plotEsp(0,1, -.1)
plotEsp(1,1, .1)

legend("topright", c("syn", "nonsyn"), col = colors2[c(5,1)], pch = 16, bg = "white", cex=.8, bty = "white" )
dev.off()



### Plot the effect size ###
effectsMF<-read.csv("Output/BetaReg/BetaReg_MF_effectsize.csv",stringsAsFactors = F)
effectsMF<-effectsMF[effectsMF$factor!="(Intercept)",]
effectsMF$factor<-factor(effectsMF$factor, levels=effectsMF$factor)
effectsMF$percent<-effectsMF$Effect*100

ggplot(effectsMF, aes(factor,percent)) +
        geom_hline(yintercept = 0, size=0.5, color="gray40")+
        geom_bar(stat="identity", width=0.8,size=0.5,color=colors2[4], fill=paste0(colors2[4],"66"))+
        theme_test() +theme_classic()+
        theme(axis.text=element_text(size=10.5, color="gray20"), axis.title.y=element_text(size=12))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
        labs(x="", y="Estimated effects (%)")
ggsave("Output/BetaReg/BetaReg.effectSizes.pdf", width = 6.7, height = 3.4)




effects<-read.csv("Output/BetaReg/BetaReg.effects.csv",stringsAsFactors = F)
effects$factor<-factor(effects$factor, levels=effects$factor[1:17])

effects$percent<-effects$percent*100
effects$percent_SC<-effects$percent_SC*100

#Horizontal
ggplot(effectsMF, aes(factor,percent)) +
        geom_hline(yintercept = 0, size=0.5, color="gray40")+
        geom_bar(stat="identity", width=0.8,size=0.8,color=colors2[4], fill=paste0(colors2[4],"B3"))+
        theme_test() +theme_classic()+
        theme(axis.text=element_text(size=11), axis.title.y=element_text(size=12))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
         labs(x="", y="Estimated effects (%)")

ggsave("Output/BetaReg/BetaReg.effectSizes.pdf", width = 6.7, height = 3)

##effect size for SC
ggplot(effects, aes(factor,percent_SC)) +
        geom_bar(stat="identity", color=colors2[5], fill=paste0(colors2[5],"CC"))+
        theme_test() +
        theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
        theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
        labs(x="", y="Estimated effects (%)")

ggsave("Output/SummaryFig.Filtered/BetaReg.effects.SC.pdf", width = 9, height = 6)


##effect size for MF & SC
colnames(effects)[4:5]<-c("MF","SC")
effects2<-melt(effects[,c(2,4,5)], id.var=c("factor"))
colnames(effects2)[2:3]<-c("Type","percent")

ggplot(effects2, aes(factor,percent, group=Type, fill=Type)) +
        scale_fill_manual(values=c(paste0(colors2[4],"CC"),paste0(colors2[1],"CC")), labels=c("Mut freq","Sel coeff"))+
        scale_color_manual(values=c(colors2[4], colors2[1]))+
        geom_bar(stat="identity",color="gray",position=position_dodge(width=.5))+
        theme_test() +
        theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
        theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
        labs(x="", y="Estimated effects (%)")+
        theme(legend.title = element_blank())

ggsave("Output/SummaryFig.Filtered/BetaReg.effects.MF-SC2.pdf", width = 9, height = 5)

