# run BetaReg first (11.1 & 11.2)
library(plotrix)
library(sfsmisc)
library(colorspace)
library(reshape2)
source("Rscripts/baseRscript.R")

colors2<-qualitative_hcl(6, palette="Dark3")
cols4.60<-paste0(cols4,"66")

####  Make Plots  #####

brData<-read.csv("Output/BetaReg/BetaReg.Data.Shape.csv", stringsAsFactors = F, row.names = 1)
d<-read.csv("Output/MutFreq.filtered/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,]

#Results file: modcoeff data (model.g2)
modcoef<-read.csv("Output/BetaReg/BetaReg_BestModel.csv",stringsAsFactors = F)
rownames(modcoef)<-modcoef$X
modcoef<-modcoef[,-1]
modcoef<-modcoef[c("(Intercept)","t","c","g","CpG","Nonsyn","bigAAChange","t:Nonsyn","c:Nonsyn","g:Nonsyn"),]



####
# plot syn vs nonsyn estiamted and observed mf
source("Rscripts/BetaRegPlotFunction2.R")
pdf("Output/GLM/Betareg.plot.pdf", width = 5.5,height = 5)
col.par <- "gray95"
makePlotAll(main="")
plotMFs(0, paste0(colors2[5],"66"),  -.1)
plotMFs(1, paste0(colors2[1],"33"),  .1)

plotEsp(0,1, -.1)
plotEsp(1,1, .1)

legend("topright", c("syn", "nonsyn"), col = colors2[c(5,1)], pch = 16, bg = "white", cex=.8, bty = "white" )
dev.off()


######
### Plot the effect size ###
effects<-read.csv("Output/BetaReg/BetaReg.effects.csv",stringsAsFactors = F)
effects$factor<-factor(effects$factor, levels=effects$factor[1:17])

effects$percent<-effects$percent*100
effects$percent_SC<-effects$percent_SC*100

#Horizontal
ggplot(effects, aes(factor,percent)) +
        geom_bar(stat="identity", color=colors2[4], fill=paste0(colors2[4],"CC"))+
        theme_test() +
        theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
        theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
        labs(x="", y="Estimated effects (%)")

ggsave("Output/SummaryFig.Filtered/BetaReg.effects.pdf", width = 9, height = 6)

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

