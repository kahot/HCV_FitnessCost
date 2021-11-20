# Plot the beta regression results
library(colorspace)
library(reshape2)
library(ggplot2)
library(sfsmisc)
source("Rscripts/BetaRegPlotFunction.R")
colors2<-qualitative_hcl(6, palette="Dark3")

# Read the mut freq. data
d<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",row.names = 1)
d<-d[d$pos>=342,]

# Beta regression results file (model coefficients)
modcoef<-read.csv("Output/BetaReg/BetaReg_MF_effectsize.csv",stringsAsFactors = F)

coef.vals <- modcoef[,1]
coef.pSE.vals <- coef.vals + modcoef[,2]
coef.mSE.vals <- coef.vals - modcoef[,2]

# plot syn vs nonsyn, estimated and observed MF

#png("Output/BetaReg/Betareg.plot.png", width = 4,height = 4, res=300, units = "in")
pdf("Output/BetaReg/Betareg.plot.pdf", width = 4.2,height = 4.3)
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

