# Create the ranking ordered figure

library(dplyr)
library(colorspace)
source("Rscripts/baseRscript.R")
colors2<-qualitative_hcl(6, palette="Dark3")

plotter <- function(d){
        par(mar = c(4.5, 4.5, 2, 2))
        main.dat <- d
        wheresthebreak <- 4
        remap = 1       
        dataset <- main.dat[complete.cases(main.dat), ]
        toPlot <- dataset[,"Freq"]
        toPlot <- toPlot[!is.na(toPlot)]
        colVect <- rep(0, nrow(dataset))
        colVect[dataset$TypeOfSite == "nonsyn"] <- colors2[1]
        colVect[dataset$TypeOfSite == "syn"] <- colors2[5]
        colVect[dataset$TypeOfSite == "stop"] <- "black"
        plot(5, type = "n", log = "y", axes = FALSE, xlim = c(0, length(toPlot[!is.na(toPlot)])), 
             ylim = c(0.0001, 0.1),  
             ylab = "Average mutation frequency", xlab = "Mutations ordered by average mutation frequency",
             cex.lab = 1.3)
        for(i in 1:4){
                abline(h = 1:10 * 10^(-i), col = "gray70")
        }
        abline(h = 10^-4, col = "gray70")
        eaxis(side = 2, at = 10^((-1):(-5)))
        box()
        
        cexval <- 1.5
        toPlot[toPlot == 0] <- 10^-(wheresthebreak)
        points(1:length(toPlot), sort(toPlot), col = colVect[order(toPlot)], pch = "|", cex = cexval)
        axis(1)
        legend("bottomright", c("Synonymous", "Non-synonymous", "Nonsense"), col = c(colors2[5], colors2[1], "black"), pch = "|", bg = "white", pt.cex = cexval)
}


TS<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)
#coding region only
TS<-TS[TS$pos>341,196:ncol(TS)]
colnames(TS)[colnames(TS)=="Type"]<-"TypeOfSite"
colnames(TS)[colnames(TS)=="mean"]<-"Freq"

pdf("Output/SummaryFig/MutFreq-Ordered.pdf", height = 5.5, width = 9)
plotter(TS)
dev.off() 
        
