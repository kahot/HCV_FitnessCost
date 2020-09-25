makePlot.axisbreak <- function(main){
        
        plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.00004, .1), axes = FALSE, ylab = "Mutation frequency", xlab = "Mutation type", main = main, log = "y")
        
        col.par <- "gray95"
        #    abline(h = seq(0, .15, by = .025), col = col.par)
        abline(h = (1:9)*10^(-1), col = col.par)
        abline(h = (1:9)*10^(-2), col = col.par)
        abline(h = (1:9)*10^(-3), col = col.par)
        abline(h = (1:9)*10^(-4), col = col.par)
        abline(h = (1:2)*10^(-5), col = col.par)
        abline(h = (1)*10^(-4.4), col = col.par)
        box()
        
        axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
        axis(1, at = 1:4, axislabs)
        eaxis(2, at = 10^c(-1*(0:4)))
        axis(2, at = 10^(-4.4), label = c("0"), las = 2)
        axis.break(2,2*10^-(4.5),style="slash")
        
        
}


makePlot1 <- function(main){
        plot(0, type = "n", xlim = c(.5, 2.5), ylim = c(0.0004, .1), axes = FALSE, 
                ylab = "", xlab = "", main = main, log = "y")
                mtext("Mutation frequency",2, line=2.5 )
        #    abline(h = seq(0, .15, by = .025), col = col.par)
        abline(h = (1:9)*10^(-1), col = col.par)
        abline(h = (1:9)*10^(-2), col = col.par)
        abline(h = (1:9)*10^(-3), col = col.par)
        abline(h = (1:9)*10^(-4), col = col.par)
        #abline(h = (1:2)*10^(-5), col = col.par)
        #abline(h = (1)*10^(-3.4), col = col.par)
        box()

        axislabs <- c("Syn","Nosnyn")
        axis(1, at = 1:2, axislabs, cex=.5)
        eaxis(2, at = 10^c(-1*(0:4)))
        axis(2, at = 10^(-3.4), label = c("0"), las = 2)
        axis.break(2, 2*10^-(3.65),style="slash")
}

makePlot2 <- function(main){
        plot(0, type = "n", xlim = c(.5, 2.5), ylim = c(0.0004, .1), axes = FALSE, 
             ylab = "", xlab = "", main = main, log = "y")
        
        #abline(h = seq(0, .15, by = .025), col = col.par)
        abline(h = (1:9)*10^(-1), col = col.par)
        abline(h = (1:9)*10^(-2), col = col.par)
        abline(h = (1:9)*10^(-3), col = col.par)
        abline(h = (1:9)*10^(-4), col = col.par)
        #abline(h = (1:2)*10^(-5), col = col.par)
        #abline(h = (1)*10^(-3.4), col = col.par)
        box()
        
        axislabs <- c("Syn","Nosnyn")
        axis(1, at = 1:2, axislabs)
        eaxis(2, at = 10^c(-1*(0:4)))
        axis(2, at = 10^(-3.4), label = c("0"), las = 2)
        axis.break(2, 2*10^-(3.65),style="slash")
        
}


makeDataFrame1 <- function(NsOrNo , CpGorNo, Core, HVR){
        
        atcg.mat <- as.data.frame(matrix(data = 0, ncol = nrow(modcoef), nrow = 4))
        #ATCG elements
        diag(atcg.mat[1:4, 1:4]) <- 1
        #Reserve the first column for intercept
        atcg.mat[,1] <- 1
        #synonymous or nonsynonymous mutation?
        atcg.mat[,6] <- NsOrNo
        #Core gene?
        atcg.mat[,8] <- Core 
        #Core gene?
        atcg.mat[,10] <- HVR
        #CpG mutation or not
        atcg.mat[1:2,5] <- CpGorNo
        atcg.mat[3:4,5] <- 0
        
        #nonysynonymous interactions with a t c g
        atcg.mat[16:18] <- atcg.mat[,2:4] * atcg.mat[,6]
        
        #names
        names(atcg.mat) <- c(rownames(modcoef) )
        
        return(as.matrix(atcg.mat))
} 


plotPoint <- function(NT, NsOrNo, CpGorNo, Core, HVR,colVal, offset){
        
        pchpar <- 21
        cexpar <- 1
        arrowlwd <- 1.5
        arrowlen <-.08
        
        setUpDat <- makeDataFrame1(NsOrNo, CpGorNo, Core, HVR)[,rownames(modcoef)]
        
        if (NT=="A") setUpDat<-setUpDat[1,]
        if (NT=="T") setUpDat<-setUpDat[2,]
        if (NT=="C") setUpDat<-setUpDat[3,]
        if (NT=="G") setUpDat<-setUpDat[4,]
        
        #compute
        #%*% = matrix manipulation
        pointToPlot <- exp(setUpDat %*% coef.vals)
        arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
        arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)
        
        #plot
        if (NsOrNo == 0) xpos<-1
        if (NsOrNo == 1) xpos<-2
        arrows(xpos + offset, arrowToPlotLow, xpos + offset, arrowToPlotHigh, length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = colVal, bg="black")
        points(xpos + offset, pointToPlot, col =colVal, bg = "black",  pch = pchpar, cex = .8)
        
}


plotPoint2 <- function(NT, NsOrNo, CpGorNo, Core, HVR,colVal, offset){
        
        pchpar <- 21
        cexpar <- 1
        arrowlwd <- 1.5
        arrowlen <-.08
        
        setUpDat <- makeDataFrame1(NsOrNo, CpGorNo, Core, HVR)[,rownames(modcoef)]
        
        if (NT=="A") setUpDat<-setUpDat[1,]
        if (NT=="T") setUpDat<-setUpDat[2,]
        if (NT=="C") setUpDat<-setUpDat[3,]
        if (NT=="G") setUpDat<-setUpDat[4,]
        
        #compute
        #%*% = matrix manipulation
        pointToPlot <- exp(setUpDat %*% coef.vals)
        arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
        arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)
        
        #don't plot estimates for CpG sites at Cs or Gs
        restrict <- 1:4
        if(CpGorNo == 1){
                restrict <- 1:2
        }
        
        #plot
        arrows(2 + offset, arrowToPlotLow, 2 + offset, arrowToPlotHigh, length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = colVal)
        points(2 + offset, pointToPlot, col =colVal, bg = substr(colVal,start=1, stop=7),  pch = pchpar, cex = .8)
        #arrows((1:4)[restrict] + offset, arrowToPlotLow[restrict], (1:4)[restrict] + offset, arrowToPlotHigh[restrict], length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = "black")
        
}
nucord <- c("a", "t", "c", "g")

plotDots <- function(NT, NsOrNo, CpGorNo, Core, HVR,colVal, offsetval){
        
        pchpar <- 16
        cexpar <- .6
        if (Core==1) genename<-"Core"
        if (HVR==1) genename<-"HVR1"
        nsinds <- which(d$Type == c("syn", "nonsyn")[NsOrNo + 1])
        cpginds <- which(d$makesCpG == CpGorNo)
        geneinds <- which(d$gene == genename)
        if (NT=="A") i<-1
        if (NT=="T") i<-2
        if (NT=="C") i<-3
        if (NT=="G") i<-4
        ntinds<-which(d$ref==nucord[i])
        
        datInds <- intersect(intersect(intersect(nsinds, cpginds), geneinds),ntinds)
        if (NsOrNo==0) xinds <- rep(1, times=length(datInds))
        if (NsOrNo==1) xinds <- rep(2, times=length(datInds))
        #xinds <- as.numeric(xinds)
        yinds <- d$mean[datInds]
        xjit <- rnorm(length(xinds), offsetval, .03) #jitter
        
        yinds[yinds == 0] <- 10^(-4.4)
        
        points(xinds + xjit, yinds, col = colVal,  pch = pchpar, cex = cexpar)
}
