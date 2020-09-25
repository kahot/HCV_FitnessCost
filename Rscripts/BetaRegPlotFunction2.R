nucord <- c("a", "t", "c", "g")


makePlot.axisbreak <- function(main){
        
        plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.00004, .1), axes = FALSE, ylab = "Mutation frequency", xlab = "Mutation type", main = main, log = "y")
        
        col.par <- "gray95"
        #    abline(h = seq(0, .15, by = .025), col = col.par)
        abline(h = (1:9)*10^(-1), col = col.par)
        abline(h = (1:9)*10^(-2), col = col.par)
        abline(h = (1:9)*10^(-3), col = col.par)
        abline(h = (1:9)*10^(-4), col = col.par)
        #abline(h = (1:2)*10^(-5), col = col.par)
        #abline(h = (1)*10^(-4.4), col = col.par)
        box()
        
        axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
        axis(1, at = 1:4, axislabs)
        eaxis(2, at = 10^c(-1*(0:4)))
        axis(2, at = 10^(-4.4), label = c("0"), las = 2)
        axis.break(2,2*10^-(4.5),style="slash")
        
        
}


makePlotAll <- function(main){
        plot(0, type = "n", xlim = c(.5, 4.5), ylim = c(0.0004, .1), axes = FALSE, 
             ylab = "", xlab = "", main = main, log = "y")
        mtext("Mutation frequency",2, line=2.5 )
        #    abline(h = seq(0, .15, by = .025), col = col.par)
        abline(h = (1:9)*10^(-1), col = col.par)
        abline(h = (1:9)*10^(-2), col = col.par)
        abline(h = (1:9)*10^(-3), col = col.par)
        abline(h = (1:9)*10^(-4), col = col.par)
        abline(v=c(1:3)+0.5, col="gray70", cex=.7)
        box()
        
        axislabs <- c(expression("A" %->% "G"), expression("T" %->% "C"), expression("C" %->% "T"), expression("G" %->% "A"))
        axis(1, at = 1:4, axislabs, cex=.5)
        eaxis(2, at = 10^c(-1*(0:4)))
        axis(2, at = 10^(-3.4), label = c("0"), las = 2)
        axis.break(2, 2*10^-(3.65),style="slash")
        
}


makeDataFrame_all <- function(NsOrNo){
        
        atcg.mat <- as.data.frame(matrix(data = 0, ncol = nrow(modcoef), nrow = 4))
        #ATCG elements
        diag(atcg.mat[1:4, 1:4]) <- 1
        #Reserve the first column for intercept
        atcg.mat[,1] <- 1
        #synonymous or nonsynonymous mutation?
        atcg.mat[,6] <- NsOrNo
        #nonysynonymous interactions with a t c g
        atcg.mat[7:10] <- atcg.mat[,2:4] * atcg.mat[,6]
        
        #names
        names(atcg.mat) <- c(rownames(modcoef) )
        
        return(as.matrix(atcg.mat))
} 



plotMFs <- function(NsOrNo, colVal,offsetval){
        pchpar <- 16
        cexpar <- .6
        datInds <- which(d$Type == c("syn", "nonsyn")[NsOrNo + 1])
        
        xinds <- as.character(d$ref[datInds])
        for(i in 1:4){
                xinds[xinds == nucord[i]] <- i
        }
        xinds <- as.numeric(xinds)
        yinds <- d$mean[datInds]
        xjit <- rnorm(length(xinds), offsetval, .03) #jitter
        
        yinds[yinds == 0] <- 10^(-4.4)
        
        points(xinds + xjit, yinds, col = colVal,  pch = pchpar, cex = cexpar)
}

plotEsp <- function(NsOrNo, colVal, offset){
        
        pchpar <- 21
        cexpar <- 1
        arrowlwd <- 1.5
        arrowlen <-.08
        
        setUpDat <- makeDataFrame_all(NsOrNo)[,rownames(modcoef)]
        
        #compute
        #%*% = matrix manipulation
        pointToPlot <- exp(setUpDat %*% coef.vals)
        arrowToPlotLow <- exp(setUpDat %*% coef.pSE.vals)
        arrowToPlotHigh <- exp(setUpDat %*% coef.mSE.vals)
        
        arrows((1:4)+ offset, arrowToPlotLow, (1:4) + offset, arrowToPlotHigh, length = arrowlen, code = 3, angle = 90, lwd = arrowlwd, col = "black")
        points((1:4)+ offset, pointToPlot, col ="black", bg = substr(colVal,start=1, stop=7),  pch = pchpar, cex = .8)
        
}


