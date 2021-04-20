source("Rscripts/baseRscript.R")
require('nloptr')

DF<-read.csv("Output/SelCoeff/SC.csv", row.names = 1, stringsAsFactors = F) 
DF<-DF[DF$pos>=342,]

geller<-read.csv("Output/Geller/Geller.MutRates.Summary_updated.csv", stringsAsFactors = F, row.names = 1)
geller<-geller[!(geller$Mutation=="AA"|geller$Mutation=="UU"|geller$Mutation=="CC"|geller$Mutation=="GG"),]
geller$Mutation<-gsub("U","T",geller$Mutation)


runSub <- function(DF, plotMe = TRUE){
        freqColName <- "mean"
        process <- function(cnames){
                toReturn <- c()
                for(i in 1:length(cnames)){
                        toReturn[i] <- strsplit(cnames[i], "V")[[1]][2]
                }
                return(as.numeric(toReturn))
        }
        
        ss<-DF$TSmutrate/DF[,freqColName]
        
        #Now, we need to fit a gamma distribution
        
        #http://stats.stackexchange.com/questions/160304/fit-gamma-distribution-to-dataset-with-r
        GammaNLL <- function(pars, data){
                alpha <- pars[[1]]
                theta <- pars[[2]]
                return (-sum(dgamma(data, shape = alpha, scale = theta, log = TRUE)))
        }
        
        
        require(nloptr) #non-linear optimization
        Fit <- nloptr(x0 = c(1, 1), eval_f = GammaNLL, lb = c(0,0), data = ss, opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e15, xtol_rel = 1e-15,  xtol_abs = 1e-15))
        
        #Fit <- nloptr(x0 = c(1, 1), eval_f = GammaNLL, lb = c(0,0), data = ss,opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e15, xtol_rel = 1e-10,  xtol_abs = 1e-10))
        
        return( list(shape = Fit$solution[1], scale = Fit$solution[2], ests = ss) )
}


subDat <- DF
yrange <- 1:nrow(DF)
runit <- runSub(subDat)

obs.scale <- runit$scale
obs.shape <- runit$shape

ptable <- data.frame(matrix(nrow = 3, ncol  = 3))
colnames(ptable)<-c("no_of_sites","Scale_parameter","Shape_parameter")
rownames(ptable)<-c("Parameter","upCI","lowCI")
ptable[1, 1] <- nrow(subDat)
ptable[1, 2] <- round(obs.scale, 6)
ptable[1, 3] <- round(obs.shape, 3)

numreps <- 1000
btstrp.vals <- matrix(data = NA, nrow = numreps, ncol = 2)

for(k in 1:numreps){
        subDat.bs <- subDat[sample(yrange, length(yrange), replace = TRUE),]
        runit.bs <- runSub(subDat.bs, plotMe = FALSE)
        btstrp.vals[k,] <- c(runit.bs$shape, runit.bs$scale)
}

write.csv(btstrp.vals, "Output/SelfCoeff/Gamma.btstrp.values.csv")

#Create confidence intervals
shape.CI <- quantile(btstrp.vals[,1], c(.025, .975) )
scale.CI <- quantile(btstrp.vals[,2], c(.025, .975) )

ptable[2:3, 2] <- round(scale.CI, 6)
ptable[2:3, 3] <- round(shape.CI, 3)

write.table(ptable, "Output/SelfCoeff/Gamma.parameters.csv")
