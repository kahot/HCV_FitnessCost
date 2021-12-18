#Simulate mutation freq at each site based on the mutation rate and fitness cost (s)
#Adapted from Theys et al 2018 paper

library(plotrix)
library(FSA)
library(colorspace)
source("Rscripts/Pcorrection.R")

#Read the data that have observed frequencies (and transpose it)
Ts<-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv", row.names = 1)
Ts<-Ts[Ts$pos>=342,] #coding regions only 
Ts<-Ts[,c(which(colnames(Ts)=="pos"), 1:195)]
TsFreq<-data.frame(t(Ts))
colnames(TsFreq)<-TsFreq[1,]
TsFreq<-TsFreq[-1,]

#Read the data that have estiamted SC and mutation rates
SC<-read.csv("Output/SelCoeff/SC.csv", row.names = 1)

#Color assignment for plots (syn. vs. nonsyn vs. nonsense)
colo2<-qualitative_hcl(6, palette="Dark3")
SC$Color<-colo2[1]
SC$Color[SC$Type=="syn"]<-colo2[5]
SC$Color[SC$Type=="stop"]<-colo2[2]


######
#Get the sample sizes for each patient (average read depth) and store them in readDepth
Reads<-read.csv("Output/ReadDepth_All.csv", row.names = 1)
Reads<-Reads[Reads$pos>=342,]
readDepth<-colMeans(Reads[,2:196], na.rm=T)
readDepth<-round(unname(readDepth))

#dir.create("Output/Simulation/")

#set the population size
Ne=100000

listPvalues<-c()
simulatedValues<-data.frame(pos=SC$pos)

# run first with raw p-values. Then run again, adding the adjusted p-value to the figures.
run=1
pdf ("Output/Simulation/comparisonData.vs.Simulations_10000_all.pdf") 

#Add the adjusted p-values to the figures,
run=2
pdf ("Output/Simulation/comparisonData.vs.Simulations_10000_adjustedPval2.pdf") 
simRe<-read.csv("Output/Simulation/Simulation100k_results_summary_statistics.csv", row.names = 1)


for (i in 1:nrow(SC)){    
        system("rm ./SimData/Data*")
        site<-SC$pos[i]
        print(site)

        realFreqs<-TsFreq[,paste(site)]
        #Mut rates and sel coefficients
        mu=round(SC$TSmutrate[SC$pos==site],8)
        print(mu)
        cost=round(SC$EstSC[SC$pos==site],6)
        print(cost)
        
        numoutputs = 195 #(max num of patients)
        theta = mu*Ne
        
        if (TRUE){
            #create data through simulations
            seed =100
            #make script 
            x<-"#!/bin/bash"
            x<-c(x,paste("mu=",mu,sep=""))
            outputfrequency=min(c(2*Ne,ceiling(5/cost)))
            x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
            x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
            x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
            x<-c(x,paste("cost=",cost,sep=""))
            x<-c(x,paste("for seed in",seed))
            
            if (Ne == 5000) sentence = paste("\" | ./Scripts/Viralevolution_5000   >./Output/Simulation/SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            if (Ne == 50000)sentence = paste("\" | ./Scripts/Viralevolution_50000 >./Output/Simulation/SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            if (Ne == 100000)sentence = paste("\" | ./Scripts/Viralevolution_100000 >./Output/Simulation/SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
            
            x<-c(x,"do",
                 "echo \"", "$seed", "$mu", "$cost",
                 "$output_every_Xgen", "$numgen_inN", "$start_output", sentence, 
                 "done")
            
            write(x,file="./Scripts/tempscript.sh")
            system("chmod 775 Scripts/tempscript.sh") #giving the permission maybe tricky (run this on terminal)
            system("./Scripts/tempscript.sh")     #Run tempscript.sh 
               # if an error '/usr/local/lib/libgsl.0.dylib' (no such file) shows up-> created sim link to the most updated file "ln -s libgsl.23.dylib libgsl.0.dylib"
        }
        #Read the data
        filename=paste("Output/Simulation/SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
        Freqs<-read.csv(filename,sep="\t")$freq
        system("rm Output/Simulation/SimData/Data*")
        
        #Get real sample Freqs (sample from pop freqs)
        #use the real n 
        numoutputs<-length(which(!is.na(realFreqs)))
        #sample simulated frequencies based on the sample depth
        simFreqs<-rbinom(numoutputs,readDepth[!is.na(realFreqs)],Freqs)/readDepth[!is.na(realFreqs)]
        realFreqs<-realFreqs[!is.na(realFreqs)]
        
        #create a summary table
        simulatedValues$numberofZero[i]<-length(simFreqs[simFreqs==0])
        simulatedValues$Median[i]<-median(simFreqs)
        simulatedValues$Mean[i]<-mean(simFreqs)
        simulatedValues$Max[i]<-max(simFreqs)
        simulatedValues$RealnumberofZero[i]<-length(realFreqs[realFreqs==0])
        simulatedValues$RealMedian[i]<-median(realFreqs)
        simulatedValues$RealMean[i]<-mean(realFreqs)
        simulatedValues$RealMax[i]<-max(realFreqs)
        simulatedValues$numberofSample[i]<-numoutputs
        
        #run a statistical test and obtain p-values
        pval<-wilcox.test(simFreqs,realFreqs)[[3]]
        simulatedValues$rawP[i]<-pval
        listPvalues<-c(listPvalues,pval)
        
        #Plot
        br<-seq(-10^-7,1.01,by=0.01)
        br2<-seq(0,1,by=0.02)
        xlims=c(0,0.01+max(c(simFreqs,realFreqs),na.rm=TRUE))
        ylims = c(0,max(length(which(realFreqs<0.01)),length(which(simFreqs<0.01))))
        
        par(mfrow=c(2,1))
        maxyheight=50

        hist(rep(0,maxyheight),breaks=br,xlim=xlims,  
             ylab="", xlab = "", main=paste("Single-SFS, site",site),
             xaxt="n", yaxt="n",col=SC$Color[SC$pos==site])
        hist(c(rep(0,min(maxyheight-10,length(which(realFreqs<0.01)))),realFreqs[which(realFreqs>=0.01)]), breaks=br, col=SC$Color[SC$pos==site], add=TRUE)
        axis(2,labels = c(0,10,20,30,40,max(maxyheight,length(which(realFreqs<0.01)))), 
             at = c(0,10,20,30,40,maxyheight), las=1)
        if (length(which(realFreqs<0.01))>=maxyheight){
            axis.break(axis=2,breakpos=maxyheight-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
            rect(-0.001, maxyheight-12, 0.0101, maxyheight-8,col="white",border="white")
        }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
        
        axis(1, at=br2,labels=br2, col.axis="black", las=2)
        mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
        mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
        abline(v=mean(realFreqs,na.rm=TRUE),col=2)
        
        hist(rep(0,maxyheight),breaks=br,xlim=xlims, #ylim = ylims, 
             ylab="", xlab = "", main=paste("Simulated data, cost",round(cost,4)),
             xaxt="n", yaxt="n",col="pink")
        hist(c(rep(0,min(maxyheight-10,length(which(simFreqs<0.01)))),simFreqs[which(simFreqs>=0.01)]), breaks=br, col="pink", add=TRUE)
        axis(2,labels = c(0,10,20,30,40,max(maxyheight,length(which(simFreqs<0.01)))), 
             at = c(0,10,20,30,40,maxyheight), las=1)
        if (length(which(simFreqs<0.01))>=maxyheight){
            axis.break(axis=2,breakpos=maxyheight-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
            rect(-0.001, maxyheight-12, 0.0101, maxyheight-8,col="white",border="white")
        }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
        
        if (run==2) text(x=mean(xlims),y=maxyheight/2,paste("corrected P-value=",round(simRe$Holm[simRe$pos==site])))
        if (run==1) text(x=mean(xlims),y=maxyheight/2,paste("corrected P-value=",round(pval)))

        axis(1, at=br2,labels=br2, col.axis="black", las=2)
        mtext("Observed frequency", side=1, line=4, cex.lab=1,las=1, col="black")
        mtext("Number of patients", side=2, line=2.5, cex.lab=1,las=3, col="black")
        abline(v=mean(simFreqs),col=2)
        
}
dev.off()

#add correction to p-values
simRe<-Pcorrection(simulatedValues)
#based on Holm's correction, how many are significant/non-significant?
simRe$Significance<-ifelse(simRe$Holm<0.05, "Y","N")

aggregate(simRe[,c(2:10)], by=list(simRe$Significance), mean)
#  Group.1 numberofZero      Median        Mean        Max  RealMedian    RealMean    RealMax numberofSample
#1       N     12.16319 0.003373085 0.004643216 0.02431215 0.003020725 0.004649867 0.05660691       141.2133
#2       Y     17.10122 0.005812292 0.007223891 0.03178292 0.003412552 0.007132051 0.09283524       146.6736

nrow(simRe[simRe$Significance=="Y",]) #588 
588/7957 #7.39%

write.csv(simRe, "Output/Simulation/Simulation100k_results_summary_statistics.csv")



### Create example of figures for the manuscript
simRe<-read.csv("Output/Simulation/Simulation100k_results_summary_statistics.csv", row.names = 1)

#find significant sites 
simRe<-merge(simRe, SC[,c("pos","Type","ref")], by="pos")
ex<-simRe[simRe$Type=="syn"&simRe$Significance_holm=="Y",]
ex2<-simRe[simRe$Type=="stop"&simRe$rawP>0.05,]

# Picked 5528:5531
simF<-list()
realF<-list()
positions<-paste(5528:5531) # for nonsignificant sites
positions<-c("587","629") #for significant sites
summary<-data.frame(pos=positions)

#Run the simulations for those sites only and save the simulated frequencies
Ne=100000
for (i in 1:length(positions)){
    site<-positions[i]
    print(site)
    realFreqs<-TsFreq[,paste(site)]
    
    #Mut rates and sel coefficients
    mu=round(SC$TSmutrate[SC$pos==site],8)
    cost=round(SC$EstSC[SC$pos==site],6)
    
    numoutputs = 195
    theta = mu*Ne
    
    if (TRUE){
        #create data through simulations
        seed =100
        #make script 
        x<-"#!/bin/bash"
        x<-c(x,paste("mu=",mu,sep=""))
        outputfrequency=min(c(2*Ne,ceiling(5/cost)))
        x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
        x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
        x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
        x<-c(x,paste("cost=",cost,sep=""))
        x<-c(x,paste("for seed in",seed))
        
        sentence = paste("\" | ./Scripts/Viralevolution_100000 >./Output/Simulation/SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
        x<-c(x,"do",
             "echo \"", "$seed", "$mu", "$cost",
             "$output_every_Xgen", "$numgen_inN", "$start_output", sentence, 
             "done")
        
        write(x,file="./Scripts/tempscript.sh")
        system("chmod 775 Scripts/tempscript.sh")
        system("./Scripts/tempscript.sh")     #Run tempscript.sh 
    }
    #Read the data
    filename=paste("Output/Simulation/SimData/Data_T_", theta, "_cost_", cost,".txt",sep="")
    Freqs<-read.csv(filename,sep="\t")$freq
    #system("rm Output/Simulation/SimData/Data*")
    
    #Get sample Freqs (sample from pop freqs)
    numoutputs<-length(which(!is.na(realFreqs)))
    simFreqs<-rbinom(numoutputs,readDepth[!is.na(realFreqs)],Freqs)/readDepth[!is.na(realFreqs)]
    realFreqs<-realFreqs[!is.na(realFreqs)]
    realF[[i]]<- realFreqs
    simF[[i]]<-simFreqs
    summary$Median[i]<-median(simFreqs)
    summary$Mean[i]<-mean(simFreqs)
    summary$RealMedian[i]<-median(realFreqs)
    summary$RealMean[i]<-mean(realFreqs)
    summary$numberofSample[i]<-numoutputs
}

summary<-merge(summary, SC[,c("pos","EstSC","Type","Color","ref")], by="pos")
summary<-merge(summary, simRe[,c("pos","Holm")], by="pos")    
#Save the summary
write.csv(summary, "Output/Simulation/ExampleSites5528.csv")

# Make plots
pdf("Output/Simulation/Examples_plots5528.pdf", width = 4, height = 5)

positions<-paste(5528:5531)
b=0.005
br<-seq(-10^-7,1.01,by=b)
br2<-seq(0,1,by=0.05)
xlims=c(0,0.16)
maxyheight=50
par(mfrow=c(2,1))
for (i in 1: length(positions)){
    simf<-simF[[i]]
    realf<-realF[[i]]
    if (i==1) title<-"Synonymous site: G5528A"
    if (i==2) title<-"Nonsense site: C5529T"
    if (i==3) title<-"Nonsynonymous site:A5530G"
    if (i==4) title<-"Synonymous site: G5531A"
    
    par(mar = c(1, 4, 4, 2))
    
    ylims = c(0,max(length(which(realf<b)),length(which(simf<b))))
    
    hist(c(rep(0,min(maxyheight,length(which(realf<b)))),realf[which(realf>=b)]),
         ylab="", xlab = "", main=title, xaxt="n", yaxt="n",xlim=xlims,
         breaks=br, col=summary$Color[i])
    axis(2,labels = c(0,10,20,30,40,max(maxyheight,length(which(realf<b)))), 
         at = c(0,10,20,30,40,maxyheight), las=1)
    if (length(which(realFreqs<b))>=maxyheight){
        axis.break(axis=2,breakpos=maxyheight-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
        rect(-0.001, maxyheight-11, 0.15, maxyheight-9,col="white",border="white")
    }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
    
    axis(1, at=br2,labels=br2, col.axis="black", las=2)
    text(x=max(xlims)/2, y=maxyheight/2+10,paste("Observed data"))
    text(x=max(xlims)/2, y=maxyheight/2+2, col="gray20", paste0("s=", round(summary$EstSC[i],5)),cex=0.9)
    abline(v=mean(realFreqs,na.rm=TRUE),col=2)
    
    par(mar = c(4, 4, 1, 2))
    hist(c(rep(0,min(maxyheight,length(which(simf<b)))),simf[which(simf>=b)]), 
         breaks=br, xlim=xlims,col="pink",  ylab="", xlab = "",main='',xaxt="n", yaxt="n")
    axis(2,labels = c(0,10,20,30,40, max(maxyheight,length(which(simf<b)))), 
         at = c(0,10,20,30,40,maxyheight), las=1)
    if (length(which(simFreqs<b))>=maxyheight){
        axis.break(axis=2,breakpos=maxyheight-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
        rect(-0.001, maxyheight-11, 0.22, maxyheight-9,col="white",border="white")
    }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
    
    axis(1, at=br2,labels=br2, col.axis="black", las=2)
    text(x=max(xlims)/2,y=maxyheight/2+10,paste("Simulated data"))
    text(x=max(xlims)/2,y=maxyheight/2+2,col="gray20",paste0("adjusted P-value=",round(summary$Holm[i],2)),cex=0.9)
    abline(v=mean(simFreqs),col=2)
}
dev.off()



## Significant sites
pdf("Output/Simulation/Examples_plotsSig.pdf", width = 4, height = 5)
positions<-c("587","629")
b=0.005
br<-seq(-10^-7,1.01,by=b)
br2<-seq(0,1,by=0.05)
xlims=c(0,0.16)
maxyheight=50
par(mfrow=c(2,1))
for (i in 1: length(positions)){
    simf<-simF[[i]]
    realf<-realF[[i]]
    if (i==1) title<-"Synonymous site: T587C"
    if (i==2) title<-"Nonsense site: G629A"
    
    par(mar = c(1, 4, 4, 2))
    
    ylims = c(0,max(length(which(realf<b)),length(which(simf<b))))
    
   hist(c(rep(0,min(maxyheight,length(which(realf<b)))),realf[which(realf>=b)]),
         ylab="", xlab = "", main=title, xaxt="n", yaxt="n",xlim=xlims,
         breaks=br, col=summary$Color[i])
    axis(2,labels = c(0,10,20,30,40,max(maxyheight,length(which(realf<b)))), 
         at = c(0,10,20,30,40,maxyheight), las=1)
    if (length(which(realFreqs<b))>=maxyheight){
        axis.break(axis=2,breakpos=maxyheight-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
        rect(-0.001, maxyheight-11, 0.15, maxyheight-9,col="white",border="white")
    }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
    
    axis(1, at=br2,labels=br2, col.axis="black", las=2)
    text(x=max(xlims)/2, y=maxyheight/2+10,paste("Observed data"))
    text(x=max(xlims)/2, y=maxyheight/2+2, col="gray20", paste0("s=", round(summary$EstSC[i],5)),cex=0.9)
    abline(v=mean(realFreqs,na.rm=TRUE),col=2)
    
    par(mar = c(4, 4, 1, 2))
    hist(c(rep(0,min(maxyheight,length(which(simf<b)))),simf[which(simf>=b)]), 
         breaks=br, xlim=xlims,col="pink",  ylab="", xlab = "",main='',xaxt="n", yaxt="n")
    axis(2,labels = c(0,10,20,30,40, max(maxyheight,length(which(simf<b)))), 
         at = c(0,10,20,30,40,maxyheight), las=1)
    if (length(which(simFreqs<b))>=maxyheight){
        axis.break(axis=2,breakpos=maxyheight-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
        rect(-0.001, maxyheight-11, 0.22, maxyheight-9,col="white",border="white")
    }else{axis(2,labels = maxyheight-10,at=maxyheight-10,las=1)}
    
    axis(1, at=br2,labels=br2, col.axis="black", las=2)
    text(x=max(xlims)/2,y=maxyheight/2+10,paste("Simulated data"))
    text(x=max(xlims)/2,y=maxyheight/2+2,col="gray20",paste0("adjusted P-value=",round(summary$Holm[i],2)),cex=0.9)
    abline(v=mean(simFreqs),col=2)
}
dev.off()

