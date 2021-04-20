library(zoo)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggthemes)
library(sfsmisc)
library(colorspace)
library(cowplot)

source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")
scColors<-c("#AEE5F6","#4477AA","#EEBAB9","#EE6677","#228833")

colors2<-qualitative_hcl(6, palette="Dark3")
col2_light<-qualitative_hcl(6, palette="Set3")
scColors2<-c("#EC4A4D","#0055EC")



###########
df<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
#coding regions only
df<-df[df$pos>=342,]

#add the gene info
genes<-read.csv("Data/HCV_annotations2.csv",stringsAsFactors = F)
genes$Gene[6]<-"NS1"
gene.vector<-c()
for (i in 1:(nrow(genes)-1)){
        gene.vector<-c(gene.vector, rep(genes$Gene[i],times=genes$start[i+1]-genes$start[i]))
}
genetable<-data.frame("pos"=c(1:length(gene.vector)))
genetable$gene<-gene.vector
end<-df$pos[nrow(df)]
genetable<-genetable[genetable$pos>=342&genetable$pos<=end,]

sc<-merge(df, genetable, by="pos")

endnuc<- sc$pos[nrow(sc)]  # from 11.Mut.Freq calculation. If using the last row -> Trans$pos[nrow(Trnas)]
n<-data.frame("pos"=c(342:endnuc))

SCs<-merge(n,sc,by="pos",all.x=T) 
range(SCs$EstSC,na.rm = T) #0.0001257571 0.0101180076

## Plot across the genome ###
# create a table for syn only and nonsyn only to calculate separate rolling average
ns<-SCs[SCs$Type=="nonsyn",]
syn<-SCs[SCs$Type=="syn",]
ns<-merge(n, ns, by="pos",all.x=T) 
syn<-merge(n, syn, by="pos",all.x=T) 

ns.roll50<-rollmean(ns$EstSC, k=50, na.rm=T, align="center")
SCs$ns.roll50<-c(rep(NA, times=25),ns.roll50,c(rep(NA, times=24)))

syn.roll50 <-rollmean(syn$EstSC, k=50, na.rm=T,align="center")
SCs$syn.roll50<-c(rep(NA, times=25),syn.roll50,c(rep(NA, times=24)))

# plot syn and nonsyn with different colors

ylow1=0.00008
yhigh=0.025
pdf("Output/SelCoeff/SC_acrossGenome.pdf",width=14,height=5)
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow1,yhigh),xlim=c(265,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

for (i in 1:endnuc){
        if (is.na(SCs$Type[i])==T) next
        if (SCs$Type[i]=="stop") next
        if (SCs$Type[i]=="syn") {c=colors2[5]}
        if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("c","g")) {c=colors2[1]}
        if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("a","t")) {c=colors2[6]}
        points(SCs$pos[i],SCs$EstSC[i],pch=21,col='gray30',lwd=0.3, bg=paste0(c,"B3"),cex=.4)
}

ylow<-0.000064
for (j in 2:(nrow(genes)-1)){
        xleft<-genes$start[j]
        xright<-genes$start[j+1]
        
        if ((j==4|j==6|j==9)){
                rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+80, 1.44*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
                #mtext(paste0(genes$Gene[j]),side= 1, line=-0.1, at= xleft+80, col="black", cex=0.8)
        }
        else if (j==12){
                rect(xleft,ylow,genes$end[j],1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xleft+600,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)
        }
        else{rect(xleft,ylow,xright,1.3*ylow,density = NULL, angle = 45,col="white",border ="gray40")
                text(xright-(xright-xleft)/2,1.15*ylow,paste0(genes$Gene[j]),col="black", cex=0.8)}
}


# roling average of 50

lines(ns.roll50~pos,data=SCs, col=scColors2[1],lwd=1)
lines(syn.roll50~pos,data=SCs, col=scColors2[2],lwd=1)

abline(v=genes$end, col="gray80", lwd=.5)
abline(v=genes$end, col="gray80", lwd=.5)
#Add legend
legpos=300; legposV=0.032
rect(legpos, 0.42*legposV, (legpos+1000), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.8,pch=21,bg=colors2[5],col=1,cex=1)
text((legpos+150),legposV*0.8,"Syn",adj=0, cex=1)
points((legpos+100),legposV*0.64,pch=21,bg=colors2[6],col=1,cex=1)
text((legpos+150),legposV*0.64,"Nonsyn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.51,pch=21,bg=colors2[1],col=1,cex=1)
text((legpos+150),legposV*0.51,"Nonsyn, C/G",adj=0, cex=1)

box()
dev.off()


