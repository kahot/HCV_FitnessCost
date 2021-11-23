# Plot the selection coefficient estimation results (Fig. 5)
library(ggplot2)
library(reshape2)
library(sfsmisc)
library(colorspace)
library(cowplot)
source("Rscripts/label_scientific.R")
colors2<-qualitative_hcl(6, palette="Dark3")
scaleFUN <- function(x) sprintf("%.2f", x)
col2_light<-qualitative_hcl(6, palette="Set3")



#### Read the SC data
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


# 1. Fig 5A: plot across the genome ###
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




pdf("Output/SelCoeff/SC_acrossGenome2.pdf",width=14,height=5)


## Create the outline
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(EstSC~pos, data=SCs,t="n",log='y',yaxt='n',xlab='Genome position',ylab="Estimated selection coefficient",
     ylim=c(ylow1,yhigh),xlim=c(265,8618))
eaxis(side = 2, at = 10^((0):(-(5))), cex=2)
for(k in 1:5){abline(h = 1:10 * 10^(-k), col = "gray80")}

for (i in 1:endnuc){
    if (is.na(SCs$Type[i])==T) next
    if (SCs$Type[i]=="stop") next
    if (SCs$Type[i]=="syn") {c=scCol[2]}
    if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("c","g")) {c=scCol[1]}
    if (SCs$Type[i]=="nonsyn"&SCs$ref[i]%in%c("a","t")) {c=scCol[3]}
    points(SCs$pos[i],SCs$EstSC[i],pch=21,col='gray30',lwd=0.3, bg=paste0(c,"99"),cex=.5)
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
lines(ns.roll50~pos,data=SCs, col="white",lwd=1.3)

lines(ns.roll50~pos,data=SCs, col=scColors2[1],lwd=1)

lines(syn.roll50~pos,data=SCs, col="white",lwd=1.3)
lines(syn.roll50~pos,data=SCs, col=scColors2[2],lwd=1)

abline(v=genes$end, col="gray80", lwd=.5)
abline(v=genes$end, col="gray80", lwd=.5)
#Add legend
legpos=300; legposV=0.032
rect(legpos, 0.42*legposV, (legpos+1000), 1.05*legposV, density = NULL, angle = 45,col=alpha("white",1))
points((legpos+100),legposV*0.8,pch=21,bg=scCol[2],col=1,cex=1)
text((legpos+150),legposV*0.8,"Syn",adj=0, cex=1)
points((legpos+100),legposV*0.64,pch=21,bg=scCol[1],col=1,cex=1)
text((legpos+150),legposV*0.64,"Nonsyn, A/T",adj=0, cex=1)
points((legpos+100),legposV*0.51,pch=21,bg=scCol[3],col=1,cex=1)
text((legpos+150),legposV*0.51,"Nonsyn, C/G",adj=0, cex=1)

box()




#dev.off()









###########
#Plot summary of 1) sel coef by type by nucleotide       
#2 no CPG creating mutations
sc3<-sc[sc$makesCpG==0,]
k=1
transSC<-list()
for (i in c("a","t","c","g")) {
    for (type in c("syn","nonsyn")){
        datavector<-sc3$EstSC[sc3$Type==type & sc3$ref==i]
        nt<-toupper(i)
        vname<-paste0(nt,".",type)
        dat<-data.frame(base=rep(nt,times=length(datavector)),
                        type=rep(type, times=length(datavector)), S.C.=datavector)
        transSC[[k]]<-dat
        names(transSC)[k]<-vname
        k=k+1
    }
}

scdata2<-do.call(rbind, transSC)

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 

scdata2$base<-factor(scdata2$base, levels=c("A","T", "C", "G"))
scdata2$type<-factor(scdata2$type, levels=c("syn","nonsyn"))

p2<-ggplot(scdata2,aes(x=base, y=S.C., fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
    scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
    labs(x="Nucleotide",y="Estimated selection coefficient")+
    scale_fill_manual(values=colors2[c(5,1)]) + theme_bw()+
    theme(legend.title = element_blank())+theme(axis.text.x = element_text(size =10, color=1), axis.title.x = element_blank())+
    theme(axis.text.y = element_text(size =10), panel.grid.major.x=element_blank())+
    geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
    scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))

#ggsave("Output/SelCoeff/SC.byNT_noCpG.pdf", width = 6,height = 4)


#Add the mutation frequency summary plot for contrast (Fig. 5C)
CVFiles3<-list.files("Output1A/Overview3/",pattern="overview3.csv")
s<-length(HCVFiles3)

Ts <-read.csv("Output/MutFreq/Filtered.Ts.Q35.csv",stringsAsFactors = F,row.names=1)
#coding region only
Ts<-Ts[Ts$pos>=342,c("pos","ref","mean","Type","makesCpG") ]
Ts<-Ts[Ts$makesCpG==0,]

#Plot summary of mutation frequency by type by nucleotide       
k=1
transMF<-list()
for (i in c("a","t","c","g")) {
    for (type in c("syn","nonsyn")){
        datavector<-Ts$mean[Ts$Type==type & Ts$ref==i]
        nt<-toupper(i)
        vname<-paste0(nt,".",type)
        dat<-data.frame(base=rep(nt,times=length(datavector)),
                        type=rep(type, times=length(datavector)), MF=datavector)
        transMF[[k]]<-dat
        names(transMF)[k]<-vname
        k=k+1
    }
}

mfdata<-do.call(rbind, transMF)
mfdata$base<-factor(mfdata$base, levels=c("A","T", "C", "G"))
mfdata$type<-factor(mfdata$type, levels=c("syn","nonsyn"))

z=rep(c(0.7,0.3),times=4)
x<-1:4
ybreaks<- c(1:10 * 10^c(-4),1:10 * 10^c(-3),1:10 * 10^c(-2),1:10 * 10^c(-1)) 
col2_light1<-qualitative_hcl(6, palette="Set2")

p3<-ggplot(mfdata,aes(x=base, y=MF, fill=factor(type)))+geom_boxplot(outlier.alpha =.4, outlier.color = "gray60")+
    scale_y_continuous(trans = 'log10',breaks=c(0.0001,0.001,0.01), minor_breaks=ybreaks, labels=label_scientific2)+
    labs(x="",y="Mutation frequency")+
    scale_fill_manual(values=col2_light1[c(5,1)]) + theme_bw()+
    theme(legend.title = element_blank()) +theme(axis.text.x = element_text(size =10, color=1))+
    theme(axis.text.y = element_text(size =10), axis.title.y= element_text(size =12), panel.grid.major.x=element_blank())+
    geom_vline(xintercept = c(1:3)+0.5, color="gray60")+
    scale_x_discrete(breaks=c("A","T","C","G"),labels=c(expression(A%->%G),expression("T"%->%C),expression(C%->%"T"),expression(G%->%A)))
#ggsave("Output/MutFreq/MF.byNT.noCpG.pdf", width = 5,height = 4)

library(gridBase)
library(grid)

plot.new()          
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(0,1,0,1)) ## create new vp with margins, you play with this values 
print(p2, p3, vp=vp1)

plot.new()          
vps <- baseViewports()
pushViewport(vps$figure) 
vp2 <-plotViewport(c(0,1,0,1)) ## create new vp with margins, you play with this values 
print(p3, vp=vp2)
