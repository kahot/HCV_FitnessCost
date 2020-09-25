library(ggplot2)
library(colorspace)
library(purrr)

cols2<-c("#66CCEE","#EE667799" ,"#22883399")
colors2<-qualitative_hcl(6, palette="Dark3")

# mut freq
mutfreq<-read.csv("Output1A/MutFreq.filtered/Filtered.Ts.Q35.csv", stringsAsFactors = F, row.names = 1)

#NCBI HCV1A sequences + 195 consensus (423 +195 = 618 sequences, coding sequnces only)
entropy<-read.table(paste0("Data/1A_combind_Logo.txt"))
colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
entropy$pos<-342:(nrow(entropy)+341)

compare<-merge(entropy,mutfreq, by="pos")
results<-cor.test(compare$Entropy,compare$mean, method = "spearman")
print(results)
#      rho 
#-0.677469 
#P < 2.2e-16

pdf(paste0("Output1A/SummaryFig/Entropy-MF_618_new.pdf"))
plot(compare$mean,-(compare$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
mtext("Within-host diversity (ave. mutation frequency)",1,2.2, cex=1.2)
mtext("Between-host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(compare$Entropy)~compare$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)

text(x=max(compare$mean, na.rm=T)-0.003,y=-1.2, labels=expression(paste(rho, " = 0.677***")),cex=1.1)
dev.off()



comparison2<-merge(entropy,mvf, by="pos")
results<-cor.test(comparison2$Entropy,comparison2$mean, method = "spearman")
print(results)
#rho 
#-0.7875385 

pdf(paste0("Output1A/SummaryFig.Filtered/Entropy-MVF.618.pdf"))
plot(comparison1$mean,-(comparison1$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
mtext("In vivo diversity (mean mutation frequency)",1,2.2, cex=1.2)
mtext("Among host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(comparison1$Entropy)~comparison1$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)
if (results[[3]]>=0.05) star<-""
if (results[[3]]<0.05&results[[3]]>=0.01) star<-"*"
if (results[[3]]<0.01&results[[3]]>=0.001) star<-"**"
if (results[[3]]<0.001) star<-"***"
paste0(rho,star)
text(x=max(comparison1$mean, na.rm=T)-0.015,y=-1.2, labels=expression(paste(rho, " = 0.788***")),cex=1.1)
dev.off()


#195 vs 618 entropy correlation:
ent195<-read.table("Data/HCV1A_logo_data.txt")
colnames(ent195)<-c("pos","A","C","G","T","Entropy195","Low","High","Weight")
ent195$pos<-342:(nrow(ent195)+341)
ent<-merge(ent195, compare[,c('pos','Entropy')], by='pos')

results<-cor.test(ent$Entropy,ent$Entropy195, method = "spearman")
print(results)
#p-value < 2.2e-16
#rho 
#-0.210369 
plot(-ent$Entropy195,-ent$Entropy, ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
