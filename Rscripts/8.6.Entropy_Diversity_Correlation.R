#Calculate the diversity levels at each site and compare between-host and within-host.

library(ggplot2)
library(colorspace)
library(purrr)

cols2<-c("#66CCEE","#EE667799" ,"#22883399")
colors2<-qualitative_hcl(6, palette="Dark3")

#Read the diversity (MVF) data
div<-read.csv("Output/MutFreq/Filtered.MVF.Q35.csv", stringsAsFactors = F, row.names = 1)
div<-div[,c("pos","mean")]

# Read the between-host diversity (-entropy) data file
# -contains NCBI HCV1A sequences + 195 consensus (423 +195 = 618 sequences, coding sequences only)
entropy<-read.table(paste0("Data/HCV1A_entropy.txt"))
colnames(entropy)<-c("pos","A","C","G","T","Entropy","Low","High","Weight")
entropy$pos<-342:(nrow(entropy)+341)
ent<-entropy[,c("pos","Entropy")]
compare<-merge(ent,div, by="pos")
results<-cor.test(-compare$Entropy,compare$mean, method = "spearman")
print(results)
#      rho 
#0.6908502 
#P < 2.2e-16

pdf(paste0("Output/SummaryFig/Entropy-MF.pdf"))
plot(compare$mean,-(compare$Entropy), ylab="", xlab="",
     pch=16,col=colors2[5],cex=0.6, cex.axis=1.2)
mtext("Within-host diversity (average MVF)",1,2.2, cex=1.2)
mtext("Between-host diversity (-Shannon's entropy)",2,2.2, cex=1.2)
abline(lm(-(compare$Entropy)~compare$mean), col = "gray70")
rho<-as.numeric(results[[4]])
rho<-format(round(rho,3), nsmall=3)

text(x=max(compare$mean, na.rm=T)-0.003,y=-1.2, labels=expression(paste(rho, " = 0.690***")),cex=1.1)
dev.off()

#### 
# Null hypothesis of correlation =1. Test the deviations od the differences from 0

#normalize both values for exact comparison
norm<-function(x){(x-min(x))/(max(x)-min(x))}

compare$nor.entropy<-norm(-compare$Entropy)
compare$nor.diversity<-norm(compare$mean)
compare$diff<-compare$nor.entropy-compare$nor.diversity
null<-rep(0, times=nrow(compare))
res<-wilcox.test(compare$diff,null)
print(res)

plot(compare$diff, pch=".")


##Compare with other genotypes
#Compare Div(1a) -Entropy(1b)
entropy1B<-read.table("Data/HCV1B_entropy.txt")
colnames(entropy1B)<-c("org.pos.1B","A","C","G","T","Entropy1B","Low","High","Weight")
entropy1B$org.pos.1B<-342:(nrow(entropy1B)+341)

positions<-read.csv("Data/MergedPositionInfo.csv", row.names = 1)

ent1b<-merge(entropy1B, positions, by="org.pos.1B")
compare2<-merge(compare, ent1b[,c("org.pos.1A","Entropy1B")], by.x="pos", by.y="org.pos.1A",)
cor.test(-compare2$Entropy1B,compare2$mean, method = "spearman")
#     rho 
#0.5817138 
#p-value < 2.2e-16


cor.test(-compare2$Entropy2,compare2$mean, method = "spearman")
#rho 
#0.5561916
#p-value < 2.2e-16

#Compare Div(1a) -Entropy(3a)
entropy3A<-read.table("Data/HCV3A_entropy.txt")
colnames(entropy3A)<-c("org.pos.3A","A","C","G","T","Entropy3A","Low","High","Weight")
entropy3A$org.pos.3A<-340:(nrow(entropy3A)+339)
ent3a<-merge(entropy3A, positions, by="org.pos.3A")
compare3<-merge(compare, ent3a[,c("org.pos.1A","Entropy3A")], by.x="pos", by.y="org.pos.1A",)
cor.test(-compare3$Entropy3A,compare3$mean, method = "spearman")

#S = 3.9702e+10, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#0.5268077 




ent1a3a<- read.csv("~/programs/HCV_project/Output_all/Diversity/Entropy_1A-3A.csv")
entropy3A<-ent1a3a[,c("merged.pos","org.pos.1A.x","org.pos.3A.x","Entropy1", "Entropy2")]

compare3<-merge(compare, entropy3A, by.x="pos", by.y="org.pos.1A.x",)

results<-cor.test(-compare3$Entropy2,compare3$mean, method = "spearman")
print(results)
#     rho 
#0.4977853
#p-value < 2.2e-16

