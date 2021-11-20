# Create the Distribution of Fitness Effects plot (Fig. S6)
library(ggplot2)
library(colorspace)
library(cowplot)


#Read the selection coefficient estiamtes file
sc<-read.csv("Output/SelCoeff/SC.csv", stringsAsFactors = F, row.names = 1)
sc<-sc[df$pos>=342,]

colors2<-qualitative_hcl(6, palette="Dark3")
scCols<-c("#E16A86","#009ADE")

## A

dt1<-sc[sc$ref=="a",]
A<-ggplot(dt1, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt1, Type=="nonsyn"),  color="gray60", alpha = 0.5,fill=scCols[1])+
        geom_histogram(data=subset(dt1, Type=="syn"),color="gray60", alpha = 0.5, fill=scCols[2])+
        scale_fill_manual(values=scCols)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        theme_bw()+ylab("Number of sites")+xlab("")+ggtitle("A")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt1$EstSC[dt1$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt1$EstSC[dt1$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)+
        theme(legend.position = "none")

dt3<-sc[sc$ref=="t",]
T_<-ggplot(dt3, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt3, Type=="nonsyn"),  color="gray60", alpha = 0.5,fill=scCols[1])+
        geom_histogram(data=subset(dt3, Type=="syn"),color="gray60", alpha = 0.5, fill=scCols[2])+
        scale_fill_manual(values=scCols)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        theme_bw()+ylab("Number of sites")+xlab("")+ggtitle("T")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt3$EstSC[dt3$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt3$EstSC[dt3$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)+
        theme(legend.position = "none")

dt4<-sc[sc$ref=="c" & sc$Type!='stop',]
mean1<-mean(dt4$EstSC[dt4$Type=="syn"], na.rm = T)
mean2<-mean(dt4$EstSC[dt4$Type=="nonsyn"], na.rm = T)

C<-ggplot(dt4, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt4, Type=="nonsyn"),  color="gray60", alpha = 0.5,fill=scCols[1])+
        geom_histogram(data=subset(dt4, Type=="syn"),color="gray60", alpha = 0.5, fill=scCols[2])+
        
        scale_fill_manual(values=scCols)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        theme_bw()+ylab("Number of sites")+xlab("")+ggtitle("C")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt4$EstSC[dt4$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt4$EstSC[dt4$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)+
        theme(legend.position = "none")


dt2<-sc[sc$ref=="g"&sc$Type!='stop',]

G<-ggplot(dt2, aes(x=EstSC, fill=Type))+
        geom_histogram(data=subset(dt2, Type=="nonsyn"), aes(x=EstSC, fill=Type), color="gray60", alpha = 0.5)+
        geom_histogram(data=subset(dt2, Type=="syn"),aes(x=EstSC,fill=Type), color="gray60", alpha = 0.5)+
        scale_x_continuous(trans = 'log10',label=label_scientific2, limits = c(0.00005,.1))+
        scale_fill_manual(name='', values=scCols,labels=c("Nonsyn","Syn"))+
        scale_color_manual(name='', values=scCols,labels=c("Nonsyn","Syn"))+
        theme_bw()+ylab("Number of sites")+xlab("")+ggtitle("G")+theme(plot.title = element_text(hjust = 0.5, size=12))+
        geom_vline(aes(xintercept=mean(dt2$EstSC[dt2$Type=="syn"], na.rm = T)),color="blue", linetype="dashed", size=.5)+
        geom_vline(aes(xintercept=mean(dt2$EstSC[dt2$Type=="nonsyn"], na.rm = T)),color="red", linetype="dashed", size=.5)


plot_grid(A,T_,C,G, nrow = 1, ncol=4, rel_widths=c(1,1,1,1.5))
#grid.arrange(A,T_,C,G,  nrow = 1)
ggsave("Output/SelCoeff/SC_DFEplot.pdf", width=9, height=2.5)
