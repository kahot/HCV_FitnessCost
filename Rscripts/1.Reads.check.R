##Quality filter SAM files mapped to each sample's consensus.
#Reads are mapped separately for merged (me) and unmerged (un) files
#Sam files are avilable at Figshare (download them to 'Output/sam/').

library(stringr)
library(ape)
library(seqinr)
library(e1071)
library(msa)
library(car)
library(readtext)

dir.create("Output/HammingDistanceFiltering/")
suppressMessages(library(msa))


#read sam file
HCVsams<-list.files("Output/sam/",recursive = F,pattern="sam")
coln<-c('QNAME','Flag','RefName','Pos','MapQ','cigar','MRNM','Mpos','isize','seq','Qual','tag1','tag2','tag3','tag4','tag5','tag6')

ham.distance<-list()
filenames<-list()
sam_indels<-list()
ham.distance.indel<-list()
filtered_sam<-list()
sum<-list()

for (i in 1:length(HCVsams)){
        print(i)
        print(HCVsams[i])
        sam<-read.table(paste0("Output/sam/",HCVsams[i]),skip=3, col.names=coln, sep = "\t",fill=T, comment.char="",quote= "")
        sam<-sam[,1:11]
        
        #only the mapped reads
        sam<-subset(sam, MapQ>0&MapQ<61)
        sam$seq<-as.character(sam$seq)
        print(paste("# of mapped reads: ",nrow(sam)))
        mapped.reads<-nrow(sam)
        file.name<-substr(paste(HCVsams[i]),start=1,stop=7)
        fname2<-substr(paste(HCVsams[i]),start=1,stop=10)
        
        #hist(sam$MapQ,main=paste('Histogram of Mapping Quliaty for ',fname2))
        
        #read the consensus seq of each sample
        consensus<-read.dna(paste0("Output/Consensus/",file.name,"_consensus.fasta"), format = "fasta",as.character=TRUE)
        #cheange bases in the consensus sequence to upper case
        consensus<-as.character(sapply(consensus, function(v) {
                if (is.character(v)) return(toupper(v))
                else return(v)
                }))
        #replace ? with N
        consensus<-recode(consensus,'"?"="N"')

        #removed the rows containing indels as it substantially increases the hamming distance
        cigars<-as.character(sam$cigar)
        indel<-grepl("I|D|N",cigars)
        sam2<-sam[!indel,]  #reads without any indels
        no_indels<-nrow(sam2)
        print(paste("# of reads w/o indels: ", no_indels))

        #create a list with reads with indels only
        sam_indels<-sam[indel,]
        with_indels<-nrow(sam_indels)

        #calculate a hamming distance of each read (without indels) to its consensus
        H<-c()
        for (j in 1:nrow(sam2)) {
                reads<-sam2$seq[j]
                reads<-unlist(strsplit(reads,""))
                if ((sam2$Pos[j]+length(reads)-1)>length(consensus)) {
                        cut<-(sam2$Pos[j]+length(reads)-1)-length(consensus)
                        reads<-reads[1:(length(reads)-cut)]
                        ref=consensus[sam2$Pos[j]:length(consensus)]
                }       
                else {ref<-consensus[sam2$Pos[j]:(sam2$Pos[j]+length(reads)-1)]
                }
                m<-cbind(ref,reads)
                H[j]<-hamming.distance(m[,1],m[,2])
                #hist(H,main=paste(file.name),xlab='humming distance',breaks=50)
        }
        
        
        filename<-paste0("Output/HammingDistanceFiltering/",fname2,".pdf")
        pdf(filename, width =10, height = 5)
        par(mfrow=c(1,2))
        par(mar = c(5,4,4,2))
        plot(table(H),main=paste(fname2),xlab='Hamming distance',ylab='Counts')

        ham.distance[[i]]<-H
        names(ham.distance[i])<-paste(fname2)
        filenames[[i]]<-paste(fname2)

        Hindel<-c()
        for (k in 1:nrow(sam_indels)){
                read<-sam_indels$seq[k]
                read.s<-strsplit(read,"")
                read.s<-read.s[[1]]
                if ((sam_indels$Pos[k]+length(read.s)-1)>length(consensus)) {
                cut<-(sam_indels$Pos[k]+length(read.s)-1)-length(consensus)
                read.s<-read.s[1:(length(read.s)-cut)]
                ref=consensus[sam_indels$Pos[k]:length(consensus)]
                }        else {ref<-consensus[sam_indels$Pos[k]:(sam_indels$Pos[k]+length(read.s)-1)]
                }

                ref<-paste(ref, sep="", collapse="")
                dna<-DNAStringSet(c(paste(ref),paste(read)))
                names(dna)<-c(paste0('ref',k),paste0(file.name,".read",paste(k)))
                invisible(capture.output(align<-msa(dna)))
                aligned<-DNAStringSet(unlist(align))
                m2<-cbind(paste(aligned[1]),paste(aligned[2]))
                m2<-strsplit(m2,"")
                m2.m<-do.call(cbind,m2)

                Hindel[k]<-hamming.distance(m2.m[,1],m2.m[,2])-1
        }
        plot(table(Hindel),main=paste(fname2,"w/ indels"),xlab='Hamming distance',ylab='Counts')
        ham.distance.indel[[i]]<-Hindel
        names(ham.distance.indel[i])<-paste(fname2)


        #eliminate the reads with hamming distance>10 and save as a new datatable
        Large.ham<-(H>10)
        sam_re<-sam2[which(Large.ham==F),]
        no_indels_removed<-length(Large.ham[Large.ham==T])
        print(paste("# of removed reads w/o indels: ",no_indels_removed))

        Large.ham2<-(Hindel>10)
        sam_re2<-sam_indels[which(Large.ham2==F),]
        with_indels_removed<-length(Large.ham2[Large.ham2==T])
        print(paste("# of removed reads with indels: ",with_indels_removed))
        sam_RC<-rbind(sam_re,sam_re2)
        write.table(sam_RC, paste0("Output/sam2/", fname2,"-filtered.sam"),sep="\t", quote=F,row.names=F,col.names=F)
        
        sum[[i]]<-data.frame(fname2,mapped.reads,no_indels, with_indels,no_indels_removed,with_indels_removed)
        dev.off()
}


        
#summary of hamming distnaces without indels
hamm.sum<-do.call(rbind,ham.distance)
rown<-do.call(rbind, filenames)
rownames(hamm.sum)<-rown
hammingDist.summary<-data.frame(t(hamm.sum))
write.csv(hammingDist.summary,paste0("Output/HammingDistanceFiltering/HammingDistSummary_",Sys.Date(),".csv"))

#with indels
hamm.sum_indel<-do.call(rbind,ham.distance.indel)
rownames(hamm.sum_indel)<-rown
hammingDist.summary_indel<-data.frame(t(hamm.sum_indel))
write.csv(hammingDist.summary_indel,paste0("Output/HammingDistanceFiltering/HammingDistSummary_indels_",Sys.Date(),".csv"))

output.summary<-data.frame(do.call(rbind,sum))
colnames(output.summary)<-c("Sample ID","# mapped reads", "# reads w/o indels ", 
                            "# reads w/ indels","removed reads w/o indels", "removed w/ indels")
write.csv(output.summary,paste0("Output/HammingDistanceFiltering/Filter_Summary2_",Sys.Date(),".csv"))
