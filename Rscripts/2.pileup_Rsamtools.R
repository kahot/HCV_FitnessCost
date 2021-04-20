# Read BAM files in and convert them to frequency tables using Rsamtools

library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#List of bam files to be processed:
bamfiles<-list.files("Output/bam2/",pattern="bam$")

#read each bam file, convert it to a freq. table using 'pileup' of Rsamtools, and save as a csv file. 
for (i in 1:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0(paste0("Output/bam2/",bam),'.bai')
        bf<-BamFile(paste0("Output/bam2/",bam), index=index)
        
        file.name<-paste(bam)
        file.name<-substr(file.name,start=1,stop=10 )
        p_param <- PileupParam(max_depth=60000,distinguish_strands=FALSE,include_insertions=TRUE)
        result<-pileup(bf, pileupParam = p_param)
        summary<-pileupFreq(result)
        
        print(file.name)
        write.csv(summary, file=paste0("Output/CSV/",file.name,".csv",collapse=""))

}

