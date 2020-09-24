library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#number of sampels to process
bamfiles<-list.files("Output1A/bam2/",pattern="bam$")
for (i in 1:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0(paste0("Output1A/bam2/",bam),'.bai')
        bf<-BamFile(paste0("Output1A/bam2/",bam), index=index)
        
        file.name<-paste(bam)
        file.name<-substr(file.name,start=1,stop=10 )
        p_param <- PileupParam(max_depth=60000,distinguish_strands=FALSE,include_insertions=TRUE)
        result<-pileup(bf, pileupParam = p_param)
        summary<-pileupFreq(result)
        
        print(file.name)
        write.csv(summary, file=paste0("Output1A/CSV/",file.name,".csv",collapse=""))

}

