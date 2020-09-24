#Step 1. Create bash executable files to process (filter, trim and map) the fastq fies
library(stringr)
library(R.utils)


#create the bash files to run bbmap and bwa
# read the template command text file:
cmmd<-readLines("Data/template/D7xxxx1_Q35.sh")
cmmd2<-readLines("Data/template/D7xxxx2.sh")
cmmd3<-readLines("Data/template/D7xxxx3.sh")

cm<-readLines("Data/template/D7xxxx1_Q35_2.sh")

#choose the fastq files to be prrocessed
fq<-list.files("Data/unzipped/",pattern="fastq") 

#create vector of odd numbers:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
for (i in 1:length(fq2)){
        #choose the paired reads fastq files
        fa1<-fq2[i]
        fa2<-gsub(pattern="R1",replace="R2",x=fa1)
        fname<-substr(fa1,start=1,stop=7)
        new<-gsub(pattern="D75000-HCV_S11_L001_R1_001.fastq", replace=paste0(fa1),x=cm)
        new<-gsub(pattern="D75000-HCV_S11_L001_R2_001.fastq", replace=paste0(fa2),x=new)
        new<-gsub(pattern="D75000",replace=paste0(fname),x=new)
        writeLines(new, con=paste0("Bashscripts/Bash1_Q35/",fname,".sh"))
        
        new2<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd2)
        writeLines(new2, con=paste0("Bashscripts/bash2_Q35/",fname,"2.sh"))
    
        new3<-gsub(pattern="D75000",replace=paste0(fname),x=cmmd3)
        writeLines(new3, con=paste0("Bashscripts/bash3_Q35/",fname,"3.sh"))
        
        
}
