#Merge two freq. tables (mapped and unmapped CSVs) created in #2 for each sample,
#and convert the frquency counts into mutation freqeunceis ('SeqFile')

library(tidyverse)
library(plyr)
source("Rscripts/baseRscript.R")

#List of csv files
HCVFiles<-list.files("Output/CSV/",pattern="un.csv")
MergeFiles<-list.files("Output/CSV/", pattern="me.csv")

start<-264
end<-8800
no<-data.frame("pos"=c(start:end))

for (i in 1:length(MergeFiles)){
        print(i)
        id<-substr(paste(HCVFiles[i]),start=1,stop=7)
        print(id)
        merge<-read.csv(paste("Output/CSV/",MergeFiles[i],sep=""))
        merge<-merge[,-c(1,2,8,9)]
        colnames(merge)<-c("pos","mA","mC","mG","mT", "mDel","mIns")
        unmerge<-read.csv(paste("Output/CSV/",HCVFiles[i],sep=""))
        unmerge<-unmerge[,-c(1,2,8,9)]
        colnames(unmerge)<-c("pos","uA","uC","uG","uT", "uDel","uIns")

        SeqData<-join(merge,unmerge,by="pos")
        SeqData[is.na(SeqData)] <- 0
        SeqData$A<-SeqData$mA+SeqData$uA
        SeqData$C<-SeqData$mC+SeqData$uC
        SeqData$G<-SeqData$mG+SeqData$uG
        SeqData$t<-SeqData$mT+SeqData$uT
        SeqData$Del<-SeqData$mDel+SeqData$uDel
        SeqData$Ins<-SeqData$mIns+SeqData$uIns
        SeqData<-SeqData[,c("pos","A","C","G","t","Del","Ins")]
        SeqData$TotalReads<-rowSums(SeqData[2:5],na.rm=T)
        SeqData$TotalReads_indels<-rowSums(SeqData[2:7],na.rm=T)
        colnames(SeqData)[2:7]<-c("a","c","g","t","deletion","insertion")
        
        #determine the majority nucleotide base at each site
        SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
        
        #read the H77 refrence sequence:
        reference<-read.dna("Data/HCVref.fasta", format = "fasta",as.character=TRUE)
        ref.code<-reference[coding.start:coding.end]
        SeqData<-merge(no,SeqData,by="pos",all.x=T)
        SeqData$ref<-ref.code[1:length(SeqData[,1])]
        
        #check that the right position is read  in the right reading frame
        print(seqinr::translate(SeqData$MajNt[79:110]))
        print(seqinr::translate(SeqData$MajNt[271:300]))
        
        SeqData$transition.maj<-NA
        SeqData$transition.ref<-NA
        for (j in 1:nrow(SeqData)) SeqData$transition.maj[j]<-transition(SeqData$MajNt[j])        
        for (j in 1:nrow(SeqData)) SeqData$transition.ref[j]<-transition(SeqData$ref[j])
        
        #rearrange the columns
        SeqData<-SeqData[,c("a","c","g","t","deletion","insertion","pos","TotalReads","TotalReads_indels","MajNt","ref","transition.maj","transition.ref")]
        
        #determine Transition mutation freq of every site.
        for (k in 1:nrow(SeqData)){
                if (is.na(SeqData$MajNt[k])) {
                        SeqData$freq.Ts[k]<-NA #transition mutations
                        SeqData$freq.Ts.ref[k]<-NA
                        
                        SeqData$freq.transv[k]<-NA #transversion mutations
                        SeqData$freq.transv.ref[k]<-NA
                        SeqData$freq.transv1[k]<-NA
                        SeqData$freq.transv2[k]<-NA
                        SeqData$freq.transv1.ref[k]<-NA
                        SeqData$freq.transv2.ref[k]<-NA
                        
                        SeqData$freq.mutations.ref[k]<-NA #all mutations
                        SeqData$freq.mutations[k]<-NA
                        SeqData$freq.mutations.indels[k]<-NA
                        SeqData$freq.mutations.indels.ref[k]<-NA
                }
                else {
                MajNum <- SeqData [k,paste0(SeqData$MajNt[k])]
                MutNum1<- SeqData [k,paste0(SeqData$transition.maj[k])]
                WTNum <- SeqData [k,paste0(SeqData$ref[k])]
                MutNum2<- SeqData [k,paste0(SeqData$transition.ref[k])]
                
                SeqData$freq.Ts[k]<-MutNum1/SeqData$TotalReads[k]
                SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
                
                
                #mutation frequencies of all transversion mutataions
                if (SeqData$MajNt[k]=="a"|SeqData$MajNt[k]=='g'){
                        TrvMutNum<-SeqData[k,"c"]+SeqData[k,"t"]}
                if (SeqData$MajNt[k]=="c"|SeqData$MajNt[k]=="t"){
                        TrvMutNum<-SeqData[k,"a"]+SeqData[k,"g"]}
                SeqData$freq.transv[k]<-TrvMutNum/SeqData$TotalReads[k]
                if (SeqData$ref[k]=="a"|SeqData$ref[k]=='g'){
                        TrvMutNum2<-SeqData[k,"c"]+SeqData[k,"t"]}
                if (SeqData$ref[k]=="c"|SeqData$ref[k]=="t"){
                        TrvMutNum2<-SeqData[k,"a"]+SeqData[k,"g"]}
                SeqData$freq.transv.ref[k]<-TrvMutNum2/SeqData$TotalReads[k]
                
                #Frequenceis for specific transversion mutations (1 & 2)
                Tvs1Num<-SeqData[k,paste0(transv1(SeqData$MajNt[k]))]
                Tvs2Num<-SeqData[k,paste0(transv2(SeqData$MajNt[k]))]
                SeqData$freq.transv1[k]<-Tvs1Num/SeqData$TotalReads[k]
                SeqData$freq.transv2[k]<-Tvs2Num/SeqData$TotalReads[k]
                Tvs1rNum<-SeqData[k,paste0(transv1(SeqData$ref[k]))]
                Tvs2rNum<-SeqData[k,paste0(transv2(SeqData$ref[k]))]
                SeqData$freq.transv1.ref[k]<-Tvs1rNum/SeqData$TotalReads[k]
                SeqData$freq.transv2.ref[k]<-Tvs2rNum/SeqData$TotalReads[k]
                
                
                #Frequencies of all SNPs (no indels)
                AllMutNum<-SeqData$TotalReads[k]-MajNum
                AllMutNum2<-SeqData$TotalReads[k]-WTNum
                
                SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
                SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
                
                #with indels
                SeqData$freq.mutations.indels[k]<-AllMutNum/SeqData$TotalReads_indels[k]
                SeqData$freq.mutations.indels.ref[k]<-AllMutNum2/SeqData$TotalReads_indels[k]
 
                }
        }
        
        write.csv(SeqData,paste0("Output/SeqDataQ35/SeqData_",id,".csv"))
}


