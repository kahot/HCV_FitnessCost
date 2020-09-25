#Script to create an overview files with mutation frequency data and their associated features 
#(mut freq/features based on Ref NT, Results=Overview2)

library(dplyr)
library(tidyverse)

source("Rscripts/baseRscript.R")
dir.create("Output1A/Overview2/")

#Get the file names (SeqData files)
HCVFiles_SeqData<-list.files("Output1A/SeqDataQ35/",pattern="SeqData")

#create a Overview file for each sample
Overview<-list()
for (i in 1:length(HCVFiles_SeqData)){   
        #for (i in 1:1){
        id<-substr(paste(HCVFiles_SeqData[i]),start=9,stop=15)
        print(id)
        OverviewDF<-read.csv(paste0("Output1A/SeqDataQ35/",HCVFiles_SeqData[i]),stringsAsFactors=FALSE)
        OverviewDF<-OverviewDF[,-1]

        ref<-read.dna("Data/HCVref.fasta", format = "fasta",as.character=TRUE)
        #replace ? with N
        ref<-ref[262:8800]
 
        TypeOfSite<-c()
        TypeOfSite.tv1<-c()
        TypeOfSite.tv2<-c()
        for (codon in 1:(nrow(OverviewDF)/3)) {#for each codon in the sequence
                positions <- c(codon*3-2,codon*3-1, codon*3)
                WTcodon <- OverviewDF$ref[positions]  
                if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
                        WTcodon<-c('n','n','n')
                        mutant1codon<-c('n','n','n')
                        mutant2codon<-c('n','n','n')
                        mutant3codon<-c('n','n','n')}
                else{                        
                        mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  #If the first position has transistion mutation, it's labeld as mutatnt1codon.
                        mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
                        mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
                        
                        #transversion mutation to 'a' or 'c'
                        mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
                        mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
                        #transversion mutation to 'g' or 't'
                        mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
                        mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
                }
                
                
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
                
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))
                
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))                
        } # This creates a vector showing if there is a transition/transversion mutation at a particular codon, &
        # whehter the mutation will be Syn, nonSyn, or Stop codon.
        
        OverviewDF$Type<-TypeOfSite[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv1<-TypeOfSite.tv1[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv2<-TypeOfSite.tv2[1:length(OverviewDF$pos)]
        
        Overview[[i]]<-OverviewDF[,-c(1:6)]
        names(Overview)[i]<-id   
}



###############################
#Mut rates from Geller paper 

## 1. Using 12 different Mutation Frequencies baed on from/to nucleotides
mutrates<-read.csv("Data/Geller.mutation.rates.csv", stringsAsFactors = F, row.names = 1)

Overview_summary<-list()
for (i in 1:length(Overview)){
        OverviewDF<-Overview[[i]]
        id<-substr(paste(HCVFiles_SeqData[i]),start=9,stop=15)
        
        OverviewDF$TSmutrate[OverviewDF$ref=="a"]<-mutrates$mut.rate[mutrates$mutation=="AG"]
        OverviewDF$TSmutrate[OverviewDF$ref=="c"]<-mutrates$mut.rate[mutrates$mutation=="CU"]
        OverviewDF$TSmutrate[OverviewDF$ref=="g"]<-mutrates$mut.rate[mutrates$mutation=="GA"]
        OverviewDF$TSmutrate[OverviewDF$ref=="t"]<-mutrates$mut.rate[mutrates$mutation=="UC"]

        OverviewDF$TVSmutrate.tv1[OverviewDF$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutation=="AC"])
        OverviewDF$TVSmutrate.tv1[OverviewDF$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutation=="CA"])
        OverviewDF$TVSmutrate.tv1[OverviewDF$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutation=="GC"])
        OverviewDF$TVSmutrate.tv1[OverviewDF$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutation=="UA"])
        
        OverviewDF$TVSmutrate.tv2[OverviewDF$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutation=="AU"])
        OverviewDF$TVSmutrate.tv2[OverviewDF$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutation=="CG"])
        OverviewDF$TVSmutrate.tv2[OverviewDF$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutation=="GU"])
        OverviewDF$TVSmutrate.tv2[OverviewDF$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutation=="UG"])
        
        OverviewDF$TVSmutrate.tvs[OverviewDF$ref=="a"]<-mean(mutrates$mut.rate[mutrates$mutation=="AU"],mutrates$mut.rate[mutrates$mutation=="AC"])
        OverviewDF$TVSmutrate.tvs[OverviewDF$ref=="c"]<-mean(mutrates$mut.rate[mutrates$mutation=="CG"],mutrates$mut.rate[mutrates$mutation=="CA"])
        OverviewDF$TVSmutrate.tvs[OverviewDF$ref=="g"]<-mean(mutrates$mut.rate[mutrates$mutation=="GU"],mutrates$mut.rate[mutrates$mutation=="GC"])
        OverviewDF$TVSmutrate.tvs[OverviewDF$ref=="t"]<-mean(mutrates$mut.rate[mutrates$mutation=="UG"],mutrates$mut.rate[mutrates$mutation=="UA"])
        

        for (k in 1:length(OverviewDF$pos)){
                OverviewDF$EstSelCoeff[k] <- EstimatedS(OverviewDF$TSmutrate[k],OverviewDF[k,colnames(OverviewDF)=='freq.Ts.ref'])
                OverviewDF$EstSelCoeff_transv[k] <- EstimatedS(OverviewDF$TVSmutrate.tvs[k],OverviewDF[k,colnames(OverviewDF)=='freq.transv.ref'])
                OverviewDF$EstSelCoeff_trans1[k] <- EstimatedS(OverviewDF$TVSmutrate.tv1[k],OverviewDF[k,colnames(OverviewDF)=='freq.transv1.ref'])
                OverviewDF$EstSelCoeff_trans2[k] <- EstimatedS(OverviewDF$TVSmutrate.tv2[k],OverviewDF[k,colnames(OverviewDF)=='freq.transv2.ref'])
                
                if (k%%3==1){
                        if (is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1])|is.na(OverviewDF$MajNt[k+2])) { 
                            OverviewDF$MajAA[k]<-"NA"
                            OverviewDF$WTAA[k]<-"NA"
                            OverviewDF$MUTAA[k]<-"NA"
                            OverviewDF$TVS1_AA[k]<-"NA"
                            OverviewDF$TVS2_AA[k]<-"NA"
                            }
                        else { 
                            OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k,k+1,k+2)])
                            OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$ref[c(k,k+1,k+2)])
                            OverviewDF$MUTAA[k] = seqinr::translate(c(transition(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
                            OverviewDF$TVS1_AA[k] = seqinr::translate(c(transv1(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
                            OverviewDF$TVS2_AA[k] = seqinr::translate(c(transv2(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
                            }
                } 
                if (k%%3==2){
                        if (is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1]))  {OverviewDF$MajAA[k]<-"NA"
                        OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-1,k,k+1)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$ref[c(k-1,k,k+1)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transition(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transv1(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transv2(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))}
                }
                if (k%%3==0){
                        if (is.na(OverviewDF$MajNt[k-2])|is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k]))  {  OverviewDF$MajAA[k]<-"NA"
                        OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$MajAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-2,k-1,k)])
                        OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$ref[c(k-2,k-1,k)])
                        OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transition(OverviewDF$ref[k])))
                        OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transv1(OverviewDF$ref[k])))
                        OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transv2(OverviewDF$ref[k])))}
                }
                
        }

        #Add whether AA change is drastic & makes CpG
        OverviewDF$bigAAChange<-0
        OverviewDF$bigAAChange.tv1<-0
        OverviewDF$bigAAChange.tv2<-0
        OverviewDF$makesCpG <- 0
        OverviewDF$makesCpG.tvs <- 0
        OverviewDF$makesCpG.tv1 <- 0
        OverviewDF$makesCpG.tv2 <- 0
        
        for(j in 2:nrow(OverviewDF)-1){
                WT <- amCat(OverviewDF[j,'WTAA'])
                MUT <- amCat(OverviewDF[j,'MUTAA'])
                MUT1<-amCat(OverviewDF[j,'TVS1_AA'])
                MUT2<-amCat(OverviewDF[j,'TVS2_AA'])
                
                if (WT != MUT) OverviewDF$bigAAChange[j] <- 1
                if (WT != MUT1) OverviewDF$bigAAChange.tv1[j] <- 1
                if (WT != MUT2) OverviewDF$bigAAChange.tv2[j] <- 1
                
                trip <- OverviewDF$ref[c(j-1, j,j+1)]
                if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                        next
                else{
                        if (trip[1] == "c" & trip[2] == "a" ) OverviewDF$makesCpG[j] <- 1 
                        if (trip[2] == "t" & trip[3] == "g")  OverviewDF$makesCpG[j] <- 1
                        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tvs[j] <- 1
                        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tvs[j] <- 1
                        
                        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tv2[j] <- 1                                
                        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tv1[j] <- 1
                        }
        }
        
                
        write.csv(OverviewDF,paste0("Output1A/Overview2/",id,"overview2.csv"))
        
        Overview_summary[[i]]<-OverviewDF
        print(id)
}        


