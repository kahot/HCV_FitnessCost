#!/bin/bash


#1  Align the filtered reads using bwa to the consensus sequence. 

bwa index -p D75000 -a is ~/programs/HCV/Output/Consensus/D75000_consensus.fasta

bwa mem -t 8 -k 15 D75000 ~/programs/HCV/Output/merged/D75000_merged.fq  > ~/programs/HCV/Output/Bash/D75000/D75000_ConMapped.merge.sam
bwa mem -t 8 -k 15 D75000 ~/programs/HCV/Output/unmerged/D75000_unmerged.fq  > ~/programs/HCV/Output/Bash/D75000/D75000_ConMapped.un.sam


#2 Hard clipped the aligned reads
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output/Bash/D75000/D75000_ConMapped.merge.sam > ~/programs/HCV/Output/sam/D75000_me_clipped.sam
java -jar ~/programs/jvarkit/dist/biostar84452.jar ~/programs/HCV/Output/Bash/D75000/D75000_ConMapped.un.sam > ~/programs/HCV/Output/sam/D75000_un_clipped.sam

rm ~/programs/HCV/Output/Bash/D75000/D75000_ConMapped.merge.sam ~/programs/HCV/Output/Bash/D75000/D75000_ConMapped.un.sam

