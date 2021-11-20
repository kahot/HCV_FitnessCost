#!/bin/bash

head -n 3 ~/programs/HCV/Output/sam/D75000_me_clipped.sam > h3
cat h3 ~/programs/HCV/Output/sam2/D75000_me-filtered.sam > ~/programs/HCV/Output/sam2/D75000_mefiltered2.sam

#convert to bam file
samtools view -S -b ~/programs/HCV/Output/sam2/D75000_mefiltered2.sam > ~/programs/HCV/Output/bam2/D75000_me_filtered.bam

#sort bam file
samtools sort ~/programs/HCV/Output/bam2/D75000_me_filtered.bam -o ~/programs/HCV/Output/bam2/D75000_me.sort.bam

#index the bam file
samtools index ~/programs/HCV/Output/bam2/D75000_me.sort.bam ~/programs/HCV/Output/bam2/D75000_me.sort.bam.bai

rm ~/programs/HCV/Output/sam2/D75000_mefiltered2.sam ~/programs/HCV/Output/bam2/D75000_me_filtered.bam


head -n 3 ~/programs/HCV/Output/sam/D75000_un_clipped.sam > h3
cat h3 ~/programs/HCV/Output/sam2/D75000_un-filtered.sam > ~/programs/HCV/Output/sam2/D75000_unfiltered2.sam

#convert to bam file
samtools view -S -b ~/programs/HCV/Output/sam2/D75000_unfiltered2.sam > ~/programs/HCV/Output/bam2/D75000_un_filtered.bam

#sort bam file
samtools sort ~/programs/HCV/Output/bam2/D75000_un_filtered.bam -o ~/programs/HCV/Output/bam2/D75000_un.sort.bam

#index the bam file
samtools index ~/programs/HCV/Output/bam2/D75000_un.sort.bam ~/programs/HCV/Output/bam2/D75000_un.sort.bam.bai

rm ~/programs/HCV/Output/sam2/D75000_unfiltered2.sam ~/programs/HCV/Output/bam2/D75000_un_filtered.bam
