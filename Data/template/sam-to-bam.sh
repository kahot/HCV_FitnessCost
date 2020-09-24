#!/bin/bash

samtools view -S -b *filtered.sam > *filtered.bam

#6. sort the bam file
samtools sort *filtered.bam -o *.sort.bam

#7. index the bam file
samtools index *.sort.bam *.sort.bam.bai






#5. convert sam to bam
samtools view -S -b D75335/D75335_BWAmapped30.sam > D75335/D75335_BWAmapped30.bam


#6. sort the bam file
samtools sort D75335/D75335_BWAmapped30.bam -o D75335/D75335_BWAmapped30.sorted.bam

#index the bam file
samtools index D75335/D75335_BWAmapped30.sorted.bam D75335/D75335_BWAmapped30.sorted.bam.bai

