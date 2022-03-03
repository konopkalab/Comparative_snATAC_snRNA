#!/bin/bash

# Load modules
module load samtools
module load bedtools

dir=consensus_peaks

cd $dir
mkdir overlap_table

# Sort reads
sort -k1,1 -k2,2n merged_final.bed > merged_lifted_human_sorted.bed

# Overlap reads and peaks
bedtools intersect -sorted -wo -a human1_reads.bed -b merged_lifted_human_sorted.bed | cut -f5-8 > human1_peak_ov.bed

bedtools intersect -sorted -wo -a human2_reads.bed -b merged_lifted_human_sorted.bed | cut -f5-8 > human2_peak_ov.bed

bedtools intersect -sorted -wo -a human3_reads.bed -b merged_lifted_human_sorted.bed | cut -f5-8 > human3_peak_ov.bed

bedtools intersect -sorted -wo -a human4_reads.bed -b merged_lifted_human_sorted.bed | cut -f5-8 > human4_peak_ov.bed






