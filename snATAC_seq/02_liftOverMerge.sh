#!/bin/bash

# Load modules
module load bedtools/2.29.0

# Set directories
dir=merge_liftover
cd $dir
chainfile1=chain_files/panTro5ToHg38.over.chain.gz
chainfile2=chain_files/hg38ToPanTro5.over.chain.gz
chainfile3=chain_files/rheMac10ToHg38.over.chain.gz
chainfile4=chain_files/hg38ToRheMac10.over.chain.gz

# Get rid of mitochondrial peaks
grep -v 'chrM' human_consensus.bed > tmp.bed
mv tmp.bed human_consensus_peaks.bed
grep -v 'chrM' chimp_consensus.bed > tmp.bed
mv tmp.bed chimp_consensus.bed
grep -v 'chrM' macaque_consensus.bed > tmp.bed
mv tmp.bed macaque_consensus.bed

# This function performs liftover while keeping original peaks in the file.
# It removes multiple hits and peaks that have more than 2-fold difference between original and lifted.
liftfunc(){

# Liftover
~/programs/liftOver $1 $2 -multiple -minMatch=0.5 -noSerial $3 "unlifted_"$3

# Get rid of multiples
sort -k4,4 $3 | uniq -u -f 3 > "tmp_"$3

# Add insert sizes
sed -i -e "s/:/\t/g" -e "s/-/\t/g" "tmp_"$3
awk 'BEGIN {OFS="\t"} {$7=$3-$2} 1' "tmp_"$3 > tmp1.bed
awk 'BEGIN {OFS="\t"} {$8=$6-$5} 1' tmp1.bed > tmp2.bed

# Remove lines with more than two-fold difference in lifted regions
awk '($8/$7 < 2) && ($8/$7 > 0.5)' tmp2.bed | cut -f1-3 > $3
awk '($8/$7 < 2) && ($8/$7 > 0.5)' tmp2.bed | cut -f4-6 > "INCOORDS_"$3


rm tmp1.bed tmp2.bed "tmp_"$3
}

# Liftover chimp and macaque to human. Lifted files are named as 'initialSpecies_liftedSpeciesL1.bed'
liftfunc chimp_consensus.bed $chainfile1 chimp_humanL1.bed
liftfunc macaque_consensus.bed $chainfile3 macaque_humanL1.bed

# Merge chimp, macaque and human if at least 200bp close to each other.
# This is because we had macs2 use 200bp window for each cutsite.
cat chimp_humanL1.bed macaque_humanL1.bed human_consensus.bed | sort -k1,1 -k2,2n > concat.bed
bedtools merge -d 200 -i concat.bed > merged.bed

# Filter any peaks that do not map to nuclear chromosomes in human
grep -v 'Un\|random\|chrM' merged.bed > tmp.bed
mv tmp.bed merged.bed

# Reciprocal liftover of merged peaks to chimp
liftfunc merged.bed $chainfile2 merged_chimpL1.bed
liftfunc merged_chimpL1.bed $chainfile1 chimpL1_mergedL1.bed

# Add chimp coordinates as last column
sed -i 's/\t/-/g' INCOORDS_chimpL1_mergedL1.bed
paste chimpL1_mergedL1.bed INCOORDS_chimpL1_mergedL1.bed > chimpL1_mergedL1_coord.bed

# Get rid of peaks not mapping main chromosomes
grep -v 'NW' chimpL1_mergedL1_coord.bed > tmp.bed
mv tmp.bed chimpL1_mergedL1_coord.bed

# Reciprocal liftover of merged peaks to macaque
liftfunc merged.bed $chainfile4 merged_macaqueL1.bed
liftfunc merged_macaqueL1.bed $chainfile3 macaqueL1_mergedL1.bed

# Add macaque coordinates as last column
sed -i 's/\t/-/g' INCOORDS_macaqueL1_mergedL1.bed
paste macaqueL1_mergedL1.bed INCOORDS_macaqueL1_mergedL1.bed > macaqueL1_mergedL1_coord.bed

# Get rid of peaks not mapping main chromosomes
grep -v 'NW' macaqueL1_mergedL1_coord.bed > tmp.bed
mv tmp.bed macaqueL1_mergedL1_coord.bed

# Intersect initial peaks with reciprocally lifted peaks to both chimp and macaque.
# Keep only if overlap is > 50% in both initial and lifted peaks.
bedtools intersect -a merged.bed \
		   -b macaqueL1_mergedL1_coord.bed -f 0.5 -F 0.5 -wa -wb | cut -f1-3,7 | awk '!seen[$1,$2,$3]++' > merged_macfilt.bed

bedtools intersect -a merged.bed \
		   -b chimpL1_mergedL1_coord.bed -f 0.5 -F 0.5 -wa -wb | cut -f1-3,7 | awk '!seen[$1,$2,$3]++' > merged_chimpfilt.bed

cut -f1-3 merged_chimpfilt.bed > tmpchimp.bed
cut -f1-3 merged_macfilt.bed > tmpmac.bed

bedtools intersect -a tmpchimp.bed \
		   -b tmpmac.bed -f 0.5 -F 0.5 -wa -wb | awk '!seen[$1,$2,$3]++' | cut -f1-3 > merged_final.bed

# Filter out the blacklisted regions of human genome
bedtools intersect -v -a merged_final.bed \
			-b GENOME_REGIONS/ENCFF419RSJ_blacklistGrCH38.bed > tmp.bed

mv tmp.bed merged_final.bed

# Grab species coordinates of final merged file.
sed 's/\t/-/g' merged_final.bed > finalcoord.bed
cut -f1-3 merged_chimpfilt.bed | sed 's/\t/-/g' > tmp.bed
paste tmp.bed merged_chimpfilt.bed > chimpcoord.bed
awk 'FNR==NR{a[$1];next}($1 in a){print}' finalcoord.bed chimpcoord.bed > merged_final_chimp.bed

cut -f1-3 merged_macfilt.bed | sed 's/\t/-/g' > tmp.bed
paste tmp.bed merged_macfilt.bed > macaquecoord.bed
awk 'FNR==NR{a[$1];next}($1 in a){print}' finalcoord.bed macaquecoord.bed > merged_final_macaque.bed


# END OF SCRIPT

