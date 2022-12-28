rm(list = ls())
require("rphast")
library(ape)
library(dplyr)
library(parallel)
library(Biostrings)

####
## PREPARE DATA
####

# Take each chromosome as argument from bash script
args = commandArgs(trailingOnly = TRUE)
chr = as.character(args[1])
print(chr)

# Load MSA
align = readRDS(paste0('HAR_30WAY/MSA/', chr, '_anthropoids_msa.RDS'))

# Read MSA offset
ofsDF = read.table(paste0("UCSC_MULTIZ30_2017/head_", chr, ".maf"), fill = T)
ofs = as.numeric(ofsDF[2,3])

# Load tree topology
treetop = read.table("UCSC_MULTIZ30_2017/treeforR.txt")
treetop = treetop$V1

# Keep only great apes
tree = read.tree(text = treetop)
dropspecies1 = c('tarSyr2', 'micMur3', 'proCoq1', 'eulMac1', 'eulFla1', 'otoGar3', 'mm10', 'canFam3', 'dasNov3')
dropspecies2 = c('papAnu3', 'rhiBie1', 'rhiRox1', 'ponAbe2', 'nasLar1', 'panPan2')
dropspecies3 = c('rheMac8', 'macFas5', 'macNem1', 'cerAty1', 'chlSab2', 'manLeu1', 'colAng1', 'calJac3', 'saiBol1', 'cebCap1', 'aotNan1')

dropspecies = Reduce(union, list(dropspecies1, dropspecies2, dropspecies3))

tree = drop.tip(tree, dropspecies)
treetop = write.tree(tree)
align2 = align[tree$tip,]

# Load CREs
pks = read.feat('merged_lifted_human_sorted.bed')
pks = pks[pks[,1] == chr,]

pks$seqname = "hg38"

# Adjust CREs by offset
pks2 = pks
pks2[,4] = pks2[,4] - ofs
pks2[,5] = pks2[,5] - ofs

####
## EXTRACT CREs FROM ALIGNMENT
####

# Loop through peaks and extract from the multi-species alignment
seq1 = seq(1, nrow(pks2), 250)
seq2 = c(seq1[2:length(seq1)] - 1, nrow(pks2))
qL = list()
for(j in 1:length(seq1)){

	qL[[j]] = mclapply(seq1[j]:seq2[j], mc.cores = 23, function(x){
				st = pks2[x,4]
				end = pks2[x,5]
				sub.msa(align2, start.col = st, end.col = end, pointer.only = F)})
	print(j)
}

saveRDS(unlist(qL, recursive = F), paste0(chr, '_extractedSequences.RDS'))


####
## FIND HS SUBSTITUTIONS
####

# Loop through the extracted sequences and find human-specific SNCs
seq1 = seq(1, nrow(pks2), 250)
seq2 = c(seq1[2:length(seq1)] - 1, nrow(pks2))
combdfL2 = list()
for(j in 1:length(seq1)){

	qL = extSeq[seq1[j]:seq2[j]]

		combdfL = list()
		for(i in 1:length(qL)){

			# Indice of the peak
			inds = seq1[j]:seq2[j]
			ind = inds[i]
			
			# Position of the peak
			st = pks[ind,4]
			end = pks[ind,5]

			# Extract the alignment for each CRE
			q = qL[[i]]
			write.msa(q, format = 'FASTA', file = paste0('tmp/', 'tmp', st, '.fasta'))
			fastaFile = readBStringSet(paste0('tmp/', 'tmp', st, '.fasta'))
			seq_name = names(fastaFile)
			sequence = paste(fastaFile)

			# Alignment as data frame
			df = lapply(sequence, function(x){strsplit(x, '') %>% unlist}) %>% do.call(rbind, .) %>% as.data.frame
			rownames(df) = seq_name
			
			# Keep only human-specific nucleotide changes
			tmp1 = apply(df, 2, function(x){(sum(x[1] != x[2:4]) == 3) &
								(sum(grepl('A|C|T|G', x[2:4])) == 3)})

			tmp2 = apply(df, 2, function(x){length(table(x[2:4])) == 1})
			df2 = rbind(tmp1, tmp2)
			sncInd = which(colSums(df2) == 2)

			# Move to the next one if there are no SNCs
			if(length(sncInd) == 0){next}

			# Create data frame of SNCs with position
			sncInd = as.numeric(sncInd)
			ucscPos = sapply(sncInd, function(x){st+x-1})

			combdf = data.frame(chr = chr, pos = ucscPos,
						human = as.character(unique(df[1,sncInd])),
						ancestral = as.character(unique(df[2:4,sncInd])),
						CRE = paste0(chr, ':', pks[ind,4], '-', pks[ind,5]))

			combdfL[[i]] = combdf

			if(i%%50==0){print(i)}
		}
		
		combdfL2[[j]] = do.call(rbind, combdfL)
		print(paste0('CRE# ', seq2[j]))
}

allSNCs = do.call(rbind, combdfL2)
saveRDS(allSNCs, paste0(chr, '_SNCs.RDS'))


