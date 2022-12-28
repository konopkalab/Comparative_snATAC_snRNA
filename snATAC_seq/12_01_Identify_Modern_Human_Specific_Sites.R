rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
library(Seurat)
library(dplyr)
library(tidyverse)
library(rio)
library(GenomicRanges)
library(liftOver)
library(data.table)
library('bedr')
source("utility_functions.R")

# Load old variants (hg19). Allele frequencies are updated with 1000G Phase 3.
hvars = fread('HumanDerived_SNC_bothgq30.all_combined_maxsco_ranked_UPDATED_PHASE3_1000G.tsv', header = T) %>% as.data.frame

# Loop over chromosomes to keep the variants common in all neanderthals
chrs = c(paste0('chr', 1:22), 'chrX')
cmnNeanL = list()
for(i in 1:length(chrs)){

	# Keep only if ancestral in neanderthals
	hvarsSub = hvars[hvars[,1] == gsub('chr', '', chrs[i]),]
	tmp = hvarsSub[,7] %>% gsub('/', ',', .) %>% strsplit(., ',') %>% do.call(rbind, .)
	tmp = cbind(hvarsSub[,5], tmp[,1:2])
	cons = apply(tmp, 1, function(x){ (x[2] == x[3]) & (x[1] != x[2]) }) %>% which
	hvarsSubNean = hvarsSub[cons,]
	hvarsSubNean$neanAllele = tmp[cons, 2] 

	# Load Other Neanderthals
	vind = fread(paste0('Vindija/', chrs[i], '_mq25_mapab100.vcf.gz'))
	char = fread(paste0('Chagyrskaya/', chrs[i], '.noRB.vcf.gz'))

	# Keep the same positions identified in all individuals
	cmnPos = Reduce(intersect, list(hvarsSubNean$Pos, vind$POS, char$POS))
	hvarsSubNeanCMN = hvarsSubNean[hvarsSubNean$Pos %in% cmnPos,]
	vindVar = vind[vind$POS %in% cmnPos, ]
	charVar = char[char$POS %in% cmnPos, ]

	# Match the order of the positions
	hvarsSubNeanCMN = hvarsSubNeanCMN[match(cmnPos, hvarsSubNeanCMN$Pos), ]
	vindVar = vindVar[match(cmnPos, vindVar$POS), ]
	charVar = charVar[match(cmnPos, charVar$POS), ]

	# Keep the positions that are the same across all three neanderthals
	cmnNeanL[[i]] = hvarsSubNeanCMN[hvarsSubNeanCMN$neanAllele == vindVar[,5] & hvarsSubNeanCMN$neanAllele == charVar[,5],]

	print(i)
}


allComb = do.call(rbind, cmnNeanL)
saveRDS(allComb, 'CommonNeanderthal_N3_Prufer2014_Table.RDS')

####
## Also restrict by the Denisovan
####

# Read neanderthals
allCombN = readRDS('CommonNeanderthal_N3_Prufer2014_Table.RDS')

# Keep only if consistent between neanderthal & denisovans
tmp = hvars[,7] %>% gsub('/', ',', .) %>% strsplit(., ',') %>% do.call(rbind, .)
cons = apply(tmp, 1, function(x){length(table(x)) == 1}) %>% which
hvarsCons = hvars[cons,]

# Keep only the ones >90% frequency
hvarsCons = hvarsCons[hvarsCons[, 'Human_Maj_1000G_PHASE3'] > 0.9,]
allCombN$id = paste0(allCombN[,1], ':', allCombN[,2])
hvarsCons$id = paste0(hvarsCons[,1], ':', hvarsCons[,2])
allCombND = allCombN[allCombN$id %in% hvarsCons$id,]

saveRDS(allCombND, 'CommonNeanderthal_N3_D1_Prufer2014_Table.RDS')
















