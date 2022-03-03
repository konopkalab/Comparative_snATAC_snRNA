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
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

####
## MODERN VARIANT DATASET (Prufer, 2013) https://www.nature.com/articles/nature12886
####
hvars = fread('~/workdir/reference_genomes/ARCHAIC_HUMAN/MODERN_OVER_ARHAIC/HumanDerived_SNC_bothgq30.all_combined_maxsco_ranked_UPDATED_PHASE3_1000G.tsv', header = T) %>% as.data.frame

# Keep only the ones >90% frequency
hvars = hvars[hvars[, 'Human_Maj_1000G_PHASE3'] > 0.9,]

hvarsbed = hvars[, c(1:2)]
hvarsbed$end = hvarsbed[,2] + 1
colnames(hvarsbed) = c('chr', 'start', 'end')
hvarsbed$chr = gsub('^', 'chr', hvarsbed$chr)

hvarsbed = makeGRangesFromDataFrame(hvarsbed)

# LiftOver to hg38
hg19hg38 = import.chain('hg19ToHg38.over.chain')
hvars_hg38 = liftOver(hvarsbed, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))

# Sort
hvars_hg38 = bedr.sort.region(hvars_hg38, verbose = F)
saveRDS(hvars_hg38, 'HumanDerived_SNC_bothgq30.all_combined_maxsco_ranked_HG38_MIN09.RDS')

####
## ENRICHMENT BY REGRESSION
####

# Load species-specific CREs and background per cell-type
dars = fread('COMBINED_Signpeaks_DARs.txt', header = T)
dars = as.data.frame(dars)
bcg = fread('COMBINED_Allpeaks_DARs.txt', header = T)
bcg = as.data.frame(bcg)

bcg_bed = bcg$peak %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))
bcg = cbind(bcg, bcg_bed)


# ALL CELL TYPES #

# Get CREs
bcgAll = bcg[,c('V1','V2','V3')] %>% add_column(peak = paste0(.[,1], ':', .[,2], '-', .[,3]))
bcgAll = bcgAll[!duplicated(bcgAll$peak),]
hsAll = bcg[bcg$Evolution == 'Human_Specific', 'peak']
csAll = bcg[bcg$Evolution == 'Chimp_Specific', 'peak']
bcgAll = bedr.sort.region(bcgAll, verbose = F)

# Overlap CREs with modern variants
bcg_bedr <- bedr(input = list(a = bcgAll, b = hvars_hg38), method = "intersect", params = "-loj", verbose = F)
bcg_bedr$modern_ov = ifelse(bcg_bedr$start.b > -1, 1, 0)
bcg_bedr$length = bcg_bedr[,3] - bcg_bedr[,2]
bcg_bedr$length = log2(bcg_bedr$length)

saveRDS(bcg_bedr, 'AllOverlap.RDS')

# Sum number of modern variant overlap per CRE
bcgAgg = aggregate(modern_ov ~ peak, bcg_bedr, FUN = sum)
bcgAgg$length = bcg_bedr[match(bcgAgg$peak, bcg_bedr$peak), 'length']

bcgAgg$is_hs = ifelse(bcgAgg$peak %in% hsAll, 'HS', 'NS')
bcgAgg$is_cs = ifelse(bcgAgg$peak %in% csAll, 'CS', 'NS')
bcgAgg$is_hs = factor(bcgAgg$is_hs, levels = c('NS', 'HS'))
bcgAgg$is_cs = factor(bcgAgg$is_cs, levels = c('NS', 'CS'))

humanAll <- glm.nb(modern_ov ~ length + is_hs, data = bcgAgg)
chimpAll <- glm.nb(modern_ov ~ length + is_cs, data = bcgAgg)


# Check enrichment after accounting for the length of CREs
ctypes = names(table(dars$cluster))
human_h1 = list()
chimp_c1 = list()
bcgSaves = list()
bcg_bedrL = list()
require(pscl)
require(MASS)
require(boot)
for(i in 1:length(ctypes)){

	# Get CREs
	bcg_ct = bcg[bcg$cluster == ctypes[i], c('V1','V2','V3')] %>% add_column(peak = paste0(.[,1], ':', .[,2], '-', .[,3]))
	hs_ct = bcg[bcg$cluster == ctypes[i] & bcg$Evolution == 'Human_Specific', 'peak']
	cs_ct = bcg[bcg$cluster == ctypes[i] & bcg$Evolution == 'Chimp_Specific', 'peak']

	bcgAgg_ct = bcgAgg[bcgAgg$peak %in% bcg_ct$peak,]
	bcgAgg_ct$is_hs = ifelse(bcgAgg_ct$peak %in% hs_ct, 'HS', 'NS')
	bcgAgg_ct$is_cs = ifelse(bcgAgg_ct$peak %in% cs_ct, 'CS', 'NS')
	bcgAgg_ct$is_hs = factor(bcgAgg_ct$is_hs, levels = c('NS', 'HS'))
	bcgAgg_ct$is_cs = factor(bcgAgg_ct$is_cs, levels = c('NS', 'CS'))

	human_h1[[i]] <- glm.nb(modern_ov ~ length + is_hs, data = bcgAgg_ct)
	chimp_c1[[i]] <- glm.nb(modern_ov ~ length + is_cs, data = bcgAgg_ct)

	# For Table
	bcgSave = bcg_bedr[bcg_bedr$modern_ov > 0,]
	bcgSave = dplyr::select(bcgSave, -c('width', 'modern_ov', 'strand'))
	bcgSave$is_hs = ifelse(bcgSave$peak %in% hs_ct, 'HS', 'NS')
	bcgSave$is_cs = ifelse(bcgSave$peak %in% cs_ct, 'CS', 'NS')
	colnames(bcgSave) = c('CRE_Chromosome', 'CRE_Start', 'CRE_End', 'ModernVariant_Chromosome', 'ModernVariant_Start', 'ModernVariant_End', 'Length', 'CRE', 'is_hs', 'is_cs')
	bcgSave$CellType = ctypes[i]
	bcgSaves[[i]] = bcgSave

	print(i)
}

human_h1[[length(ctypes) + 1]] = humanAll
chimp_c1[[length(ctypes) + 1]] = chimpAll

# Save overlap as database
bcgHvarOv = do.call(rbind, bcgSaves)
rio::export(bcgHvarOv, 'Peak_ModernVar_Overlap.xlsx')

# Save and plot enrichments
saveRDS(human_h1, 'Human_Peaks_NegBinom_LengthLog2.RDS')
saveRDS(chimp_c1, 'Chimp_Peaks_NegBinom_LengthLog2.RDS')

human_pval = sapply(human_h1, function(x){coef(summary(x))[3,4]})
chimp_pval = sapply(chimp_c1, function(x){coef(summary(x))[3,4]})

OR_human = sapply(human_h1, function(x){coef(summary(x))[3,1]})
OR_chimp = sapply(chimp_c1, function(x){coef(summary(x))[3,1]})

presdf = data.frame(human_pval = human_pval, chimp_pval = chimp_pval, Coef_Human = OR_human, Coef_Chimp = OR_chimp, CellType = c(ctypes, 'All'))
presdf$log10FDR_Human = -log10(p.adjust(presdf$human_pval, method = 'BH'))
presdf$log10FDR_Chimp = -log10(p.adjust(presdf$chimp_pval, method = 'BH'))
presdf[is.na(presdf$log10FDR_Human) | is.nan(presdf$log10FDR_Human), 'log10FDR_Human'] = 0
presdf[is.na(presdf$log10FDR_Chimp) | is.nan(presdf$log10FDR_Chimp), 'log10FDR_Chimp'] = 0

Coefmelt = melt(measure.vars = c('Coef_Human', 'Coef_Chimp'), variable.name = 'Var1', value.name = 'Coef', presdf[,c(3,4,5)])
fdrmelt = melt(measure.vars = c('log10FDR_Human', 'log10FDR_Chimp'), variable.name = 'Var1', value.name = 'log10FDR', presdf[,c(6,7,5)])
toplot = cbind(Coefmelt, fdrmelt)

toplot$Species = gsub('.*_', '', toplot$Var1)
colnames(toplot) = make.unique(colnames(toplot))
toplot$condition = ifelse(toplot$Coef > 0, 'Enrichment', 'Depletion')

library(ggrepel)
toplot$CellType = gsub('Astro', 'Astrocyte', toplot$CellType)
toplot$CellType = gsub('Micro', 'Microglia', toplot$CellType)
pdf('Modern_Human_Enrichment_PEAKS_NegBinomREG_LengthLog2_LowCoverage.pdf', width = 10, height = 6)
ggscatter(toplot, x = 'log10FDR', y = 'CellType', color = 'condition', size = 'log10FDR', palette = c('blue', 'red')) +
geom_label_repel(data = toplot[toplot$log10FDR > 1.3,], aes(label = CellType), nudge_y = 0.02, nudge_x = -0.05) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
facet_wrap(~Species)
dev.off()


