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
require(pscl)
require(MASS)
require(boot)
library('bedr')
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

####
## MODERN VARIANT DATASET FROM 12_01
####

hvarsCons = readRDS('CommonNeanderthal_N3_D1_Prufer2014_Table.RDS')

# Keep only the ones >90% frequency
hvarsCons = hvars[hvars$Human_Maj_1000G_PHASE3 > 0.9, ]

# Prepare for overlap
hvarsbed = hvarsCons[, c(1:2)]
hvarsbed$end = hvarsbed[,2] + 1
colnames(hvarsbed) = c('chr', 'start', 'end')
hvarsbed$chr = gsub('^', 'chr', hvarsbed$chr)

hvarsbed = makeGRangesFromDataFrame(hvarsbed)

# LiftOver to hg38
hg19hg38 = import.chain('~/workdir/reference_genomes/chain_files/forR/hg19ToHg38.over.chain')
hvars_hg38 = liftOver(hvarsbed, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))

# Sort
hvars_hg38 = bedr.sort.region(hvars_hg38, verbose = F)
saveRDS(hvars_hg38, 'CommonNeanderthal_N3_Prufer2014_Table_HG38_MIN09_N3_D1.RDS')

####
## ENRICHMENT BY REGRESSION
####

# Load background peaks
dars = readRDS('~/workdir/pr3/atacseq/human_chimp_macaque/12_DARs/PSEUDOBULK_DARs_ALL.RDS')
dars$peak = dars$Gene
dars$Major = dars$CellType
dars[grepl('^L[0-9]', dars$CellType), 'Major'] = 'Excitatory'
dars[grepl('Upper|SST|PVALB|VIP|LAMP5', dars$CellType), 'Major'] = 'Inhibitory'

# All CREs in BED format
bcg = dars

bcg_bed = bcg$peak %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))
bcg = cbind(bcg, bcg_bed)
hsAll = bcg[bcg$Evolution == 'Human_Specific', 'peak']
table(bcg$CellType)

# Get CREs in BED format (remove duplicate entries)
bcgAll = bcg[,c('V1','V2','V3')] %>% add_column(peak = paste0(.[,1], ':', .[,2], '-', .[,3]))
bcgAll = bcgAll[!duplicated(bcgAll$peak),]
bcgAll = bedr.sort.region(bcgAll, verbose = F)

# Overlap CREs with modern variants
bcg_bedr = bedr(input = list(a = bcgAll, b = hvars_hg38), method = "intersect", params = "-loj", verbose = F)
bcg_bedr$modern_ov = ifelse(bcg_bedr$start.b > -1, 1, 0)
bcg_bedr$length = bcg_bedr[,3] - bcg_bedr[,2]

# Sum number of modern variant overlap per CRE
bcgAgg = aggregate(modern_ov ~ peak, bcg_bedr, FUN = sum)
bcgAgg$length = bcg_bedr[match(bcgAgg$peak, bcg_bedr$peak), 'length']
bcgAgg$is_hs = ifelse(bcgAgg$peak %in% hsAll, 'HS', 'NS')
bcgAgg$is_hs = factor(bcgAgg$is_hs, levels = c('NS', 'HS'))

# Overall enrichment
humanAll = glm.nb(modern_ov ~ length + is_hs, data = bcgAgg)

####
## ENRICH ALL CELL TYPES
####

darsEXC = dars

ctypes = names(table(dars$CellType))
human_h1 = list()
bcgSaves = list()
subOvMODERNL = list()
for(i in 1:length(ctypes)){

	# Get CREs
	bcg_ct = bcg[bcg$CellType == ctypes[i], c('V1','V2','V3')] %>% add_column(peak = paste0(.[,1], ':', .[,2], '-', .[,3]))
	hs_ct = bcg[bcg$CellType == ctypes[i] & bcg$Evolution == 'Human_Specific', 'peak']

	bcgAgg_ct = bcgAgg[bcgAgg$peak %in% bcg_ct$peak,]
	bcgAgg_ct$is_hs = ifelse(bcgAgg_ct$peak %in% hs_ct, 'HS', 'NS')
	bcgAgg_ct$is_hs = factor(bcgAgg_ct$is_hs, levels = c('NS', 'HS'))

	human_h1[[i]] <- glm.nb(modern_ov ~ length + is_hs, data = bcgAgg_ct)

	# For table -- take only cre-modern variant overlaps
	subOvMODERN = bcg_bedr[bcg_bedr$peak %in% bcg_ct$peak,]
	subOvMODERN = subOvMODERN[subOvMODERN$start.b != '-1',]
	subOvMODERN$CellType = ctypes[i]
	subOvMODERN$strand = NULL
	subOvMODERN$modern_ov = NULL
	subOvMODERN$end.b = NULL
	subOvMODERN$width = NULL
	subOvMODERN$length = subOvMODERN[,3] - subOvMODERN[,2]
	subOvMODERN$is_hs = ifelse(subOvMODERN$peak %in% hs_ct, 'HS', 'NS')
	colnames(subOvMODERN) = c('CRE_Chr', 'CRE_Start', 'CRE_End', 'CRE', 'Modern_Chr', 'Modern_Pos', 'Length', 'CellType', 'is_hs')
	subOvMODERNL[[i]] = subOvMODERN

	print(i)
}


# Save overlap as database
subOvMODERN_DF = do.call(rbind, subOvMODERNL)
rio::export(subOvMODERN_DF, 'Peak_ModernVar_Overlap_ALL.xlsx')


# Save and plot enrichments
human_pval = sapply(human_h1, function(x){coef(summary(x))[3,4]})
OR_human = sapply(human_h1, function(x){coef(summary(x))[3,1]})

presdf = data.frame(human_pval = human_pval, Coef_Human = OR_human, CellType = c(ctypes))
presdf$log10FDR_Human = -log10(p.adjust(presdf$human_pval, method = 'BH'))
presdf[is.na(presdf$log10FDR_Human) | is.nan(presdf$log10FDR_Human), 'log10FDR_Human'] = 0

presdf$Major = 'Inhibitory'
presdf[grepl('^L[0-9]', presdf$CellType), 'Major'] = 'Excitatory'
presdf[grepl('Oli|Mic|OPC|Ast', presdf$CellType), 'Major'] = 'Glia'
presdf$Major = factor(presdf$Major, levels = c('Glia', 'Inhibitory', 'Excitatory'))

pdf('Modern_Human_Enrichment_ALL_N3_D1.pdf', width = 8, height = 6)
ggscatter(presdf, x = 'log10FDR_Human', y = 'Coef_Human', color = 'Major', palette = c('Blues'), size = 5) +
geom_label_repel(data = presdf[presdf$log10FDR_Human > 1.3,],
		aes(x=log10FDR_Human, y = Coef_Human, label = CellType), size = 7, label.size = NA,
		fontface = 'bold', fill = alpha(c("white"),0)) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
theme(text = element_text(size=20, face = 'bold'), legend.pos = 'right')
dev.off()





