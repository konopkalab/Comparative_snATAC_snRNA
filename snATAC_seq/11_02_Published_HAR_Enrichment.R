rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggpubr)
library(reshape2)
library(data.table)
library(Signac)
library(Seurat)
library(dplyr)
library(tidyverse)
library(rio)
library(GenomicRanges)
library(liftOver)
library(ggrepel)
library('bedr')
source("utility_functions.R")

####
## MERGE HAR DATASETS
####

# Read HARs
doan = rio::import('Doan_2016_Compendium.xlsx')
gittel = rio::import('Gittelman_2015_hg19.xlsx')

# Prepare to merge
doan = doan[, c(1:3)]
gittel = gittel[, c(1:3)]
colnames(gittel) = colnames(doan)

doan = makeGRangesFromDataFrame(doan)
gittel = makeGRangesFromDataFrame(gittel)

# LiftOver HARs to hg38
hg19hg38 = import.chain('chain_files/forR/hg19ToHg38.over.chain')

doan_hg38 = liftOver(doan, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))
gittel_hg38 = liftOver(gittel, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))

# Merge all HARs
doan_hg38 = bedr.sort.region(doan_hg38, verbose = F)
gittel_hg38 = bedr.sort.region(gittel_hg38, verbose = F)

allHars = rbind(doan_hg38, gittel_hg38)
allHars = bedr.merge.region(allHars) # 3267 -> 3176
rio::export(allHars %>% dplyr::rename(width = 'names'), 'Doan_Combined_HARs_hg38.xlsx')

# Read HARs
allHars = rio::import('Doan_Combined_HARs_hg38.xlsx')
allHars = allHars[, c(1:3)]

####
## OVERLAP PEAKS AND HARs
####

# DARs (Differential CREs)

# Load background peaks
dars = readRDS('PSEUDOBULK_DARs_ALL.RDS')
dars$peak = dars$Gene

# All CREs
allpks = dars

# Evolved CREs
dars = dars[dars$Evolution != 'NS',]

# Overlap with all CREs
bcg_bed = allpks$peak %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))
bcg_bed = cbind(bcg_bed, peak = allpks$peak %>% gsub(':|-', '_', .))
bcg_bed = bcg_bed[!duplicated(bcg_bed$peak),]
bcg_bed = bedr.sort.region(bcg_bed, verbose = F)
bcg_bed = bedr(input = list(a = bcg_bed, b = allHars), method = "intersect", params = c("-loj"), verbose = F)

bcg_bed$start.b = as.numeric(bcg_bed$start.b)
bcg_bed$end.b = as.numeric(bcg_bed$end.b)
bcg_bed$ov = pmin(bcg_bed$V3, bcg_bed$end.b) - pmax(bcg_bed$V2, bcg_bed$start.b)
bcg_bed$ov = ifelse(bcg_bed$ov < 0, 0, bcg_bed$ov)

# Assign species specific CREs
hsdars = dars[dars$Evolution == 'Human_Specific', 'peak'] %>% gsub(':|-', '_', .) %>% unique

bcg_bed$is_hs = ifelse(bcg_bed$peak %in% hsdars, 'HS', 'NS')
bcg_bed$is_hs = factor(bcg_bed$is_hs, levels = c('NS', 'HS'))

# Length of CRE
bcg_bed$length = bcg_bed$V3 - bcg_bed$V2

####
## ENRICHMENT BY LOGISTIC REGRESSION
####

# Prepare for enrichment. Binarize the HAR overlap (i.e CREs either overlap with HAR (1) or don't (0).)
bcg_bed = bcg_bed[!duplicated(bcg_bed$peak),]
bcg_bed$ov = ifelse(bcg_bed$ov > 0, 1, 0)

# Test enrichment on all CREs. Fit length as covariate. Collect results in a data frame
hres = glm(ov ~ length + is_hs, data = bcg_bed, family = 'binomial')

human_pval = coef(summary(hres))[3,4]
OR_human = coef(summary(hres))[3,1]

####
## Run enrichment per major cell type
####

# Assign major cell type
allpks$Major = allpks$CellType
allpks[grepl('^L[0-9]', allpks$CellType), 'Major'] = 'Excitatory'
allpks[grepl('Upper|SST|PVALB|VIP|LAMP5', allpks$CellType), 'Major'] = 'Inhibitory'

dars$Major = dars$CellType
dars[grepl('^L[0-9]', dars$CellType), 'Major'] = 'Excitatory'
dars[grepl('Upper|SST|PVALB|VIP|LAMP5', dars$CellType), 'Major'] = 'Inhibitory'


# Loop over cell types and run enrichment per cell type
ctypes = names(table(allpks$Major))
hresCt = list()
cresCt = list()
subOvHARL = list()
for(i in 1:length(ctypes)){

	# Assign CREs as human-specific versus non-specific per cell type
	tested = allpks[allpks$Major == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)
	subOv = bcg_bed[bcg_bed$peak %in% tested,]

	hsubdars = dars[dars$Evolution == 'Human_Specific' & dars$Major == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)

	subOv$is_hs = ifelse(subOv$peak %in% hsubdars, 'HS', 'NS')
	subOv$is_hs = factor(subOv$is_hs, levels = c('NS', 'HS'))

	hresCt[[i]] = glm(ov ~ length + is_hs, data = subOv, family = 'binomial')

	# For table -- take only cre-har overlaps
	subOvHAR = subOv[subOv$start.b != '-1',]
	subOvHAR$CellType = ctypes[i]
	subOvHAR$ov = NULL
	subOvHAR$is_cs = NULL
	colnames(subOvHAR) = c('CRE_Chr', 'CRE_Start', 'CRE_End', 'CRE', 'HAR_Chr', 'HAR_Start', 'HAR_End', 'is_hs', 'Length', 'CellType')
	subOvHARL[[i]] = subOvHAR

	print(i)
}

# Save the final table
subOvHARDF = do.call(rbind, subOvHARL)
rio::export(subOvHARDF, 'HAR_CRE_FINAL_DF.xlsx')

# Collect results in a data frame
human_pval = sapply(hresCt, function(x){coef(summary(x))[3,4]})
OR_human = sapply(hresCt, function(x){coef(summary(x))[3,1]})
presdf = data.frame(human_pval = human_pval, Coef_Human = OR_human, Major = ctypes)

# Adjust for p-value using FDR
presdf$log10FDR_Human = -log10(p.adjust(presdf$human_pval, method = 'BH'))
presdf[is.na(presdf$log10FDR_Human) | is.nan(presdf$log10FDR_Human), 'log10FDR_Human'] = 0

# PLOT
pdf('PUBHAR_DAR_Enrichment_BINOMIAL.pdf')
ggscatter(presdf, x = 'log10FDR_Human', y = 'Coef_Human', color = 'blue', repel = T, font.label = c(20, 'bold', 'black'), label = 'Major', size = 5, palette = c('orange', 'blue')) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
theme(text = element_text(size=20, face = 'bold'), legend.pos = 'right')
dev.off()



