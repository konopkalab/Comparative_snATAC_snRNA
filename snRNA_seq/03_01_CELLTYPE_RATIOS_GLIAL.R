rm(list = ls())
library(Seurat)
library(Matrix)
library(ggplot2)
library(plyr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(rio)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(harmony)
library(GeneOverlap)
source("utility_functions.R")
set.seed(1234)

####
## OLIGO RATIOS IN GLIA -- RNA-SEQ
####

allMeta = readRDS('RNA_ALL_META.RDS')


# HUMAN
allMetaH = allMeta[allMeta$Species == 'human',]
allMetaH = allMetaH[grepl('Oli|OPC|Mic|Ast|Glia', allMetaH$broadannot),] # Keep only glia
allMetaH$CellType = droplevels(allMetaH$broadannot)

allMetaH$CellType = gsub('Oli', 'MOL', allMetaH$CellType)

allDF = allMetaH %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaH %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaH$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfH = allDF2
allDfH$Species = 'Human'

# CHIMP
allMetaC = allMeta[allMeta$Species == 'chimp',]
allMetaC = allMetaC[grepl('Oli|OPC|Mic|Ast|Glia', allMetaC$broadannot),] # Keep only glia
allMetaC$CellType = droplevels(allMetaC$broadannot)

allMetaC$CellType = gsub('Oli', 'MOL', allMetaC$CellType)

allDF = allMetaC %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaC %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaC$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfC = allDF2
allDfC$Species = 'Chimp'

# MACAQUE
allMetaM = allMeta[allMeta$Species == 'macaque',]
allMetaM = allMetaM[grepl('Oli|OPC|Mic|Ast|Glia', allMetaM$broadannot),] # Keep only glia
allMetaM$CellType = droplevels(allMetaM$broadannot)

allMetaM$CellType = gsub('Oli', 'MOL', allMetaM$CellType)

allDF = allMetaM %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaM %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaM$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfM = allDF2
allDfM$Species = 'Macaque'

# Combine
combPlot = do.call(rbind, list(allDfH, allDfC, allDfM))
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'))

rnaOLRATIOS = combPlot
rnaOLRATIOS$type = 'RNA'

####
## OLIGO RATIOS IN GLIA -- ATAC-SEQ
####

allMeta = readRDS('ATAC_ALL_META.RDS')

# For glial proportions
allMeta = allMeta[grepl('MOL|OPC|Mic|Ast', allMeta$broadannot),] # Keep only glia
allMeta$CellType = allMeta$broadannot
allMeta$CellType = gsub('Astro', 'Astrocyte', allMeta$CellType)

# HUMAN
allMetaH = allMeta[allMeta$Species == 'Human',]

allDF = allMetaH %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaH %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaH$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfH = allDF2
allDfH$Species = 'Human'


# CHIMP
allMetaC = allMeta[allMeta$Species == 'Chimp',]
allDF = allMetaC %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaC %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaC$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfC = allDF2
allDfC$Species = 'Chimp'

# MACAQUE
allMetaM = allMeta[allMeta$Species == 'Macaque',]
allDF = allMetaM %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaM %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaM$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfM = allDF2
allDfM$Species = 'Macaque'

# Combine
combPlot = do.call(rbind, list(allDfH, allDfC, allDfM))
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'))
atacOLRATIOS = combPlot
atacOLRATIOS$type = 'ATAC'


####
## COMPARE GLIA
####

tmp = rbind(atacOLRATIOS, rnaOLRATIOS)
tmp[tmp$CellType == 'Ast', 'CellType'] = 'Astrocyte'
tmp[tmp$CellType == 'Mic', 'CellType'] = 'Microglia'

comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'))

# PLOT OPC AND MOL
tmp2 = tmp[tmp$CellType %in% c('OPC', 'MOL'),]
tmp2$CellType = factor(tmp2$CellType, levels = c('OPC', 'MOL'))

pdf('BA23_CB_OLIGO_PROPORTIONS.pdf')
ggboxplot(tmp2, x = 'Species', y = 'ratio', color = 'Species',
	 ylim = c(0,1.1)) +
	facet_wrap(~CellType) +

	geom_point(data = tmp2, aes(x = factor(Species), y = ratio, colour = factor(type)), position = position_jitter(width = 0.2)) +
	scale_colour_manual(values = c('red', 'orange', 'blue', 'darkgreen', 'pink')) +

	ylab('Proportion_To_Glia') + xlab('') + theme(text=element_text(size=20), legend.position = 'right') +
	rotate_x_text(90) + 
	stat_compare_means(comparisons = comps, method = 'wilcox.test')
dev.off()

# PLOT ALL
tmp$CellType = factor(tmp$CellType, levels = c('OPC', 'MOL', 'Astrocyte', 'Microglia'))

pdf('BA23_CB_GLIAL_PROPORTIONS.pdf', width = 10)
ggboxplot(tmp, x = 'Species', y = 'ratio', color = 'Species',
	 ylim = c(0,1.1)) +
	facet_wrap(~CellType, nrow = 1) +

	geom_point(data = tmp, aes(x = factor(Species), y = ratio, colour = factor(type)), position = position_jitter(width = 0.2)) +
	scale_colour_manual(values = c('red', 'orange', 'blue', 'darkgreen', 'pink')) +

	ylab('Proportion_To_Glia') + xlab('') + theme(text=element_text(size=20), legend.position = 'right') +
	rotate_x_text(90) + 
	stat_compare_means(comparisons = comps, method = 'wilcox.test')
dev.off()


# Test the ratios
library(lmtest)

# Human - Chimp
ctypes = unique(tmp$CellType)
pvalL = list()
for(i in 1:length(ctypes)){

	totest1 = tmp[tmp$CellType == ctypes[i] & tmp$Species %in% c('Human', 'Chimp'),]
	res0 = lm(ratio ~ type + Species, totest1)
	res1 = lm(ratio ~ type, totest1)
	pvalL[[i]] = lrtest(res0, res1)[2,5]
}

hcPval = unlist(pvalL)
hcComp = data.frame(CellType = ctypes, Comparison = 'Human-Chimp', Pvalue = hcPval)


# Human - Macaque
ctypes = unique(tmp$CellType)
pvalL = list()
for(i in 1:length(ctypes)){

	totest1 = tmp[tmp$CellType == ctypes[i] & tmp$Species %in% c('Human', 'Macaque'),]
	res0 = lm(ratio ~ type + Species, totest1)
	res1 = lm(ratio ~ type, totest1)
	pvalL[[i]] = lrtest(res0, res1)[2,5]
}

hmPval = unlist(pvalL)
hmComp = data.frame(CellType = ctypes, Comparison = 'Human-Macaque', Pvalue = hmPval)


# Chimp - Macaque
ctypes = unique(tmp$CellType)
pvalL = list()
for(i in 1:length(ctypes)){

	totest1 = tmp[tmp$CellType == ctypes[i] & tmp$Species %in% c('Chimp', 'Macaque'),]
	res0 = lm(ratio ~ type + Species, totest1)
	res1 = lm(ratio ~ type, totest1)
	pvalL[[i]] = lrtest(res0, res1)[2,5]
}

cmPval = unlist(pvalL)
cmComp = data.frame(CellType = ctypes, Comparison = 'Chimp-Macaque', Pvalue = cmPval)

# Combine and save the statistics
stats = do.call(rbind, list(hcComp, hmComp, cmComp))
stats$Group = 'Glia'

rio::export(stats, '~/workdir/pr3/05_CellType_Ratios/Glia_Subtype_Abundance_Stats.xlsx')
rio::export(tmp, '~/workdir/pr3/05_CellType_Ratios/Glia_Subtype_Fractions.xlsx')


