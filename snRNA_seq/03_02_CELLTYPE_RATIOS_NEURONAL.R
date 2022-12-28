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
## EXCITATORY
####

# RNA-SEQ #

allMeta = readRDS('RNA_ALL_META.RDS')
allMeta$Species = gsub('human', 'Human', allMeta$Species)
allMeta$Species = gsub('chimp', 'Chimp', allMeta$Species)
allMeta$Species = gsub('macaque', 'Macaque', allMeta$Species)


allMetaEXC = allMeta[allMeta$broadannot == 'Excitatory',]

# HUMAN
allMetaH = allMetaEXC[allMetaEXC$Species == 'Human',]
allMetaH$CellType = allMetaH$newannot

allDF = allMetaH %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaH %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaH$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfH = allDF2
allDfH$Species = 'Human'

# CHIMP
allMetaC = allMetaEXC[allMetaEXC$Species == 'Chimp',]
allMetaC$CellType = allMetaC$newannot

allDF = allMetaC %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaC %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaC$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfC = allDF2
allDfC$Species = 'Chimp'

# MACAQUE
allMetaM = allMetaEXC[allMetaEXC$Species == 'Macaque',]
allMetaM$CellType = allMetaM$newannot

allDF = allMetaM %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaM %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaM$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfM = allDF2
allDfM$Species = 'Macaque'

# Combine
combPlot = do.call(rbind, list(allDfH, allDfC, allDfM))
rnaOLRATIOS = combPlot
rnaOLRATIOS$type = 'RNA'


# ATAC-SEQ #

allMeta = readRDS('ATAC_ALL_META.RDS')
allMeta$broadannot = gsub('Exc', 'Excitatory', allMeta$broadannot)
allMeta$broadannot = gsub('Inh', 'Inhibitory', allMeta$broadannot)
allMetaEXC = allMeta[allMeta$broadannot == 'Excitatory',]
allMetaEXC$CellType = allMetaEXC$newannot

# HUMAN
allMetaH = allMetaEXC[allMetaEXC$Species == 'Human',]

allDF = allMetaH %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaH %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaH$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfH = allDF2
allDfH$Species = 'Human'


# CHIMP
allMetaC = allMetaEXC[allMetaEXC$Species == 'Chimp',]

allDF = allMetaC %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaC %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaC$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfC = allDF2
allDfC$Species = 'Chimp'

# MACAQUE
allMetaM = allMetaEXC[allMetaEXC$Species == 'Macaque',]
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


# COMPARE SUBTYPES
tmp = rbind(atacOLRATIOS, rnaOLRATIOS)
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'))

pdf('BA23_CB_EXC_SUBTYPE_PROPORTIONS.pdf', width = 10)
ggboxplot(tmp, x = 'Species', y = 'ratio', color = 'Species',
	 ylim = c(0,0.8)) +
	facet_wrap(~CellType, nrow = 2) +

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
stats$Group = 'Excitatory'

rio::export(stats, 'Excitatory_Subtype_Abundance_Stats.xlsx')
rio::export(tmp, 'Excitatory_Subtype_Fractions.xlsx')

####
## INHIBITORY
####

# RNA-SEQ #

allMeta = readRDS('RNA_ALL_META.RDS')
allMeta$Species = gsub('human', 'Human', allMeta$Species)
allMeta$Species = gsub('chimp', 'Chimp', allMeta$Species)
allMeta$Species = gsub('macaque', 'Macaque', allMeta$Species)

allMetaEXC = allMeta[allMeta$broadannot == 'Inhibitory',]

# HUMAN
allMetaH = allMetaEXC[allMetaEXC$Species == 'Human',]
allMetaH$CellType = allMetaH$newannot

allDF = allMetaH %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaH %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaH$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfH = allDF2
allDfH$Species = 'Human'

# CHIMP
allMetaC = allMetaEXC[allMetaEXC$Species == 'Chimp',]
allMetaC$CellType = allMetaC$newannot

allDF = allMetaC %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaC %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaC$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfC = allDF2
allDfC$Species = 'Chimp'

# MACAQUE
allMetaM = allMetaEXC[allMetaEXC$Species == 'Macaque',]
allMetaM$CellType = allMetaM$newannot

allDF = allMetaM %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaM %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaM$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfM = allDF2
allDfM$Species = 'Macaque'

# Plot Transcriptome
combPlot = do.call(rbind, list(allDfH, allDfC, allDfM))
rnaOLRATIOS = combPlot
rnaOLRATIOS$type = 'RNA'

# ATAC-SEQ #

allMeta = readRDS('ATAC_ALL_META.RDS')
allMeta$broadannot = gsub('Exc', 'Excitatory', allMeta$broadannot)
allMeta$broadannot = gsub('Inh', 'Inhibitory', allMeta$broadannot)
allMetaEXC = allMeta[allMeta$broadannot == 'Inhibitory',]
allMetaEXC$CellType = allMetaEXC$newannot

# HUMAN
allMetaH = allMetaEXC[allMetaEXC$Species == 'Human',]

allDF = allMetaH %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaH %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaH$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfH = allDF2
allDfH$Species = 'Human'


# CHIMP
allMetaC = allMetaEXC[allMetaEXC$Species == 'Chimp',]

allDF = allMetaC %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaC %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaC$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfC = allDF2
allDfC$Species = 'Chimp'

# MACAQUE
allMetaM = allMetaEXC[allMetaEXC$Species == 'Macaque',]
allDF = allMetaM %>% group_by(Sample, CellType) %>% group_keys %>% as.data.frame
allDF$size = allMetaM %>% group_by(Sample, CellType) %>% group_size
allDFTotSize = table(allMetaM$Sample) %>% as.data.frame %>% dplyr::rename(Sample='Var1',totsize='Freq')
allDF2 = merge(allDF, allDFTotSize)
allDF2$ratio = allDF2$size / allDF2$totsize

allDfM = allDF2
allDfM$Species = 'Macaque'

# Plot
combPlot = do.call(rbind, list(allDfH, allDfC, allDfM))
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'))
atacOLRATIOS = combPlot
atacOLRATIOS$type = 'ATAC'



# COMPARE SUBTYPES
tmp = rbind(atacOLRATIOS, rnaOLRATIOS)
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'))

# PLOT
pdf('BA23_CB_INH_SUBTYPE_PROPORTIONS.pdf', width = 10)
ggboxplot(tmp, x = 'Species', y = 'ratio', color = 'Species',
	 ylim = c(0,0.8)) +
	facet_wrap(~CellType, nrow = 2) +

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
stats$Group = 'Inhibitory'

rio::export(stats, 'Inhibitory_Subtype_Abundance_Stats.xlsx')
rio::export(tmp, 'Inhibitory_Subtype_Fractions.xlsx')



