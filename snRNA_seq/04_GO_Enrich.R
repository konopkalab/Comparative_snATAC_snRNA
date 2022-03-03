rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(rio)
library(GeneOverlap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(grid)
library(gridExtra)
library(showtext)
source("utility_functions.R")

####
## ENRICH EXCITATORY
####

uni = import("Excitatory_AllGenes_DEGs.xlsx")
degs = import("Excitatory_SignGenes_DEGs.xlsx")
humanUP = degs[degs$Regulation == 'Human_UP',]
humanDOWN = degs[degs$Regulation == 'Human_DOWN',]

# Run Enrichment HUMAN UP
ctypes = names(table(degs$cluster))
enriched = list()
for(i in 1:length(ctypes)){

	if(length(humanUP[humanUP$cluster == ctypes[i], 'Gene']) <= 1){next}

	deg1 = humanUP[humanUP$cluster == ctypes[i], 'Gene']
	uni1 = uni[uni$cluster == ctypes[i], 'Gene']

	print(length(deg1))
	print(length(uni1))

	res = GOenrich(deg1, uni1)

	if(is.null(res) == T){next}
	if(nrow(res) == 0){next}

	res$CellType = ctypes[i]
	enriched[[i]] = res
}

rnadf = do.call(rbind, enriched)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_Human_UP_DEGs.xlsx")

# Run Enrichment HUMAN DOWN
ctypes = names(table(degs$cluster))
enriched = list()
for(i in 1:length(ctypes)){

	if(length(humanDOWN[humanDOWN$cluster == ctypes[i], 'Gene']) <= 1){next}

	deg1 = humanDOWN[humanDOWN$cluster == ctypes[i], 'Gene']
	uni1 = uni[uni$cluster == ctypes[i], 'Gene']

	print(length(deg1))
	print(length(uni1))

	res = GOenrich(deg1, uni1)

	if(is.null(res) == T){next}
	if(nrow(res) == 0){next}

	res$CellType = ctypes[i]
	enriched[[i]] = res
}

rnadf = do.call(rbind, enriched)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_Human_DOWN_EXC.xlsx")

####
## ENRICH INHIBITORY
####

uni = import("Inhibitory_AllGenes_DEGs.xlsx")
degs = import("Inhibitory_SignGenes_DEGs.xlsx")
humanUP = degs[degs$Regulation == 'Human_UP',]
humanDOWN = degs[degs$Regulation == 'Human_DOWN',]

# Run Enrichment HUMAN UP
ctypes = names(table(degs$cluster))
enriched = list()
for(i in 1:length(ctypes)){

	if(length(humanUP[humanUP$cluster == ctypes[i], 'Gene']) <= 1){next}

	deg1 = humanUP[humanUP$cluster == ctypes[i], 'Gene']
	uni1 = uni[uni$cluster == ctypes[i], 'Gene']

	print(length(deg1))
	print(length(uni1))

	res = GOenrich(deg1, uni1)

	if(is.null(res) == T){next}
	if(nrow(res) == 0){next}

	res$CellType = ctypes[i]
	enriched[[i]] = res
}

rnadf = do.call(rbind, enriched)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_Human_UP_INH.xlsx")

# Run Enrichment HUMAN DOWN
ctypes = names(table(degs$cluster))
enriched = list()
for(i in 1:length(ctypes)){

	if(length(humanDOWN[humanDOWN$cluster == ctypes[i], 'Gene']) <= 1){next}

	deg1 = humanDOWN[humanDOWN$cluster == ctypes[i], 'Gene']
	uni1 = uni[uni$cluster == ctypes[i], 'Gene']

	print(length(deg1))
	print(length(uni1))

	res = GOenrich(deg1, uni1)

	if(is.null(res) == T){next}
	if(nrow(res) == 0){next}

	res$CellType = ctypes[i]
	enriched[[i]] = res
}

rnadf = do.call(rbind, enriched)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_Human_DOWN_INH.xlsx")

####
## ENRICH GLIA
####

oli = import("Oligodendrocyte_AllGenes_DEGs.xlsx")
ast = import("Astrocyte_AllGenes_DEGs.xlsx")
mic = import("Microglia_AllGenes_DEGs.xlsx")

uni = rbind(oli, ast, mic)

olidegs = import("Oligodendrocyte_SignGenes_DEGs.xlsx")
astdegs = import("Astrocyte_SignGenes_DEGs.xlsx")
micdegs = import("Microglia_SignGenes_DEGs.xlsx")

gliDegs = rbind(olidegs, astdegs, micdegs)

humanUP = gliDegs[gliDegs$Regulation == 'Human_UP',]
humanDOWN = gliDegs[gliDegs$Regulation == 'Human_DOWN',]

# Run Enrichment HUMAN UP
ctypes = names(table(gliDegs$cluster))
enriched = list()
for(i in 1:length(ctypes)){

	if(length(humanUP[humanUP$cluster == ctypes[i], 'Gene']) <= 1){next}

	deg1 = humanUP[humanUP$cluster == ctypes[i], 'Gene']
	uni1 = uni[uni$cluster == ctypes[i], 'Gene']

	print(length(deg1))
	print(length(uni1))

	res = GOenrich(deg1, uni1)

	if(is.null(res) == T){next}
	if(nrow(res) == 0){next}

	res$CellType = ctypes[i]
	enriched[[i]] = res
}

rnadf = do.call(rbind, enriched)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_Human_UP_GLIA.xlsx")



# Run Enrichment HUMAN DOWN
ctypes = names(table(gliDegs$cluster))
enriched = list()
for(i in 1:length(ctypes)){

	if(length(humanDOWN[humanDOWN$cluster == ctypes[i], 'Gene']) <= 1){next}

	deg1 = humanDOWN[humanDOWN$cluster == ctypes[i], 'Gene']
	uni1 = uni[uni$cluster == ctypes[i], 'Gene']

	print(length(deg1))
	print(length(uni1))

	res = GOenrich(deg1, uni1)

	if(is.null(res) == T){next}
	if(nrow(res) == 0){next}

	res$CellType = ctypes[i]
	enriched[[i]] = res
}

rnadf = do.call(rbind, enriched)
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_Human_DOWN_GLIA.xlsx")


####
## PLOT
####

# PLOT LAMP5
library(showtext)

inhUP = import("GO_Enrich_Human_UP_INH.xlsx")
inhDOWN = import("GO_Enrich_Human_DOWN_INH.xlsx")

lamp5_UP = inhUP[grepl('LAMP5', inhUP$CellType),] %>% add_column(type = 'UP')
lamp5_DOWN = inhDOWN[grepl('LAMP5', inhDOWN$CellType),] %>% add_column(type = 'DOWN')

pltgns = rbind(lamp5_UP, lamp5_DOWN)

# Plot only axon / dendrite
pltgns = pltgns[grepl('axon|dendrite|pre|post', pltgns$Description),]
pltgns$log10FDR = -log10(pltgns$p.adjust)

pdf('GO_ENRICH_LAMP5.pdf', width = 15)
print( ggscatter(pltgns, x = 'log10FDR', y = 'Description', size = 6,
	color = 'CellType', xlim = c(1,max(pltgns$log10FDR)), palette = c('blue', 'orange', 'purple')) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
theme_bw() +
facet_wrap(~type, drop = T, scales = 'free', nrow = 2) +
theme(text=element_text(size=30, face = 'bold')) + ylab('') +
xlab(expression('-log'[10]*'(FDR)')))
dev.off()



# PLOT OPC2
library(showtext)

olUP = import("GO_Enrich_Human_UP_GLIA.xlsx")
olDOWN = import("GO_Enrich_Human_DOWN_GLIA.xlsx")

olUP = olUP[grepl('OPC2', olUP$CellType),] %>% add_column(type = 'UP')
olDOWN = olDOWN[grepl('OPC2', olDOWN$CellType),] %>% add_column(type = 'DOWN')

pltgns = rbind(olUP, olDOWN)

# Plot only axon / dendrite
pltgns = pltgns[grepl('channel|synapse', pltgns$Description),]
pltgns$log10FDR = -log10(pltgns$p.adjust)

pdf('GO_ENRICH_OPC2.pdf', width = 15)
print( ggscatter(pltgns, x = 'log10FDR', y = 'Description', size = 6,
	color = 'CellType', xlim = c(1,max(pltgns$log10FDR)), palette = c('red')) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
theme_bw() +
facet_wrap(~type, drop = T, scales = 'free', nrow = 2) +
theme(text=element_text(size=30, face = 'bold')) + ylab('') +
xlab(expression('-log'[10]*'(FDR)')))
dev.off()


