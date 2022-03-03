rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(rio)
library(GeneOverlap)
library(Seurat)
library(ggrepel)
source("utility_functions.R")

####
## PREPARE DATASETS
####

# DEGs #
alldegsBCG = import('COMBINED_AllGenes_DEGs.xlsx')
alldegs = import('COMBINED_SignGenes_DEGs.xlsx')

# Remove SST-NPY
humandeg = alldegs[alldegs$Evolution == 'Human_Specific',] %>% select(Gene,cluster) %>% rename(Class = 'cluster') %>% as.data.frame
chimpdeg = alldegs[alldegs$Evolution == 'Chimp_Specific',] %>% select(Gene,cluster) %>% rename(Class = 'cluster') %>% as.data.frame
macaquedeg = alldegs[alldegs$MvsLCA != 'NS',] %>% select(Gene,cluster) %>% rename(Class = 'cluster') %>% as.data.frame

humandeg$Class = factor(humandeg$Class)
chimpdeg$Class = factor(chimpdeg$Class)
macaquedeg$Class = factor(macaquedeg$Class)

human_list = split(humandeg, humandeg$Class)
chimp_list = split(chimpdeg, chimpdeg$Class)
macaque_list = split(macaquedeg, macaquedeg$Class)

human_list = lapply(human_list, function(x){x$Gene})
chimp_list = lapply(chimp_list, function(x){x$Gene})
macaque_list = lapply(macaque_list, function(x){x$Gene})

names(human_list) = paste('Human', names(human_list), sep = '_')
names(chimp_list) = paste('Chimp', names(chimp_list), sep = '_')
names(macaque_list) = paste('Macaque', names(macaque_list), sep = '_')

# Velmeshev et al., 2019 #
velm_deg = rio::import_list("Velmeshev_STable4.xlsx")
velm_deg = velm_deg[[1]]

# Reshape for enrichment
velm_asd = velm_deg[,c(1,3)] %>% as.data.frame %>% rename(DEFINITION = 'Cell type') %>% rename(Gene = 'Gene name')

velm_list = split(velm_asd, velm_asd$DEFINITION)
velm_list = lapply(velm_list, function(x){x$Gene})

# Nagy et al., 2020 #
nagy_deg = rio::import_list("Nagy_2020_STable.xlsx")
nagy_deg = nagy_deg[5:28]
nagyArranged = list()
for(i in 1:length(nagy_deg)){

	q = nagy_deg[[i]]
	ctype = gsub('.*results for ', '', colnames(q)[1]) %>% gsub(' \\(con.*', '', .)
	colnames(q) = q[2,]
	q = q[-c(1:2),]
	q$CellType = ctype
	q$p.adjust = as.numeric(q$p.adjust)
	q = q[q$p.adjust < 0.1,]
	nagyArranged[[i]] = q
}

nagyArranged = do.call(rbind, nagyArranged)
nagyArranged = nagyArranged[nagyArranged$CellType != 'Endo',]
nagyArranged$Major = 'Unk'
nagyArranged[grepl('Ex', nagyArranged$CellType), 'Major'] = 'Excitatory'
nagyArranged[grepl('Inh', nagyArranged$CellType), 'Major'] = 'Inhibitory'
nagyArranged[grepl('OPC', nagyArranged$CellType), 'Major'] = 'OPC'
nagyArranged[grepl('Ast', nagyArranged$CellType), 'Major'] = 'Astrocyte'

# Reshape for enrichment
nagyArranged2 = nagyArranged[,c('Gene', 'Major')] %>% as.data.frame %>% rename(DEFINITION = 'Major')

nagy_list = split(nagyArranged2, nagyArranged2$DEFINITION)
nagy_list = lapply(nagy_list, function(x){x$Gene})

# Mathys et al., 2019 #
mathys_deg = rio::import_list("Mathys_STable2.xlsx")
mathys_deg = mathys_deg[2:7]
mathysArranged = list()
for(i in 1:length(mathys_deg)){

	q = mathys_deg[[i]]
	colnames(q) = q[1,]
	q = q[-1,]
	q = q[, 1:9] # Keep only pathology vs ctl comparison
	q$CellType = names(mathys_deg)[i]
	colnames(q)[1] = 'Gene'
	q$DEGs.Ind.Mix.models = as.logical(q$DEGs.Ind.Mix.models)
	q = q[q$DEGs.Ind.Mix.models == 1, c('Gene', 'CellType')]
	mathysArranged[[i]] = q
}

mathysArranged = do.call(rbind, mathysArranged)
mathysArranged = mathysArranged[mathysArranged$CellType != 'Endo',]
mathysArranged$Major = 'Unk'
mathysArranged[grepl('Ex', mathysArranged$CellType), 'Major'] = 'Excitatory'
mathysArranged[grepl('In', mathysArranged$CellType), 'Major'] = 'Inhibitory'
mathysArranged[grepl('Opc', mathysArranged$CellType), 'Major'] = 'OPC'
mathysArranged[grepl('Oli', mathysArranged$CellType), 'Major'] = 'Oligodendrocyte'
mathysArranged[grepl('Mic', mathysArranged$CellType), 'Major'] = 'Microglia'
mathysArranged[grepl('Ast', mathysArranged$CellType), 'Major'] = 'Astrocyte'

# Reshape for enrichment
mathysArranged2 = mathysArranged[,c('Gene', 'Major')] %>% as.data.frame %>% rename(DEFINITION = 'Major')

export(mathysArranged2, 'mathysArranged.xlsx')

mathys_list = split(mathysArranged2, mathysArranged2$DEFINITION)
mathys_list = lapply(mathys_list, function(x){x$Gene})


####
## VELMESHEV ENRICHMENTS (ASD)
####

library(GeneOverlap)

# Assign background
bcg = unique(c(velm_asd$gene_name, unique(alldegsBCG$Gene))) %>% length

# Create list of cell type matches
velmCT = list(c('AST-FB', 'AST-PP'), c('OPC'), c('Oligodendrocytes'), c('Microglia'), c('IN-SST'), c('IN-PV'), c('IN-VIP'), c('IN-SV2C'), c('L2/3'), c('L4'), c('L5/6', 'L5/6-CC'))

ba23CT = list(grep('Ast', names(human_list)),
		grep('Human_OPC|Human_COP', names(human_list)),
		grep('Human_MOL|Human_NFOL', names(human_list)),
		grep('Human_Microglia', names(human_list)),
		grep('Human_SST', names(human_list)),
		grep('Human_PVALB', names(human_list)),
		grep('Human_VIP', names(human_list)),
		grep('Human_LAMP5', names(human_list)),
		grep('^Human_L2-3', names(human_list)),
		grep('^Human_L3-5|^Human_L4', names(human_list)),
		grep('^Human_L5-6', names(human_list)))
nms = c('AST', 'OPC', 'OL', 'MIC', 'SST', 'PV', 'VIP', 'LAMP5', 'L2-3', 'L4', 'L5-6')

resL = list()
for(i in 1:length(ba23CT)){


	velm_listSUB = velm_list[names(velm_list) %in% velmCT[[i]]]
	human_listSUB = human_list[ba23CT[[i]]]
	chimp_listSUB = chimp_list[ba23CT[[i]]]
	macaque_listSUB = macaque_list[ba23CT[[i]]]

	# Run Enrichment
	degs_list = c(human_listSUB, chimp_listSUB, macaque_listSUB)

	resgom = newGOM(degs_list, velm_listSUB, genome.size=bcg)
	pvalmat = getMatrix(resgom, name="pval")
	pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
	pvalmelt = melt(pvalmat, value.name = 'pval')
	pvalmelt$Species = gsub('_.*', '', pvalmelt$Var1)
	pvalmelt$CellType = gsub('Human_|Chimp_|Macaque_', '', pvalmelt$Var1)
	pvalmelt$log10_FDR = -log10(pvalmelt$pval)

	oddsmat = getMatrix(resgom, name="odds.ratio")
	oddsmelt = melt(oddsmat, value.name = 'odds.ratio')

	resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio)
	resdf$log10_round = round(resdf$log10_FDR, digits = 2)
	resdf$OR_round = round(resdf$OR, digits = 2)
	resdf$FDR = ifelse(resdf$log10_FDR > 2, '<0.01',
				ifelse(resdf$log10_FDR > 1.3, '<0.05', '>0.05'))
	resdf$FDR = factor(resdf$FDR, levels = c('<0.01', '<0.05', '>0.05'))

	resdf$Species = gsub('Macaque', 'MvsHC', resdf$Species)
	resdf$Species = factor(resdf$Species, levels = c('Human', 'Chimp', 'MvsHC'))

	resL[[i]] = resdf

}

resdfs = do.call(rbind, resL)
export(resdfs, 'ASD_VELMESHEV_Stats.xlsx')

pdf('Species_VELMESHEV_GeneRegulation_RNA_PUB.pdf', width = 20, height = 3)
ggscatter(resdfs, x = 'log10_round', y = 'Species', color = 'Species', size = 3, palette = c('blue', 'orange', 'darkgreen'), alpha = 0.5) +
rotate_x_text(90) + theme(text=element_text(size=20)) + ylab(expression('-log'[10]*'(FDR)')) + xlab('') + NoLegend() +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
facet_wrap(~Var2, nrow = 1)
dev.off()



## VENN DIAGRAM ##
# Load library
library(VennDiagram)
 
resdfs2 = resdfs[resdfs$log10_round > -log10(0.05),]


# VENN
set1 <- resdfs2[resdfs2$Species == 'Human', 'CellType']
set2 <- resdfs2[resdfs2$Species == 'Chimp', 'CellType']
set3 <- resdfs2[resdfs2$Species == 'MvsHC', 'CellType']

setlist = list(set1, set2, set3)
cols = c('blue', 'orange', 'darkgreen')
nms = c('Human', 'Chimpanzee', 'MvsHC')
names(setlist) = nms

library(ggvenn)
pdf('ASD_SINGLECELL_VENN.pdf')
ggvenn(setlist, fill_color = cols, stroke_size = 0, set_name_size = 5,text_size = 15, digits = 0, stroke_alpha = 0.5, show_percentage = F)
dev.off()


####
## NAGY ENRICHMENTS (MDD)
####

library(GeneOverlap)

# Assign background
bcg = unique(c(nagyArranged2$Gene, unique(alldegsBCG$Gene))) %>% length

# Create list of cell type matches
nagyCT = list('Astrocyte', 'Excitatory', 'Inhibitory', 'OPC')

ba23CT = list(grep('Ast', names(human_list)),
		grep('^Human_L2-3|^Human_L3-5|^Human_L4|^Human_L5-6', names(human_list)),
		grep('Human_SST|Human_PVALB|Human_VIP|Human_LAMP5', names(human_list)),
		grep('Human_OPC|Human_COP', names(human_list)))

nms = names(nagy_list)
resL = list()
for(i in 1:length(ba23CT)){


	nagy_listSUB = nagy_list[names(nagy_list) %in% nagyCT[[i]]]
	human_listSUB = human_list[ba23CT[[i]]]
	chimp_listSUB = chimp_list[ba23CT[[i]]]
	macaque_listSUB = macaque_list[ba23CT[[i]]]

	# Run Enrichment
	degs_list = c(human_listSUB, chimp_listSUB, macaque_listSUB)

	resgom = newGOM(degs_list, nagy_listSUB, genome.size=bcg)
	pvalmat = getMatrix(resgom, name="pval")
	pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
	pvalmelt = melt(pvalmat, value.name = 'pval')
	pvalmelt$Species = gsub('_.*', '', pvalmelt$Var1)
	pvalmelt$CellType = gsub('Human_|Chimp_|Macaque_', '', pvalmelt$Var1)
	pvalmelt$log10_FDR = -log10(pvalmelt$pval)

	oddsmat = getMatrix(resgom, name="odds.ratio")
	oddsmelt = melt(oddsmat, value.name = 'odds.ratio')

	resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio)
	resdf$log10_round = round(resdf$log10_FDR, digits = 2)
	resdf$OR_round = round(resdf$OR, digits = 2)
	resdf$FDR = ifelse(resdf$log10_FDR > 2, '<0.01',
				ifelse(resdf$log10_FDR > 1.3, '<0.05', '>0.05'))
	resdf$FDR = factor(resdf$FDR, levels = c('<0.01', '<0.05', '>0.05'))

	resdf$Species = gsub('Macaque', 'MvsHC', resdf$Species)
	resdf$Species = factor(resdf$Species, levels = c('Human', 'Chimp', 'MvsHC'))

	resL[[i]] = resdf
}

resdfs = do.call(rbind, resL)
export(resdfs, 'NAGY_Stats.xlsx')

####
## MATHYS ENRICHMENTS (ALZ)
####

library(GeneOverlap)

# Assign background
bcg = unique(c(mathysArranged2$Gene, unique(alldegsBCG$Gene))) %>% length

# Create list of cell type matches
mathysCT = list('Astrocyte', 'Excitatory', 'Inhibitory', 'Microglia', 'Oligodendrocyte', 'OPC')

ba23CT = list(grep('Ast', names(human_list)),
		grep('^Human_L2-3|^Human_L3-5|^Human_L4|^Human_L5-6', names(human_list)),
		grep('Human_SST|Human_PVALB|Human_VIP|Human_LAMP5|Human_PAX6|Human_LHX6|Human_ADARB2', names(human_list)),
		grep('Human_Microglia', names(human_list)),
		grep('Human_MOL|Human_NFOL', names(human_list)),
		grep('Human_OPC|Human_COP', names(human_list)))

nms = names(mathys_list)

resL = list()
for(i in 1:length(ba23CT)){

	mathys_listSUB = mathys_list[names(mathys_list) %in% mathysCT[[i]]]
	human_listSUB = human_list[ba23CT[[i]]]
	chimp_listSUB = chimp_list[ba23CT[[i]]]
	macaque_listSUB = macaque_list[ba23CT[[i]]]

	# Run Enrichment
	
	degs_list = c(human_listSUB, chimp_listSUB, macaque_listSUB)

	resgom = newGOM(degs_list, mathys_listSUB, genome.size=bcg)
	pvalmat = getMatrix(resgom, name="pval")
	pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
	pvalmelt = melt(pvalmat, value.name = 'pval')
	pvalmelt$Species = gsub('_.*', '', pvalmelt$Var1)
	pvalmelt$CellType = gsub('Human_|Chimp_|Macaque_', '', pvalmelt$Var1)
	pvalmelt$log10_FDR = -log10(pvalmelt$pval)

	oddsmat = getMatrix(resgom, name="odds.ratio")
	oddsmelt = melt(oddsmat, value.name = 'odds.ratio')

	resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio)
	resdf$log10_round = round(resdf$log10_FDR, digits = 2)
	resdf$OR_round = round(resdf$OR, digits = 2)
	resdf$FDR = ifelse(resdf$log10_FDR > 2, '<0.01',
				ifelse(resdf$log10_FDR > 1.3, '<0.05', '>0.05'))
	resdf$FDR = factor(resdf$FDR, levels = c('<0.01', '<0.05', '>0.05'))

	resdf$Species = gsub('Macaque', 'MvsHC', resdf$Species)
	resdf$Species = factor(resdf$Species, levels = c('Human', 'Chimp', 'MvsHC'))


	pdf(paste0(nms[i], '_Mathys_Enrich_Heat.pdf'))
	print( ggplot(resdf, aes(Species, CellType, fill = FDR))+
	 xlab('') + ylab('') +
	 geom_tile(color = "white") +
	 scale_fill_manual(values = c('<0.01'='red', '<0.05'='pink', '>0.05'='grey')) +
	 facet_wrap(~Var2) +
	 rotate_x_text(90) +
	 theme_minimal() + rotate_x_text(90) +
	 theme(text = element_text(size=20)) +
	 geom_text(aes(label = log10_round)) )
	dev.off()

	resL[[i]] = resdf
}


resdfs = do.call(rbind, resL)
export(resdfs, 'MATHYS_Stats.xlsx')


resdfs$Var2 = factor(resdfs$Var2, levels = c('Astrocyte', 'OPC', 'Oligodendrocyte', 'Microglia', 'Inhibitory', 'Excitatory'))
pdf('Species_MATHYS_GeneRegulation_RNA_PUB.pdf', width = 10, height = 3)
ggscatter(resdfs, x = 'log10_round', y = 'Species', color = 'Species', size = 3, palette = c('blue', 'orange', 'darkgreen'), alpha = 0.5) +
rotate_x_text(90) + theme(text=element_text(size=20)) + ylab(expression('-log'[10]*'(FDR)')) + xlab('') + NoLegend() +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
scale_x_continuous(breaks = c(0, 2.5, 5, 7.5)) +
facet_wrap(~Var2, nrow = 1)
dev.off()




## VENN DIAGRAM ##
# Load library
library(VennDiagram)
 
#resdfs$CellType_Disease = paste0(resdfs$CellType, '_', resdfs$Var2)
resdfs2 = resdfs[resdfs$log10_round > -log10(0.05),]


# VENN
set1 <- resdfs2[resdfs2$Species == 'Human', 'CellType']
set2 <- resdfs2[resdfs2$Species == 'Chimp', 'CellType']
set3 <- resdfs2[resdfs2$Species == 'MvsHC', 'CellType']

setlist = list(set1, set2, set3)
cols = c('blue', 'orange', 'darkgreen')
nms = c('Human', 'Chimpanzee', 'MvsHC')
names(setlist) = nms

library(ggvenn)
pdf('ALZ_SINGLECELL_VENN.pdf')
ggvenn(setlist, fill_color = cols, stroke_size = 0, set_name_size = 5,text_size = 15, digits = 0, stroke_alpha = 0.5, show_percentage = F)
dev.off()




