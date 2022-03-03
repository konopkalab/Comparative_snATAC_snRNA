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
source("utility_functions.R")


####
## PREPARE DATASETS
####

# DEGs #
alldegs = import('COMBINED_SignGenes_DEGs.xlsx')
humandeg = alldegs[alldegs$Evolution == 'Human_Specific',] %>% select(Gene,cluster) %>% rename(Class = 'cluster') %>% as.data.frame
chimpdeg = alldegs[alldegs$Evolution == 'Chimp_Specific',] %>% select(Gene,cluster) %>% rename(Class = 'cluster') %>% as.data.frame
macaquedeg = alldegs[alldegs$MvsLCA != 'NS',] %>% select(Gene,cluster) %>% rename(Class = 'cluster') %>% as.data.frame

human_list = split(humandeg, humandeg$Class)
chimp_list = split(chimpdeg, chimpdeg$Class)
macaque_list = split(macaquedeg, macaquedeg$Class)

human_list = lapply(human_list, function(x){x$Gene})
chimp_list = lapply(chimp_list, function(x){x$Gene})
macaque_list = lapply(macaque_list, function(x){x$Gene})


human_list = list(unique(humandeg$Gene))
chimp_list = list(unique(chimpdeg$Gene))
macaque_list = list(unique(macaquedeg$Gene))


# Gandal et al., 2018 #
gandal_deg = import_list("Gandal_2018/aat8127_Table_S1.xlsx")
gandal_deg = gandal_deg[[2]]

# Keep protein coding and autosomal
gandal_sign = gandal_deg[gandal_deg$gene_type == 'protein_coding' & gandal_deg$chr %in% paste0('chr', c(1:22)),]

# Reshape for enrichment
gandal_asd = gandal_sign[gandal_sign$ASD.fdr < 0.05, ] %>% select(gene_name) %>% as.data.frame %>% add_column(DEFINITION = 'ASD') %>% rename(Gene = 'gene_name')

gandal_scz = gandal_sign[gandal_sign$SCZ.fdr < 0.05, ] %>% select(gene_name) %>% as.data.frame %>% add_column(DEFINITION = 'SCZ') %>% rename(Gene = 'gene_name')

gandal_bd = gandal_sign[gandal_sign$BD.fdr < 0.05, ] %>% select(gene_name) %>% as.data.frame %>% add_column(DEFINITION = 'BD') %>% rename(Gene = 'gene_name')

gandal_df = rbind(gandal_asd, gandal_scz, gandal_bd)
gandal_list = split(gandal_df, gandal_df$DEFINITION)
gandal_list = lapply(gandal_list, function(x){x$Gene})

####
## RUN ENRICHMENTS
####

# Enrich with Gandal2018 DEGs #
library(GeneOverlap)

# Assign background
bcg = unique(c(gandal_sign$gene_name, unique(alldegs$Gene))) %>% length

# Run Enrichment
names(human_list) = paste('Human', names(human_list), sep = '_')
names(chimp_list) = paste('Chimp', names(chimp_list), sep = '_')
names(macaque_list) = paste('Macaque', names(macaque_list), sep = '_')

degs_list = c(human_list, chimp_list, macaque_list)

resgom = newGOM(degs_list, gandal_list, genome.size=bcg)
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
resdf$is_sign = ifelse(resdf$log10_FDR > 2, '<0.01',
			ifelse(resdf$log10_FDR > 1.3, '<0.05', 'NS'))
resdf$Species = factor(resdf$Species, levels = c('Human', 'Chimp', 'Macaque'))
resdf$Var2 = factor(resdf$Var2, levels = c('SCZ', 'BD', 'ASD'))

pdf('Species_GeneRegulation_Enrich_Heat.pdf', width = 10, height = 10)
ggplot(resdf, aes(Species, CellType, fill = is_sign))+
 geom_tile(color = "white") +
 scale_fill_manual(values = c('red', 'pink', 'grey')) +
 facet_wrap(~Var2) +
 rotate_x_text(90) +
  theme_minimal() +
 geom_text(aes(label = log10_round))
dev.off()

resdf$Var2 = factor(resdf$Var2, levels = c('ASD', 'SCZ', 'BD'))
pdf('Species_GeneRegulation_RNA_PUB.pdf', width = 7, height = 7)
ggscatter(resdf, y = 'log10_round', x = 'Species', color = 'Species', size = 7, palette = c('blue', 'orange', 'darkgreen')) +
rotate_x_text(90) + theme(text=element_text(size=20)) + ylab(expression('-log'[10]*'(FDR)')) + xlab('') + NoLegend() +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') + scale_y_reverse() +
facet_wrap(~Var2)
dev.off()

resdf$CellType_Disease = paste0(resdf$CellType, '_', resdf$Var2)
export(resdf, 'Bulk_Stats.xlsx')


