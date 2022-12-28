rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(Signac)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(feather)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

####
## Load Data
####

# Human -- SNARE
hum = read.csv('Analysis_lein_bdbag_2020_08_11/data/Multimodal/sncell/SNARE/human/processed/counts/counts/M1/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt', sep = '\t', comment.char = "", skip = 0, check.names = FALSE, header = T)

humAll = hum

# Marmoset -- SNARE
mar = read.csv('Analysis_lein_bdbag_2020_08_11/data/Multimodal/sncell/SNARE/marmoset/processed/counts/counts/M1/Zhang_BICCN-H_20190730_20190903_marMOp_Final_Sample_Metadata.txt', sep = '\t', comment.char = "", skip = 0, check.names = FALSE, header = T)

marAll = mar

# Human -- Transcriptomics
humTr = read.csv('Analysis_lein_bdbag_2020_08_11/data/Transcriptomics/sncell/10X/human/processed/counts/M1/Human_MTG_SSv4_Samples_Columns/Human_MTG_SSv4_Samples_Columns.csv')

# Macaque -- Transcriptomics
macTr = read_feather('Analysis_lein_bdbag_2020_08_11/data/Transcriptomics/sncell/10X/macaque/processed/counts/counts/M1/Macaque_L5_Dissection_M1_10xV3_Metadata.feather')
macTr = as.data.frame(macTr)

# Marmoset -- Transcriptomics
marTr = read_feather('Analysis_lein_bdbag_2020_08_11/data/Transcriptomics/sncell/10X/marmoset/processed/counts/counts/M1/Marmoset_M1_10xV3_Metadata.feather')
marTr = as.data.frame(marTr)


####
## OPC OR MOL TO ALL
####

# HUMAN -- SNARE
humAll$CellType = humAll$AC_cluster_label
humAll = humAll[grepl('Oli|OPC|Mic|Ast|Glia', humAll$CellType),] # Keep only glia

humAll[grepl('Oligo', humAll$CellType), 'CellType'] = 'MOL'
humAll[grepl('OPC', humAll$CellType), 'CellType'] = 'OPC'

humSnare = humAll %>% group_by(patient, CellType) %>% group_keys %>% as.data.frame
humSnare$size = humAll %>% group_by(patient, CellType) %>% group_size
humSnareTotSize = table(humAll$patient) %>% as.data.frame %>% dplyr::rename(patient='Var1',totsize='Freq')
humSnare2 = merge(humSnare, humSnareTotSize)
humSnare2$ratio = humSnare2$size / humSnare2$totsize

humSnare2[humSnare2$CellType == 'OPC', 'ratio'] %>% summary
humSnare2[humSnare2$CellType == 'MOL', 'ratio'] %>% summary


# MARMOSET
marAll$CellType = marAll$AC_cluster_label
marAll = marAll[grepl('Oli|OPC|Mic|Ast|Glia', marAll$CellType),] # Keep only glia

marAll[grepl('Oligo', marAll$CellType), 'CellType'] = 'MOL'
marAll[grepl('OPC', marAll$CellType), 'CellType'] = 'OPC'

marSnare = marAll %>% group_by(patient, CellType) %>% group_keys %>% as.data.frame
marSnare$size = marAll %>% group_by(patient, CellType) %>% group_size
marSnareTotSize = table(marAll$patient) %>% as.data.frame %>% dplyr::rename(patient='Var1',totsize='Freq')
marSnare2 = merge(marSnare, marSnareTotSize)
marSnare2$ratio = marSnare2$size / marSnare2$totsize

marSnare2[marSnare2$CellType == 'OPC', 'ratio'] %>% summary
marSnare2[marSnare2$CellType == 'MOL', 'ratio'] %>% summary


# HUMAN -- Transcriptome
humTr$CellType = humTr$cluster
humTr = humTr[humTr$donor %in% names(which(table(humTr$donor) > 1000)),]
humTr = humTr[grepl('Oli|OPC|Mic|Ast|Glia', humTr$CellType),] # Keep only glia
humTr = humTr[!(grepl('FLT1', humTr$CellType)),] # Remove endothelia

humTr[grepl('Oligo', humTr$cluster), 'CellType'] = 'MOL'
humTr[grepl('OPC', humTr$cluster), 'CellType'] = 'OPC'

humTrans = humTr %>% group_by(donor, CellType) %>% group_keys %>% as.data.frame
humTrans$size = humTr %>% group_by(donor, CellType) %>% group_size

humTransTotSize = table(humTr$donor) %>% as.data.frame %>% dplyr::rename(donor='Var1',totsize='Freq')

humTrans2 = merge(humTrans, humTransTotSize, by = 'donor')
humTrans2$ratio = humTrans2$size / humTrans2$totsize

humTrans2[humTrans2$CellType == 'OPC', 'ratio'] %>% summary
humTrans2[humTrans2$CellType == 'MOL', 'ratio'] %>% summary

# MACAQUE -- Transcriptome
macTr$CellType = macTr$cluster_label
macTr$donor = macTr$Donor_label
macTr = macTr[grepl('Oli|OPC|Mic|Ast|Glia', macTr$CellType),] # Keep only glia
macTr = macTr[!(grepl('FLT1', macTr$CellType)),] # Remove endothelia
macTr[grepl('OPALIN', macTr$CellType), 'CellType'] = 'MOL'
macTr[grepl('PDGFRA', macTr$CellType), 'CellType'] = 'OPC'

macTrans = macTr %>% group_by(donor, CellType) %>% group_keys %>% as.data.frame
macTrans$size = macTr %>% group_by(donor, CellType) %>% group_size

macTransTotSize = table(macTr$donor) %>% as.data.frame %>% dplyr::rename(donor='Var1',totsize='Freq')

macTrans2 = merge(macTrans, macTransTotSize, by = 'donor')
macTrans2$ratio = macTrans2$size / macTrans2$totsize

macTrans2[macTrans2$CellType == 'OPC', 'ratio'] %>% summary
macTrans2[macTrans2$CellType == 'MOL', 'ratio'] %>% summary


# MARMOSET -- Transcriptome
marTr$CellType = marTr$cluster_label
marTr$donor = marTr$donor_id
marTr = marTr[grepl('Oli|OPC|Mic|Ast|Glia', marTr$CellType),] # Keep only glia
marTr = marTr[!(grepl('FLT1', marTr$CellType)),] # Remove endothelia
marTr[grepl('Oligo', marTr$CellType), 'CellType'] = 'MOL'
marTr[grepl('OPC', marTr$CellType), 'CellType'] = 'OPC'

marTrans = marTr %>% group_by(donor, CellType) %>% group_keys %>% as.data.frame
marTrans$size = marTr %>% group_by(donor, CellType) %>% group_size

marTransTotSize = table(marTr$donor) %>% as.data.frame %>% dplyr::rename(donor='Var1',totsize='Freq')

marTrans2 = merge(marTrans, marTransTotSize, by = 'donor')
marTrans2$ratio = marTrans2$size / marTrans2$totsize

marTrans2[marTrans2$CellType == 'OPC', 'ratio'] %>% summary
marTrans2[marTrans2$CellType == 'MOL', 'ratio'] %>% summary



# Plot Transcriptome
humOPC_Tr = humTrans2[humTrans2$CellType == 'OPC', 'ratio']
macOPC_Tr = macTrans2[macTrans2$CellType == 'OPC', 'ratio']
marOPC_Tr = marTrans2[marTrans2$CellType == 'OPC', 'ratio']

humMOL_Tr = humTrans2[humTrans2$CellType == 'MOL', 'ratio']
macMOL_Tr = macTrans2[macTrans2$CellType == 'MOL', 'ratio']
marMOL_Tr = marTrans2[marTrans2$CellType == 'MOL', 'ratio']


trToplotOPC = data.frame(vals = c(humOPC_Tr, macOPC_Tr, marOPC_Tr),
			species = c(rep('Human', length(humOPC_Tr)),
					rep('Macaque', length(macOPC_Tr)),
					rep('Marmoset', length(marOPC_Tr))))

trToplotMOL = data.frame(vals = c(humMOL_Tr, macMOL_Tr, marMOL_Tr),
			species = c(rep('Human', length(humMOL_Tr)),
					rep('Macaque', length(macMOL_Tr)),
					rep('Marmoset', length(marMOL_Tr))))

trToplotOPC$CellType = 'OPC'
trToplotMOL$CellType = 'MOL'
combPlot = rbind(trToplotOPC, trToplotMOL)
combPlot$CellType = factor(combPlot$CellType, levels = c('OPC', 'MOL'))
comps = list(c('Human', 'Macaque'), c('Human', 'Marmoset'))

pdf('Bakken_Ratios_Transcriptome.pdf', width = 7)
ggboxplot(combPlot, x = 'species', y = 'vals', color = 'species', add = 'dotplot', ylim = c(0,1.1)) +
ylab('Proportion_To_Glia') +
xlab('') +
theme(text=element_text(size=20), legend.position = 'right') +
rotate_x_text(90) +
stat_compare_means(comparisons = comps, method = 'wilcox.test') +
facet_wrap(~CellType)
dev.off()


# Plot Snare AND Transcriptome
combPlotTR = combPlot

humOPC_Sn = humSnare2[humSnare2$CellType == 'OPC', 'ratio']
marOPC_Sn = marSnare2[marSnare2$CellType == 'OPC', 'ratio']

humMOL_Sn = humSnare2[humSnare2$CellType == 'MOL', 'ratio']
marMOL_Sn = marSnare2[marSnare2$CellType == 'MOL', 'ratio']


trToplotOPC = data.frame(vals = c(humOPC_Sn, marOPC_Sn),
			species = c(rep('Human', length(humOPC_Sn)),
					rep('Marmoset', length(marOPC_Sn))))

trToplotMOL = data.frame(vals = c(humMOL_Sn, marMOL_Sn),
			species = c(rep('Human', length(humMOL_Sn)),
					rep('Marmoset', length(marMOL_Sn))))

trToplotOPC$CellType = 'OPC'
trToplotMOL$CellType = 'MOL'
combPlot = do.call(rbind, list(trToplotOPC, trToplotMOL, combPlotTR))
combPlot$CellType = factor(combPlot$CellType, levels = c('OPC', 'MOL'))
combPlot$species = factor(combPlot$species, levels = c('Human', 'Macaque', 'Marmoset'))
comps = list(c('Human', 'Marmoset'), c('Human', 'Macaque'))

pdf('Bakken_Ratios_Snare_AND_Transcriptome_Human_vs_Marmoset_Macaque.pdf', width = 6)
ggboxplot(combPlot, x = 'species', y = 'vals', color = 'species', add = 'dotplot', ylim = c(0,1.1)) +
ylab('Proportion_To_Glia') +
xlab('') +
theme(text=element_text(size=20), legend.position = 'right') +
rotate_x_text(90) +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = c(0.9, 1)) +
facet_wrap(~CellType)
dev.off()






