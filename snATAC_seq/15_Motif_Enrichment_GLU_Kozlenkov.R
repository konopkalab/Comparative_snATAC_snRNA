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
library(Signac)
library(Seurat)
library(dplyr)
library(tidyverse)
library(rio)
library(TFBSTools)
library(motifmatchr)
library(Seurat)
library(JASPAR2018)
library(liftOver)
library(dplyr)
library(Signac)
library(ggpubr)
library(data.table)
library('bedr')
source("utility_functions.R")
library("BSgenome.Hsapiens.UCSC.hg38")

####
## CHIP-SEQ to snATAC-SEQ
####

# Load Chip-seq
koz = import_list('KOZLENKOV_2020/pnas.2011884117.sd06.xlsx')
gluHUMAN = koz[[1]][,1:3]

kozDA = import_list('KOZLENKOV_2020/pnas.2011884117.sd11.xlsx')

# Up CREs for human and chimp
kozGLU_HUMAN = kozDA[[1]][,1:3] %>% add_column(id = paste0(.[,1], .[,2], .[,3])) %>% .[!duplicated(.$id),] %>% .[,1:3]
kozGLU_CHIMP = kozDA[[2]][,1:3] %>% add_column(id = paste0(.[,1], .[,2], .[,3])) %>% .[!duplicated(.$id),] %>% .[,1:3]

# Load ATAC-seq
bcg = fread('COMBINED_Allpeaks_DARs.txt') %>% as.data.frame()
bcg = bcg[bcg$major == 'Excitatory',]

####
## DATASETS FOR MOTIF ENRICHMENT
####

# Load snRNA-seq motif object
humanobj = readRDS('motif_obj_human_all.RDS')
chimpobj = readRDS('motif_obj_chimp_all.RDS')

# Load background peaks and chimp coordinates
Idents(humanobj) = humanobj$broadannot
Idents(chimpobj) = chimpobj$broadannot

# Read each species own coordinates for the peakset
chimp_coords = read.table('merged_lifted_chimp.bed')
macaque_coords = read.table('merged_lifted_macaque.bed')

#####
## ALL CREs
#####

# BED file of all snATAC-seq CREs
bcgSub = bcg
bcg_bed = unique(bcgSub$peak) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))
bcg_bed = bedr.sort.region(bcg_bed, verbose = F)

# Overlap Chip-Seq with snATAC-seq for background
atacOv = bedr(input = list(a = gluHUMAN, b = bcg_bed), method = "intersect", params = "-loj", verbose = F)
atacOv = atacOv[atacOv[,5] != '-1',]
kozATACBCG = unique(paste0(atacOv[,4], ':', atacOv[,5], '-', atacOv[,6]))

# Overlap Chip-Seq with snATAC-seq for HS-CREs and run enrichment
atacOvHS = bedr(input = list(a = kozGLU_HUMAN, b = bcg_bed), method = "intersect", params = "-loj", verbose = F)
atacOvHS = atacOvHS[atacOvHS[,5] != '-1',]
kozATACHS = unique(paste0(atacOvHS[,4], ':', atacOvHS[,5], '-', atacOvHS[,6]))
kozATACHS = kozATACHS %>% unique %>% gsub(':', '-', .)
kozATACBCG = kozATACBCG %>% unique %>% gsub(':', '-', .)
kozenrH = FindMotifs(object = humanobj, features = kozATACHS, background = kozATACBCG)
kozenrH$FDR = p.adjust(kozenrH$pvalue, method = 'BH')

rio::export(kozenrH, 'Koz_Human_UP_Motif_Enrich_GLU.xlsx')

# Overlap Chip-Seq with snATAC-seq for CS-CREs and run enrichment
atacOvCS = bedr(input = list(a = kozGLU_CHIMP, b = bcg_bed), method = "intersect", params = "-loj", verbose = F)
atacOvCS = atacOvCS[atacOvCS[,5] != '-1',]
kozATACCS = unique(paste0(atacOvCS[,4], ':', atacOvCS[,5], '-', atacOvCS[,6]))
kozATACCS = kozATACCS %>% unique %>% gsub(':', '-', .)
kozATACBCG = kozATACBCG %>% unique %>% gsub(':', '-', .)
kozATACCS_Ccoord = chimp_coords[match(kozATACCS, chimp_coords$V1), 'V5']
kozATACBCG_Ccoord = chimp_coords[match(kozATACBCG, chimp_coords$V1), 'V5']
kozenrC = FindMotifs(object = chimpobj, features = kozATACCS_Ccoord, background = kozATACBCG_Ccoord)
kozenrC$FDR = p.adjust(kozenrC$pvalue, method = 'BH')

rio::export(kozenrC, 'Koz_Chimp_UP_Motif_Enrich_GLU.xlsx')

####
## PLOT IEGs
####

# Load enrichments
kozenrH = rio::import('Koz_Human_UP_Motif_Enrich_GLU.xlsx')
kozenrC = rio::import('Koz_Chimp_UP_Motif_Enrich_GLU.xlsx')

# Load all TFs found to be enriched in HS in snATAC-seq
hmot = rio::import('All_Human_UP_Motif_Enrich.xlsx')
hmotNAMES = hmot[hmot$cluster %in% names(table(bcg$cluster)), 'motif.name']

# Constitutively active TFs upregulating IEG-TFs
hmotNAMES1 = hmotNAMES[grepl('MEF2|CREB|SRF', hmotNAMES)] %>% unique

# IEG TFs for excitatory neurons
iegsH = readRDS('hrvatin_2018/IEGs_EXC_HUMAN.RDS')
iegsH_check = paste(iegsH, collapse = '|')
hmotNAMES2 = hmotNAMES[grep(iegsH_check, hmotNAMES)] %>% unique

hmotNAMES = c(hmotNAMES1, hmotNAMES2)


# Plot
hIEG = kozenrH[kozenrH$motif.name %in% hmotNAMES,] %>% add_column(Species = 'Human')
cIEG = kozenrC[kozenrC$motif.name %in% hmotNAMES,] %>% add_column(Species = 'Chimp')
toplot = rbind(hIEG, cIEG)
toplot$log10FDR = -log10(toplot$FDR)

hIEG = kozenrH[kozenrH$motif.name %in% hmotNAMES,] %>% add_column(Species = 'Human')
cIEG = kozenrC[kozenrC$motif.name %in% hmotNAMES,] %>% add_column(Species = 'Chimp')
toplot = rbind(hIEG, cIEG)
toplot$log10FDR = -log10(toplot$FDR)

toplot[toplot$motif.name %in% hmotNAMES1, 'type'] = 'CONSTITUTIVE_TFs'
toplot[toplot$motif.name %in% hmotNAMES2, 'type'] = 'IEG_TFs'
toplot$Species = factor(toplot$Species, levels = c('Human', 'Chimp'))

toplotIEG = toplot[toplot$type != 'CONSTITUTIVE_TFs',]
toplotCONST = toplot[toplot$type == 'CONSTITUTIVE_TFs',]

library(ggrepel)
pdf('IEG_GLU_KOZLENKOV_ENRICH_CONST_AND_IEG.pdf', width = 12, height = 5)
ggscatter(toplot, x = 'log10FDR', y = 'fold.enrichment', size = 5, color = 'grey',  xlim = c(0,15), ylim = c(0,2), alpha = 0.3) +
theme(text=element_text(size=30, face = 'bold')) + ylab('Fold Enrichment') +
xlab(expression('-log'[10]*'(FDR)')) +
geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'red') + NoLegend() +
geom_text_repel(data = toplotCONST[order(toplotCONST$log10FDR,decreasing=T),][1:5,], aes(label = motif.name), nudge_y = 0.02, nudge_x = -0.05, fontface = 'bold', size = 5) +
geom_point(data = toplotCONST, aes(label = motif.name, color = Species), size = 3) +
scale_colour_manual(values = c('blue', 'orange')) +
facet_wrap(~Species)
dev.off()


# Plot conserved TFs
cmns = readRDS('OpenRegionTfsCommon.RDS')

hCONS = kozenrH[kozenrH$motif.name %in% cmns,] %>% add_column(Species = 'Human')
cCONS = kozenrC[kozenrC$motif.name %in% cmns,] %>% add_column(Species = 'Chimp')
toplot = rbind(hCONS, cCONS)
toplot$log10FDR = -log10(toplot$FDR)

pdf('CONS_GLU_KOZLENKOV_ENRICH.pdf', width = 10, height = 7)
ggscatter(toplot, x = 'log10FDR', y = 'fold.enrichment', size = 5, color = 'Species', palette = c('orange', 'blue'), alpha = 0.5, xlim = c(0,10), ylim = c(0,1.5)) +
theme(text=element_text(size=30, face = 'bold')) + ylab('Fold Enrichment') +
xlab(expression('-log'[10]*'(FDR)')) +
geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'red') + NoLegend() +
facet_wrap(~Species)
dev.off()

