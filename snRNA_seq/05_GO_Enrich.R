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
source("utility_functions.R")

####
## PREPARE DATA
####

ourOL = readRDS('PseudoBulk_DEGs_OL.RDS')
ourOPC = readRDS('PseudoBulk_DEGs_OPC.RDS')

####
## RUN GENE ONTOLOGY
####

# Divide into up and down
hsUPOPC = ourOPC[ourOPC$Regulation == 'Human_UP', 'Gene']
hsDOWNOPC = ourOPC[ourOPC$Regulation == 'Human_DOWN', 'Gene']
bcg = rownames(ourOPC)

# Run Enrichment HUMAN UP
res = GOenrich(hsUPOPC, bcg)
res$EvoType = 'Human_UP'
res$CellType = 'OPC'
resUP = res

# Run Enrichment HUMAN DOWN
res = GOenrich(hsDOWNOPC, bcg)
res$EvoType = 'Human_DOWN'
res$CellType = 'OPC'
resDOWN = res

#allResOPC = do.call(rbind, list(resUP, resDOWN))
allResOPC = resDOWN

# COMBINE ALL
rnadf = allResOPC
rnadf$ObsOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$ObsAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'GeneRatio']) %>% as.numeric()})
rnadf$BcgOv = sapply(1:nrow(rnadf), function(x){gsub('/.*', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$BcgAll = sapply(1:nrow(rnadf), function(x){gsub('.*/', '', rnadf[x, 'BgRatio']) %>% as.numeric()})
rnadf$OddsRatio = (rnadf$ObsOv / rnadf$ObsAll) / (rnadf$BcgOv / rnadf$BcgAll)
export(rnadf, "GO_Enrich_OPC_DEGs.xlsx")

