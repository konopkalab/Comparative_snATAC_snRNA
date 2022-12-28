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
library(scran)
library(scater)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(DESeq2)
source("utility_functions.R")

####
## PREPARE DATA
####

# Take arguments from bash
args = commandArgs(trailingOnly = TRUE)
seurObjFN = args[1]
celltypeIndex = args[2]
celltypeIndex = as.numeric(celltypeIndex)

# Read the seurat object
seurObj = readRDS(seurObjFN)

seurObj$CellType = seurObj$newannot
seurObj$Sample = seurObj$orig.ident

# Load RNA
metakeep = c('CellType', 'Sample', 'Species')
seurObj@meta.data = seurObj[[metakeep]]

# Add sample metadata
smeta = import('ba23_rna_metadata.xlsx')
smeta$Species = NULL
meta = seurObj[[]]
meta = cbind(meta, smeta[match(meta$Sample, smeta$Sample_id), ])
meta$lib_batch = factor(meta$lib_batch)

# Order species so that it is Human, Chimp, Macaque
meta$Species = factor(meta$Species, levels = c("human", "chimp", "macaque"))
seurObj@meta.data = meta

# Normalize and log transform RNA
seurObj = NormalizeData(seurObj, normalization.method = "LogNormalize", scale.factor = 10000, assay = 'RNA')


####
## FIND GENES TO TEST
####

# Subset the given cell type
celltypes = unique(seurObj$CellType)
subSeur = subset(seurObj, subset = CellType %in% celltypes[celltypeIndex])

# Convert seurat to sce
DefaultAssay(subSeur) = 'RNA'
allSCE = as.SingleCellExperiment(subSeur)

# Pseudobulk SCE assay
allGroups = colData(allSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
allPseudoSCE = sumCountsAcrossCells(allSCE, allGroups)
allAggMat = allPseudoSCE@assays@data$sum
colnames(allAggMat) = colData(allPseudoSCE)[['Sample_id']]

# Find genes expressed in all 4 samples of at least one species
hExp = which(rowSums(allAggMat[, grepl('h', colnames(allAggMat))] <= 10) == 0) %>% names
cExp = which(rowSums(allAggMat[, grepl('c', colnames(allAggMat))] <= 10) == 0) %>% names
mExp = which(rowSums(allAggMat[, grepl('m', colnames(allAggMat))] <= 10) == 0) %>% names
expgns = Reduce(union, list(hExp, cExp, mExp))

####
## PAIRWISE DGE
####

# Human vs Chimp #
hcSeurat = subset(subSeur, subset = Species %in% c('human', 'chimp'))
hcSeurat$Species = droplevels(hcSeurat$Species)

# Create SCE assay
DefaultAssay(hcSeurat) = 'RNA'
hcSCE = as.SingleCellExperiment(hcSeurat)

# Pseudobulk SCE assay
hcGroups = colData(hcSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
hcPseudoSCE = sumCountsAcrossCells(hcSCE, hcGroups)
hcGroups = colData(hcPseudoSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
hcGroups$Sample_id = factor(hcGroups$Sample_id)
hcGroups$sex = factor(hcGroups$sex)
hcGroups$Species = factor(hcGroups$Species, levels = c('chimp', 'human'))
hcAggMat = hcPseudoSCE@assays@data$sum
hcAggMat = hcAggMat[expgns,]

# Run DGE
hcDGEL = DGEList(counts = hcAggMat)
hcDGEL = calcNormFactors(hcDGEL)
hcDesign = model.matrix(~Species + human_age + sex + lib_batch, data = hcGroups)
hcDGEL = estimateDisp(hcDGEL, hcDesign)

hcFit <- glmFit(hcDGEL,hcDesign)
hcLrt <- glmLRT(hcFit,coef='Specieshuman')
hcRes = topTags(hcLrt, n = nrow(hcAggMat), sort.by = 'none') %>% as.data.frame

# Human vs Macaque #
hmSeurat = subset(subSeur, subset = Species %in% c('human', 'macaque'))
hmSeurat$Species = droplevels(hmSeurat$Species)

# Create SCE assay
DefaultAssay(hmSeurat) = 'RNA'
hmSCE = as.SingleCellExperiment(hmSeurat)

# Pseudobulk SCE assay
hmGroups = colData(hmSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
hmPseudoSCE = sumCountsAcrossCells(hmSCE, hmGroups)
hmGroups = colData(hmPseudoSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
hmGroups$Sample_id = factor(hmGroups$Sample_id)
hmGroups$sex = factor(hmGroups$sex)
hmGroups$Species = factor(hmGroups$Species, levels = c('macaque', 'human'))
hmAggMat = hmPseudoSCE@assays@data$sum
hmAggMat = hmAggMat[expgns,]

# Run DGE
hmDGEL = DGEList(counts = hmAggMat)
hmDGEL = calcNormFactors(hmDGEL)
hmDesign = model.matrix(~Species + human_age + sex + lib_batch, data = hmGroups)
hmDGEL = estimateDisp(hmDGEL, hmDesign)

hmFit = glmFit(hmDGEL, hmDesign)
hmLrt = glmLRT(hmFit, coef = 'Specieshuman')
hmRes = topTags(hmLrt, n = nrow(hmAggMat), sort.by = 'none') %>% as.data.frame

# Chimp vs Macaque #
cmSeurat = subset(subSeur, subset = Species %in% c('chimp', 'macaque'))
cmSeurat$Species = droplevels(cmSeurat$Species)

# Create SCE assay
DefaultAssay(cmSeurat) = 'RNA'
cmSCE = as.SingleCellExperiment(cmSeurat)

# Pseudobulk SCE assay
cmGroups = colData(cmSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
cmPseudoSCE = sumCountsAcrossCells(cmSCE, cmGroups)
cmGroups = colData(cmPseudoSCE)[, c("Sample_id", "Species", "human_age", "sex", "lib_batch")]
cmGroups$Sample_id = factor(cmGroups$Sample_id)
cmGroups$sex = factor(cmGroups$sex)
cmGroups$Species = factor(cmGroups$Species, levels = c('macaque', 'chimp'))
cmAggMat = cmPseudoSCE@assays@data$sum
cmAggMat = cmAggMat[expgns,]

# Run DGE
cmDGEL = DGEList(counts = cmAggMat)
cmDGEL = calcNormFactors(cmDGEL)
cmDesign = model.matrix(~Species + human_age + sex + lib_batch, data = cmGroups)
cmDGEL = estimateDisp(cmDGEL, cmDesign)

cmFit = glmFit(cmDGEL, cmDesign)
cmLrt = glmLRT(cmFit, coef = 'Specieschimp')
cmRes = topTags(cmLrt, n = nrow(cmAggMat), sort.by = 'none') %>% as.data.frame

# Give unique names
colnames(hcRes) = paste0('HC_', colnames(hcRes))
colnames(hmRes) = paste0('HM_', colnames(hmRes))
colnames(cmRes) = paste0('CM_', colnames(cmRes))

# Bind all comparisons
finalRes = Reduce(cbind, list(hcRes, hmRes, cmRes)) %>% as.data.frame
finalRes$Regulation = 'NS'
finalRes$Evolution = 'NS'

# Find human and chimpanzee specific up / down regulation
finalRes[finalRes$HC_FDR < 0.05 & finalRes$HC_logFC > 0.3 &
	finalRes$HM_FDR < 0.05 & finalRes$HM_logFC > 0.3 &
	finalRes$CM_FDR > 0.1, 'Regulation'] = 'Human_UP'

finalRes[finalRes$HC_FDR < 0.05 & finalRes$HC_logFC < -0.3 &
	finalRes$HM_FDR < 0.05 & finalRes$HM_logFC < -0.3 &
	finalRes$CM_FDR > 0.1, 'Regulation'] = 'Human_DOWN'

finalRes[finalRes$HC_FDR < 0.05 & finalRes$HC_logFC < -0.3 &
	finalRes$HM_FDR > 0.1 &
	finalRes$CM_FDR < 0.05 & finalRes$CM_logFC > 0.3, 'Regulation'] = 'Chimp_UP'

finalRes[finalRes$HC_FDR < 0.05 & finalRes$HC_logFC > 0.3 &
	finalRes$HM_FDR > 0.1 &
	finalRes$CM_FDR < 0.05 & finalRes$CM_logFC < -0.3, 'Regulation'] = 'Chimp_DOWN'


# Find MvsHC genes
finalRes[finalRes$HC_FDR > 0.1 &
	finalRes$HM_FDR < 0.05 & finalRes$HM_logFC < -0.3 &
	finalRes$CM_FDR < 0.05 & finalRes$CM_logFC < -0.3, 'Regulation'] = 'Macaque_UP'

finalRes[finalRes$HC_FDR > 0.1 &
	finalRes$HM_FDR < 0.05 & finalRes$HM_logFC > 0.3 &
	finalRes$CM_FDR < 0.05 & finalRes$CM_logFC > 0.3, 'Regulation'] = 'Macaque_DOWN'

# Assign species specifically regulated genes
finalRes[finalRes$Regulation %in% c('Human_UP', 'Human_DOWN'), 'Evolution'] = 'Human_Specific'
finalRes[finalRes$Regulation %in% c('Chimp_UP', 'Chimp_DOWN'), 'Evolution'] = 'Chimp_Specific'
finalRes[finalRes$Regulation %in% c('Macaque_UP', 'Macaque_DOWN'), 'Evolution'] = 'MvsHC'

geneCtype = data.frame(Gene = rownames(finalRes), CellType = celltypes[celltypeIndex])
finalRes = cbind(geneCtype, finalRes)

saveRDS(finalRes, paste0('PseudoBulk_DEGs_', celltypes[celltypeIndex], '.RDS'))


