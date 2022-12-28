rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(Signac)
library(patchwork)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(rio)
library(scater)
library(scran)
library(scater)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library('bedr')
source("utility_functions.R")

####
## DAR ANALYSIS EXCITATORY
####

# Load ATAC-seq data
human = readRDS("human_annotated.RDS")
chimp = readRDS("chimp_annotated.RDS")
macaque = readRDS("macaque_annotated.RDS")

# Merge across species
merged = Reduce(merge, list(human, chimp, macaque))
merged$tot_acc_cicero = log2(merged$nCount_RNA)

mat = merged@assays$peaks@counts
meta = merged[[]]

# Add sample metadata
smeta = import('atac_metadata.xlsx')
smeta$Species = NULL
meta = cbind(meta, smeta[match(meta$Sample, smeta$Sample_id), ])
meta$lib_batch = factor(meta$lib_batch)

# Order species so that it is Human, Chimp, Macaque
meta$Species = factor(meta$Species, levels = c("Human", "Chimp", "Macaque"))
merged@meta.data = meta

saveRDS(merged, 'merged_atac.RDS')

####
## PSEUDOBULK DGE PER CELL TYPE
####
ctypes = unique(merged$newannot)
print(ctypes)

finalResL = list()
for(i in 1:length(ctypes)){

	# Subset the cell type
	subSeur = subset(merged, subset = newannot == ctypes[i])
	print(unique(subSeur$newannot))

	# Convert seurat to sce
	DefaultAssay(subSeur) = 'peaks'
	allSCE = as.SingleCellExperiment(subSeur)

	# Pseudobulk SCE assay
	allGroups = colData(allSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	allPseudoSCE = sumCountsAcrossCells(allSCE, allGroups)
	allAggMat = allPseudoSCE@assays@data$sum
	colnames(allAggMat) = colData(allPseudoSCE)[['Sample_id']]

	# Find CREs accessible in all 4 samples of at least one species
	hExp = which(rowSums(allAggMat[, grepl('h', colnames(allAggMat))] <= 3) == 0) %>% names
	cExp = which(rowSums(allAggMat[, grepl('c', colnames(allAggMat))] <= 3) == 0) %>% names
	mExp = which(rowSums(allAggMat[, grepl('m', colnames(allAggMat))] <= 3) == 0) %>% names
	expgns1 = Reduce(union, list(hExp, cExp, mExp))

	# Top 100k CREs
	top100k = sort(rowSums(allAggMat)) %>% tail(100000) %>% names

	# Keep either top 100k CREs or the accessible CREs
	if(length(expgns1) > length(top100k)){expgns=top100k}else{expgns=expgns1}
	print(length(expgns))

	# Human vs Chimp #
	hcSeurat = subset(subSeur, subset = Species %in% c('Human', 'Chimp'))
	hcSeurat$Species = droplevels(hcSeurat$Species)

	# Create SCE assay
	DefaultAssay(hcSeurat) = 'peaks'
	hcSCE = as.SingleCellExperiment(hcSeurat)

	# Pseudobulk SCE assay
	hcGroups = colData(hcSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	hcPseudoSCE = sumCountsAcrossCells(hcSCE, hcGroups)
	hcGroups = colData(hcPseudoSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	hcGroups$Sample_id = factor(hcGroups$Sample_id)
	hcGroups$sex = factor(hcGroups$sex)
	hcGroups$Species = factor(hcGroups$Species, levels = c('Chimp', 'Human'))
	hcAggMat = hcPseudoSCE@assays@data$sum
	hcAggMat = hcAggMat[expgns,]

	# Run differential accessibility
	hcDGEL = DGEList(counts = hcAggMat)
	hcDGEL = calcNormFactors(hcDGEL)
	hcDesign = model.matrix(~Species + human_age + sex, data = hcGroups)
	hcDGEL = estimateDisp(hcDGEL, hcDesign)

	hcFit = glmFit(hcDGEL, hcDesign)
	hcLrt = glmLRT(hcFit, coef='SpeciesHuman')
	hcRes = topTags(hcLrt, n = nrow(hcAggMat), sort.by = 'none') %>% as.data.frame

	# Human vs Macaque #
	hmSeurat = subset(subSeur, subset = Species %in% c('Human', 'Macaque'))
	hmSeurat$Species = droplevels(hmSeurat$Species)

	# Create SCE assay
	DefaultAssay(hmSeurat) = 'peaks'
	hmSCE = as.SingleCellExperiment(hmSeurat)

	# Pseudobulk SCE assay
	hmGroups = colData(hmSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	hmPseudoSCE = sumCountsAcrossCells(hmSCE, hmGroups)
	hmGroups = colData(hmPseudoSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	hmGroups$Sample_id = factor(hmGroups$Sample_id)
	hmGroups$sex = factor(hmGroups$sex)
	hmGroups$Species = factor(hmGroups$Species, levels = c('Macaque', 'Human'))
	hmAggMat = hmPseudoSCE@assays@data$sum
	hmAggMat = hmAggMat[expgns,]

	# Run differential accessibility
	hmDGEL = DGEList(counts = hmAggMat)
	hmDGEL = calcNormFactors(hmDGEL)
	hmDesign = model.matrix(~Species + human_age + sex, data = hmGroups)
	hmDGEL = estimateDisp(hmDGEL, hmDesign)

	hmFit = glmFit(hmDGEL, hmDesign)
	hmLrt = glmLRT(hmFit,coef='SpeciesHuman')
	hmRes = topTags(hmLrt, n = nrow(hmAggMat), sort.by = 'none') %>% as.data.frame

	# Chimp vs Macaque #
	cmSeurat = subset(subSeur, subset = Species %in% c('Chimp', 'Macaque'))
	cmSeurat$Species = droplevels(cmSeurat$Species)

	# Create SCE assay
	DefaultAssay(cmSeurat) = 'peaks'
	cmSCE = as.SingleCellExperiment(cmSeurat)

	# Pseudobulk SCE assay
	cmGroups = colData(cmSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	cmPseudoSCE = sumCountsAcrossCells(cmSCE, cmGroups)
	cmGroups = colData(cmPseudoSCE)[, c("Sample_id", "Species", "human_age", "sex")]
	cmGroups$Sample_id = factor(cmGroups$Sample_id)
	cmGroups$sex = factor(cmGroups$sex)
	cmGroups$Species = factor(cmGroups$Species, levels = c('Macaque', 'Chimp'))
	cmAggMat = cmPseudoSCE@assays@data$sum
	cmAggMat = cmAggMat[expgns,]

	# Run differential accessibility
	cmDGEL = DGEList(counts = cmAggMat)
	cmDGEL = calcNormFactors(cmDGEL)
	cmDesign = model.matrix(~Species + human_age + sex, data = cmGroups)
	cmDGEL = estimateDisp(cmDGEL, cmDesign)

	cmFit = glmFit(cmDGEL, cmDesign)
	cmLrt = glmLRT(cmFit,coef='SpeciesChimp')
	cmRes = topTags(cmLrt, n = nrow(cmAggMat), sort.by = 'none') %>% as.data.frame

	# Give unique names
	colnames(hcRes) = paste0('HC_', colnames(hcRes))
	colnames(hmRes) = paste0('HM_', colnames(hmRes))
	colnames(cmRes) = paste0('CM_', colnames(cmRes))

	# Bind all comparisons
	finalRes = Reduce(cbind, list(hcRes, hmRes, cmRes)) %>% as.data.frame
	finalRes$Regulation = 'NS'
	finalRes$Evolution = 'NS'

	# Find human and chimpanzee specific changes
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


	# Find Macaque versus Human and Macaque versus Chimp changes
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

	geneCtype = data.frame(Gene = rownames(finalRes), CellType = ctypes[i])
	finalRes = cbind(geneCtype, finalRes)
	
	finalResL[[i]] = finalRes
	print(i)

	saveRDS(finalRes, paste0('PSEUDOBULK_', ctypes[i], '_DARs.RDS'))
}


finalResDF = do.call(rbind, finalResL)
saveRDS(finalResDF, 'PSEUDOBULK_DARs.RDS')





