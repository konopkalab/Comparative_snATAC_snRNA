rm(list = ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(Matrix)
library(ggplot2)
library(plyr)
library(Matrix.utils)
library(stringi)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(harmony)
source("utility_functions.R")
set.seed(1234)

####
## 1- PROCESS ATAC PER SPECIES
####

human_annot = readRDS("human_annotated.RDS")

# Extract excitatory neurons and rerun UMAP with same components
human_exc = subset(human_annot, subset = broadannot %in% "Exc")
DefaultAssay(human_exc) = 'peaks'
human_exc = RunTFIDF(human_exc)
human_exc = FindTopFeatures(human_exc, min.cutoff = 'q50')
human_exc <- RunSVD(object = human_exc, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
human_exc = RunUMAP(object = human_exc, dims = 1:30, reduction = 'lsi')

pdf('human_Exc_Atac_BeforeHarmony.pdf')
DimPlot(human_exc, label = T, group.by = 'Sample')
dev.off()

human_exc <- RunHarmony(object = human_exc, group.by.vars = 'Sample', reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
human_exc = RunUMAP(object = human_exc, dims = 1:30, reduction = 'harmony')

pdf('human_Exc_Atac_AfterHarmony.pdf')
DimPlot(human_exc, label = T, group.by = 'Sample')
dev.off()

human_exc <- FindNeighbors(object = human_exc, reduction = 'harmony', dims = 1:30)
DefaultAssay(human_exc) = 'peaks'
human_exc <- FindClusters(object = human_exc, verbose = FALSE, resolution = 1)

pdf('human_Exc_Atac_Clusters.pdf')
DimPlot(human_exc, label = T)
dev.off()

pdf('human_Exc_Atac_Clusters_Depth.pdf')
ggboxplot(human_exc[[]], x = 'seurat_clusters', y=  'rip')
dev.off()


pdf('human_clusters_GLIAMARKERS.pdf', width = 20, height = 8)
VlnPlot(human_exc, c('ST18', 'MOG', 'AQP4', 'APBB1IP', 'PTPRZ1'), assay = 'RNA', pt.size = 0)
dev.off()

saveRDS(human_exc, "~/workdir/pr3/atacseq/human_chimp_macaque/11_exc_integration/human_to_human/human_reclustered.RDS")


####
## 2- LABEL TRANSFER
####

# Take human RNA
allrna  = readRDS('EXC_Integrated_Annotated.RDS')
rna_exc_seur = subset(allrna, subset = Species == 'human')

# Check common features for transfer
rna_exc_seur = SCTransform(rna_exc_seur, ncells = 1000)
rna_exc_seur = FindVariableFeatures(rna_exc_seur)
feats = intersect(VariableFeatures(rna_exc_seur), rownames(human_exc@assays$SCT@data))
length(feats)
Idents(rna_exc_seur) = rna_exc_seur$newannot

# Find anchors
transfer_anchors <- FindTransferAnchors(reference = rna_exc_seur, query = human_exc,
					dims = 1:20,
					features = feats,
					reference.assay = "SCT",
					query.assay = "SCT",
					reduction = "cca")


# Make cell type predictions based on anchors
pred_sub <- TransferData(anchorset = transfer_anchors,
			refdata = rna_exc_seur$newannot, 
			weight.reduction = human_exc[["harmony"]])

# Add predictions to meta data
human_exc$pred_sub = pred_sub$predicted.id
human_exc$pred_sub_score = pred_sub$prediction.score.max
pdf("Human_pred_EXC_AllenRef_SUB.pdf")
hist(human_exc$pred_sub_score)
dev.off()

# Plot label transfer
# Make the colors match
human_exc$pred_sub <- factor(human_exc$pred_sub, levels = levels(rna_exc_seur$newannot))
Idents(human_exc) = human_exc$pred_sub
transferPlot(human_exc, rna_exc_seur, 'LABEL_TRANSFER_human_Rna_Atac_SUB')

saveRDS(human_exc, "Exc_human_annotated.RDS")

# Cluster at higher resolution if needed
human_exc$predicted.id = human_exc$pred_sub
human_exc <- FindClusters(object = human_exc, verbose = FALSE, resolution = 1.2)
pdf('human_Exc_Atac_Clusters_res1.2.pdf')
DimPlot(human_exc, label = T)
dev.off()

ltpred(human_exc, 'HEATMAP_LABEL_TRANSFER_human_Rna_Atac_SUB_res1.2')

####
## 03 - FINAL ANNOTATION
####


mapnames = setNames(c("L2-3_1", "L3-5_RORB_1", "Mixed", "L5-6_THEMIS_1", "L4-5_RORB_2-L4-6_RORB_1", "L4-5_RORB_1",
			 "L2-3_2", "L5-6_FEZF2_1", "L2-3_3", "L3-5_RORB_2", "L4-6_RORB_2",
			 "L3-5_RORB_3", "L5-6_THEMIS_2", "L5-6_FEZF2_2-3", "L4-6_FEZF2", "L5-6_FEZF2_2-3"),
		      c(0:15))

human_exc[["newannot"]] <- mapnames[as.character(human_exc$seurat_clusters)]
human_exc$newannot = factor(human_exc$newannot)

pdf('Human_Atac_Exc_annotated_FINAL.pdf')
DimPlot(human_exc, group.by = 'newannot', label = T) + NoLegend()
dev.off()

saveRDS(human_exc, "~/workdir/pr3/atacseq/human_chimp_macaque/11_exc_integration/human_to_human/human_annotated.RDS")




