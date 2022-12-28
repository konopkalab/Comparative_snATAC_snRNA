rm(list = ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")


####
## 01- INTEGRATE ALL SPECIES
####

# Read and extract excitatory cells
humanseur = readRDS('Human_Integrated_Filtered2.RDS')
humanseur = subset(humanseur, subset = broadannot == 'Exc')

chimpseur = readRDS('Chimp_Integrated_Filtered2.RDS')
chimpseur = subset(chimpseur, subset = broadannot == 'Exc')

macaqueseur = readRDS('Macaque_Integrated_Filtered2.RDS')
macaqueseur = subset(macaqueseur, subset = broadannot == 'Exc')

# Merge
seurmerg = Reduce(merge, list(humanseur, chimpseur, macaqueseur))

# Plot before integration
options(future.globals.maxSize = 4000 * 1024^2)
seurmerg = SCTransform(seurmerg, ncells = 1000)
seurmerg <- RunPCA(seurmerg, verbose = F)
seurmerg <- RunUMAP(seurmerg, dims = 1:30)

pdf("EXC_preinteg_umap.pdf")
DimPlot(seurmerg, label = TRUE) + NoLegend()
dev.off()

pdf("EXC_preinteg_species.pdf")
DimPlot(seurmerg, group.by = 'Species')
dev.off()

pdf("EXC_preinteg_samples.pdf")
DimPlot(seurmerg, group.by = 'orig.ident', split.by = 'Species')
dev.off()

# Perform SCT
seurmerg_list = SplitObject(seurmerg, split.by = "orig.ident")
for (i in 1:length(seurmerg_list)) {
	seurmerg_list[[i]] <- SCTransform(seurmerg_list[[i]], verbose = T, ncells = 1000)
	print(i)
}


# Integration
seurmerg_feats <- SelectIntegrationFeatures(object.list = seurmerg_list, nfeatures = 2000)
commonfeats = seurmerg_feats

# Integrate and cluster the data
seurmerg_list <- PrepSCTIntegration(object.list = seurmerg_list, anchor.features = commonfeats, verbose = FALSE)
seurmerg_anchors <- FindIntegrationAnchors(object.list = seurmerg_list, normalization.method = "SCT", anchor.features = commonfeats, verbose = FALSE)
seurmerg_integrated <- IntegrateData(anchorset = seurmerg_anchors, normalization.method = "SCT", verbose = FALSE)
seurmerg_integrated <- RunPCA(seurmerg_integrated, verbose = FALSE)
seurmerg_integrated <- RunUMAP(seurmerg_integrated, dims = 1:30)

pdf("EXC_species_integrated.pdf")
DimPlot(seurmerg_integrated, group.by = 'Species')
dev.off()

seurmerg_integrated <- FindNeighbors(seurmerg_integrated, dims = 1:30)
seurmerg_integrated <- FindClusters(seurmerg_integrated, resolution = 0.5)

pdf("EXC_integrated_clusters.pdf", width = 10, height = 10)
DimPlot(seurmerg_integrated, label = T) + NoLegend()
dev.off()

stackedbarplot(seurmerg_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', 'Sample_Stacked_EXC')
stackedbarplot(seurmerg_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', 'Species_Stacked_EXC')

# Plot marker genes, number of UMI to determine potential doublets
marks = c('PCDH15', 'MOG', 'BCAS1', 'AQP4', 'GFAP', 'GAD1', 'SLC17A7', 'APBB1IP', 'FLT1', 'RBFOX3')
DefaultAssay(seurmerg_integrated) = 'RNA'
seurmerg_integrated = NormalizeData(seurmerg_integrated)

pdf("EXC_markers_vlnplot.pdf", width = 24, height = 10)
VlnPlot(seurmerg_integrated, features = marks, group.by = 'seurat_clusters', pt.size = 0)
dev.off()

pdf("EXC_markers_featureplot.pdf", width = 16, height = 14)
FeaturePlot(seurmerg_integrated, features = marks, sort = T, raster = T)
dev.off()

pdf("EXC_integrated_depth.pdf")
ggboxplot(seurmerg_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()

saveRDS(seurmerg_integrated, 'EXC_Integrated.RDS')

####
## 02- FILTER INTEGRATED
####

# Filter and redo the dimensionality reduction
seurmerg_integrated = subset(seurmerg_integrated, subset = seurat_clusters %in% c(10,13,19:21), invert = T)

DefaultAssay(seurmerg_integrated) = 'integrated'
seurmerg_integrated <- RunPCA(seurmerg_integrated, verbose = FALSE)
seurmerg_integrated <- RunUMAP(seurmerg_integrated, dims = 1:30)

pdf("EXC_species_integrated.pdf")
DimPlot(seurmerg_integrated, group.by = 'Species')
dev.off()

seurmerg_integrated <- FindNeighbors(seurmerg_integrated, dims = 1:30)
seurmerg_integrated <- FindClusters(seurmerg_integrated, resolution = 0.5)

pdf("EXC_integrated_clusters.pdf", width = 10, height = 10)
DimPlot(seurmerg_integrated, label = T) + NoLegend()
dev.off()

stackedbarplot(seurmerg_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', 'Sample_Stacked_EXC')
stackedbarplot(seurmerg_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', 'Species_Stacked_EXC')

# Plot marker genes, number of UMI to determine potential doublets
marks = c('PCDH15', 'MOG', 'BCAS1', 'AQP4', 'GFAP', 'GAD1', 'SLC17A7', 'APBB1IP', 'FLT1', 'RBFOX3')

DefaultAssay(seurmerg_integrated) = 'RNA'
seurmerg_integrated = NormalizeData(seurmerg_integrated)

pdf("EXC_markers_vlnplot.pdf", width = 24, height = 10)
VlnPlot(seurmerg_integrated, features = marks, group.by = 'seurat_clusters', pt.size = 0)
dev.off()

pdf("EXC_markers_featureplot.pdf", width = 16, height = 14)
FeaturePlot(seurmerg_integrated, features = marks, sort = T, raster = T)
dev.off()

pdf("EXC_integrated_depth.pdf")
ggboxplot(seurmerg_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()

saveRDS(seurmerg_integrated, 'EXC_Integrated_Annotated.RDS')

####
## 02- FILTER INTEGRATED - 2
####

# Filter potential doublets and redo the dimensionality reduction
seurmerg_integrated = subset(seurmerg_integrated, subset = seurat_clusters %in% c(16,17), invert = T)

DefaultAssay(seurmerg_integrated) = 'integrated'
seurmerg_integrated <- RunPCA(seurmerg_integrated, verbose = FALSE)
seurmerg_integrated <- RunUMAP(seurmerg_integrated, dims = 1:30)

pdf("EXC_species_integrated.pdf")
DimPlot(seurmerg_integrated, group.by = 'Species')
dev.off()

seurmerg_integrated <- FindNeighbors(seurmerg_integrated, dims = 1:30)
seurmerg_integrated <- FindClusters(seurmerg_integrated, resolution = 1)

pdf("EXC_integrated_clusters.pdf", width = 10, height = 10)
DimPlot(seurmerg_integrated, label = T) + NoLegend()
dev.off()

stackedbarplot(seurmerg_integrated[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', 'Sample_Stacked_EXC')
stackedbarplot(seurmerg_integrated[[]], groupx = 'seurat_clusters', groupfill = 'Species', 'Species_Stacked_EXC')

# Plot marker genes, number of UMI to determine potential doublets
marks = c('PCDH15', 'MOG', 'BCAS1', 'AQP4', 'GFAP', 'GAD1', 'SLC17A7', 'APBB1IP', 'FLT1', 'RBFOX3')

DefaultAssay(seurmerg_integrated) = 'RNA'
seurmerg_integrated = NormalizeData(seurmerg_integrated)

pdf("EXC_markers_vlnplot.pdf", width = 24, height = 10)
VlnPlot(seurmerg_integrated, features = marks, group.by = 'seurat_clusters', pt.size = 0)
dev.off()

pdf("EXC_markers_featureplot.pdf", width = 16, height = 14)
FeaturePlot(seurmerg_integrated, features = marks, sort = T, raster = T)
dev.off()

pdf("EXC_integrated_depth.pdf")
ggboxplot(seurmerg_integrated[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf("EXC_integrated_depth_individual.pdf")
ggboxplot(seurmerg_integrated[[]], x = 'orig.ident', y = 'nFeature_RNA', color = 'orig.ident') + NoLegend() +
rotate_x_text(90)
dev.off()

#seurmerg_integrated = subset(seurmerg_integrated, subset = seurat_clusters %in% c(22), invert = T)
saveRDS(seurmerg_integrated, 'EXC_Integrated_Annotated.RDS')

####
## 02 - LABEL TRANSFER USING ALLEN MTG
####

newrna_int = seurmerg_integrated

# Read allen brain excitatory data
allen_mtg_seur = readRDS("exc_forannot_processed.rds")
allen_mtg_seur$cell_type = droplevels(allen_mtg_seur$cluster)

# Check how many genes  will be used for label transfer
DefaultAssay(newrna_int) = 'integrated'
DefaultAssay(allen_mtg_seur) = 'SCT'
VariableFeatures(allen_mtg_seur) = rownames(allen_mtg_seur)
sum(VariableFeatures(allen_mtg_seur) %in% rownames(newrna_int@assays$integrated))

# Find transfer anchors
transfer_anchors <- FindTransferAnchors(reference = allen_mtg_seur, query = newrna_int,
					dims = 1:20,
					features = rownames(allen_mtg_seur),
					reference.assay = "SCT",
					query.assay = "integrated",
					reduction = "cca")

# Make cell type predictions based on anchors
pred <- TransferData(anchorset = transfer_anchors,
			refdata = allen_mtg_seur$cell_type, 
			weight.reduction = newrna_int[["pca"]])

# Add predictions to meta data
newrna_int$pred = pred$predicted.id
newrna_int$pred_score = pred$prediction.score.max

pdf("pred_exc_AllenRef_overlap2k.pdf")
hist(newrna_int$pred_score)
dev.off()

# Plot label transfer
# Make the colors match
newrna_int$pred <- factor(newrna_int$pred, levels = levels(allen_mtg_seur$cell_type))
Idents(newrna_int) = newrna_int$pred
Idents(allen_mtg_seur) = allen_mtg_seur$cell_type
transferPlot(newrna_int, allen_mtg_seur, 'LABEL_TRANSFER_Allen_to_Exc')

# Examine label transfer per cluster. Increase resolution if cluster has two distinct labels.
newrna_int$predicted.id = newrna_int$pred
ltpred(newrna_int, "PREDICTION_AllenRef_excitatory_res0.5")

DefaultAssay(newrna_int) = 'integrated'
newrna_int <- FindNeighbors(newrna_int, dims = 1:30)
newrna_int <- FindClusters(newrna_int, verbose = FALSE, resolution = 0.5)
pdf("rna_exc_int_clusters_res05.pdf")
DimPlot(newrna_int, label = T, group.by = 'seurat_clusters') + NoLegend()
dev.off()
ltpred(newrna_int, "PREDICTION_AllenRef_excitatory_res05")

# Check dendrogram
Idents(newrna_int) = newrna_int$seurat_clusters
newrna_int = BuildClusterTree(newrna_int, assay = 'integrated', dims = 1:30, reorder = T)
dendr = Tool(object = newrna_int, slot = 'BuildClusterTree')

pdf('Exc_tree.pdf', width = 15, height = 7.5)
PlotClusterTree(newrna_int)
dev.off()


####
## 03 - ANNOTATE SUBTYPES
####
mapnames = setNames(c('L3-5_RORB_1', 'L2-3_2', 'L2-3_1', 'L2-3_3', 'L3-5_RORB_2', 'L4-5_RORB_1',
			'L5-6_THEMIS_1', 'L2-3_2', 'L5-6_FEZF2_1', 'L4-6_RORB_3', 'L3-5_RORB_3',
			'L4-6_RORB_3', 'L5-6_FEZF2_2-3', 'L4-6_RORB_2', 'L4-6_FEZF2', 'L5-6_THEMIS_2'),
		      c(0:15))

newrna_int[["newannot"]] <- mapnames[newrna_int[["seurat_clusters"]][,1]]
newrna_int$newannot = factor(newrna_int$newannot)
Idents(newrna_int) = newrna_int$newannot

pdf('Exc_annotated.pdf')
DimPlot(newrna_int, group.by = 'newannot', label = T, repel = T) + NoLegend()
dev.off()


# Save annotated and filtered version
saveRDS(newrna_int, 'EXC_Integrated_Annotated.RDS')


####
## 04 - FIND CELLTYPE MARKERS PER SPECIES
####
newrna_int = readRDS('EXC_Integrated_Annotated.RDS')
humanExc = subset(newrna_int, subset = Species == 'human')
DefaultAssay(humanExc) = 'RNA'
humanExc = NormalizeData(humanExc)
Idents(humanExc) = humanExc$newannot
humanMarkers = FindAllMarkers(humanExc, only.pos = T)

chimpExc = subset(newrna_int, subset = Species == 'chimp')
DefaultAssay(chimpExc) = 'RNA'
chimpExc = NormalizeData(chimpExc)
Idents(chimpExc) = chimpExc$newannot
chimpMarkers = FindAllMarkers(chimpExc, only.pos = T)

macaqueExc = subset(newrna_int, subset = Species == 'macaque')
DefaultAssay(macaqueExc) = 'RNA'
macaqueExc = NormalizeData(macaqueExc)
Idents(macaqueExc) = macaqueExc$newannot
macaqueMarkers = FindAllMarkers(macaqueExc, only.pos = T)

saveRDS(humanMarkers, 'EXC_HumanMarkers.RDS')
saveRDS(chimpMarkers, 'EXC_ChimpMarkers.RDS')
saveRDS(macaqueMarkers, 'EXC_MacaqueMarkers.RDS')


