rm(list = ls())
library(dplyr)
library(plyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
source("utility_functions.R")

####
## READ MATRICES
####

fls_filt = list.files(path = 'CELLBENDER/', pattern = 'out_filtered.h5', recursive = T, full.names = T)

fls_filt_human = fls_filt[grepl('HUMAN', fls_filt)]

# Load ensembl data for orthologs
orthologs = read.table("hum_chimp_mac_homologs_biomart_April2021.txt",
			stringsAsFactors = F, sep = "\t", header = T)
colnames(orthologs) = c("Hum_ID", "Chimp_ID", "Mac_ID", "Gene_name")

# Find genes orthologous across all three species
orthologs = orthologs[orthologs[,1] != "" & orthologs[,2] != "" & orthologs[,3] != "",]

# Keep only protein coding among orthologous
gnsids = genes(EnsDb.Hsapiens.v86) %>% as.data.frame()
gnspr = gnsids[gnsids$seqnames %in% c(1:22, 'X', 'Y', 'MT') & gnsids$gene_biotype == "protein_coding",]
orthologs_pr = orthologs[orthologs$Hum_ID %in% gnspr$gene_id,]

seurs = list()
for(i in 1:length(fls_filt_human)){

	mat = Read10X_h5(fls_filt_human[i])
	
	# Mitochrondrial percentage before filtering
	mitgns = gnsids[grepl("^MT-", gnsids$symbol), 'symbol']
	mtsum = mat[rownames(mat) %in% mitgns,]	%>% colSums
	mtperc = mtsum / colSums(mat) * 100
	mtperc[is.na(mtperc)] = 0

	# Filter the matrix
	orthologs_pr_nomt = orthologs_pr[!grepl('^MT-', orthologs_pr$Gene_name),]
	mat = mat[rownames(mat) %in% orthologs_pr_nomt$Gene_name,]
	rownames(mat) = make.unique(rownames(mat))

	# Create seurat object per sample
	seurs[[i]] = CreateSeuratObject(counts = mat,
					project = paste0('h', i),
					min.cells = 10,
					min.features = 0)
	seurs[[i]]$percent_mt = mtperc

	seurs[[i]]$Species = 'human'
	seurs[[i]] = RenameCells(seurs[[i]], new.names = paste(colnames(seurs[[i]]), paste0('h', i), sep = '_'))
	print(i)
}

seurmerg = Reduce(merge, seurs)

####
## QUALITY CONTROL
####

# QC Plots
meta = seurmerg[[]]
meta2 = melt(meta)

pdf("MitPerc_Human.pdf", width = 10, height = 8)
ggboxplot(meta, x = 'orig.ident', y = "percent_mt", fill = 'orig.ident') +
	xlab("") +
	ylab("% mitochondria") +
	scale_y_continuous(breaks = seq(0,100,5)) +
	theme_classic() +
	theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
	rotate_x_text(90)
dev.off()

# Subset by mitochondria percentage & feature count
seurmergFilt = subset(seurmerg, subset = percent_mt < 5 & nCount_RNA > 200)

# Sample size
meta = seurmergFilt[[]]
ssize = meta[, c("orig.ident", "Species")]
ssize$cnt = 1
ssize = aggregate(. ~ orig.ident + Species, ssize, FUN=sum)

pdf("RNA_sample_sizes.pdf")
ggbarplot(data = ssize,
	 x = 'orig.ident',
	 y = 'cnt',
	 fill = 'grey',
	 ylab = "Total Cell Number",
	 xlab = "",
	 color = "grey") +
	rotate_x_text(90) +
	theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		strip.text = element_text(size=20))
dev.off()


# Plot before integration
options(future.globals.maxSize = 4000 * 1024^2)
seurmergFilt = SCTransform(seurmergFilt, ncells = 1000, vars.to.regress = 'percent_mt')
seurmergFilt <- RunPCA(seurmergFilt, verbose = FALSE)
seurmergFilt <- RunUMAP(seurmergFilt, dims = 1:30)

pdf("human_before_harmony_umap.pdf")
DimPlot(seurmergFilt, group.by = 'orig.ident')
dev.off()

# Batch correction
library(harmony)
seurmergFilt <- RunHarmony(object = seurmergFilt, group.by.vars = 'orig.ident', assay.use="SCT")

# After Harmony
seurmergFilt <- RunUMAP(seurmergFilt, dims = 1:30, reduction = 'harmony')
seurmergFilt <- FindNeighbors(seurmergFilt, dims = 1:30, verbose = FALSE, reduction = 'harmony')
seurmergFilt <- FindClusters(seurmergFilt, verbose = FALSE, resolution = 1)

pdf("human_after_harmony_umap.pdf")
DimPlot(seurmergFilt, group.by = 'orig.ident')
dev.off()

pdf("human_after_harmony_clusters.pdf")
DimPlot(seurmergFilt, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

library(tidyverse)
stackedbarplot(seurmergFilt[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', 'Sample_Stacked_Human')

# Marker Genes
DefaultAssay(seurmergFilt) = 'RNA'
seurmergFilt = NormalizeData(seurmergFilt)
marks = c('APBB1IP', 'SOX6', 'MOG', 'BCAS1', 'OPALIN', 'AQP4', 'SLC17A7', 'GAD1', 'RBFOX3', 'NRGN', 'CALM1', 'FLT1', 'PCDH15', 'LUZP2', 'GFAP', 'SLC1A2', 'SLC1A3', 'GJA1', 'EBF1', 'COBLL1')

pdf("human_markers_vlnplot.pdf", width = 22, height = 10)
VlnPlot(seurmergFilt, features = marks, group.by = 'seurat_clusters', pt.size = 0)
dev.off()

pdf("human_depth.pdf")
ggboxplot(seurmergFilt[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf("human_mito_perc.pdf")
ggboxplot(seurmergFilt[[]], x = 'seurat_clusters', y = 'percent_mt', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()


allmarks = FindAllMarkers(seurmergFilt, only.pos = T, min.pct = 0.5, logfc.threshold = 1)


saveRDS(seurmergFilt, 'CELLBENDER/Human_Integrated_Filtered.RDS')

####
## CLUSTER FILTER 1
####

seurmergFilt = readRDS('CELLBENDER/DATA/Human_Integrated_Filtered.RDS')

# Remove possible doublets / mixed clusters as well as endothelial cells
filtcl = c(14, 17, 20, 25, 29, 30, 31, 33:36)
seurmergFilt2 = subset(seurmergFilt, subset = seurat_clusters %in% filtcl, invert = T)

# Rerun the clustering
options(future.globals.maxSize = 4000 * 1024^2)
seurmergFilt2 = SCTransform(seurmergFilt2, ncells = 1000, vars.to.regress = 'percent_mt')
seurmergFilt2 <- RunPCA(seurmergFilt2, verbose = FALSE)
seurmergFilt2 <- RunUMAP(seurmergFilt2, dims = 1:30)

# Batch correction
library(harmony)
seurmergFilt2<- RunHarmony(object = seurmergFilt2, group.by.vars = 'orig.ident', assay.use="SCT")

# After Harmony
seurmergFilt2 <- RunUMAP(seurmergFilt2, dims = 1:30, reduction = 'harmony')
seurmergFilt2 <- FindNeighbors(seurmergFilt2, dims = 1:30, verbose = FALSE, reduction = 'harmony')
seurmergFilt2 <- FindClusters(seurmergFilt2, verbose = FALSE, resolution = 1)

pdf("human_after_harmony_umap.pdf")
DimPlot(seurmergFilt2, group.by = 'orig.ident')
dev.off()

pdf("human_after_harmony_clusters.pdf")
DimPlot(seurmergFilt2, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

library(tidyverse)
stackedbarplot(seurmergFilt2[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', 'Sample_Stacked_Human')

# Marker Genes
DefaultAssay(seurmergFilt2) = 'RNA'
seurmergFilt2 = NormalizeData(seurmergFilt2)

pdf("human_markers_vlnplot.pdf", width = 22, height = 10)
VlnPlot(seurmergFilt2, features = marks, group.by = 'seurat_clusters', pt.size = 0)
dev.off()


pdf("human_depth.pdf")
ggboxplot(seurmergFilt2[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()


pdf("human_mito_perc.pdf")
ggboxplot(seurmergFilt2[[]], x = 'seurat_clusters', y = 'percent_mt', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()


####
## CLUSTER FILTER 2
####

# Remove possible doublets / mixed clusters as well as endothelial cells
filtcl = c(20, 28)
seurmergFilt2 = subset(seurmergFilt2, subset = seurat_clusters %in% filtcl, invert = T)

# Rerun the clustering
options(future.globals.maxSize = 4000 * 1024^2)
seurmergFilt2 = SCTransform(seurmergFilt2, ncells = 1000)
seurmergFilt2 <- RunPCA(seurmergFilt2, verbose = FALSE)
seurmergFilt2 <- RunUMAP(seurmergFilt2, dims = 1:30)

library(harmony)
seurmergFilt2<- RunHarmony(object = seurmergFilt2, group.by.vars = 'orig.ident', assay.use="SCT")

# After Harmony
seurmergFilt2 <- RunUMAP(seurmergFilt2, dims = 1:30, reduction = 'harmony')
seurmergFilt2 <- FindNeighbors(seurmergFilt2, dims = 1:30, verbose = FALSE, reduction = 'harmony')
seurmergFilt2 <- FindClusters(seurmergFilt2, verbose = FALSE, resolution = 1)

pdf("human_after_harmony_umap.pdf")
DimPlot(seurmergFilt2, group.by = 'orig.ident')
dev.off()

pdf("human_after_harmony_clusters.pdf")
DimPlot(seurmergFilt2, group.by = 'seurat_clusters', label = T) + NoLegend()
dev.off()

library(tidyverse)
stackedbarplot(seurmergFilt2[[]], groupx = 'seurat_clusters', groupfill = 'orig.ident', 'Sample_Stacked_Human')

# Marker Genes
DefaultAssay(seurmergFilt2) = 'RNA'
seurmergFilt2 = NormalizeData(seurmergFilt2)

pdf("human_markers_vlnplot.pdf", width = 22, height = 10)
VlnPlot(seurmergFilt2, features = marks, group.by = 'seurat_clusters', pt.size = 0)
dev.off()


pdf("human_depth.pdf")
ggboxplot(seurmergFilt2[[]], x = 'seurat_clusters', y = 'nFeature_RNA', color = 'seurat_clusters') + NoLegend() +
rotate_x_text(90)
dev.off()

pdf('human_cluster_facet.pdf')
DimPlot(seurmergFilt2, group.by = 'seurat_clusters', label = T) + NoLegend() + facet_wrap(~seurat_clusters)
dev.off()

saveRDS(seurmergFilt2, 'Human_Integrated_Filtered2.RDS')


####
## BROADLY ANNOTATE
####

mapnames = setNames(c('Exc', 'Exc', 'Oli', 'Inh', 'Exc', 'Inh',
			'Inh', 'Exc', 'Exc', 'OPC', 'Ast',
			'Exc', 'Mic', 'Exc', 'Oli', 'Exc',
			'Inh', 'Inh', 'Inh', 'Exc', 'Exc',
			'Exc', 'Exc', 'Exc', 'Inh', 'Inh',
			'Inh'),
		      c(0:26))

seurmergFilt2[["broadannot"]] <- mapnames[seurmergFilt2[["seurat_clusters"]][,1]]
seurmergFilt2$broadannot = factor(seurmergFilt2$broadannot)

pdf('Human_integrated_annotated_broad.pdf')
DimPlot(seurmergFilt2, group.by = 'broadannot', label = T) + NoLegend()
dev.off()

seurmergFilt2$broadannot = factor(seurmergFilt2$broadannot, levels = c("OPC", "Oli", "Ast", "Mic", "Exc", "Inh"))
stackedbarplot(seurmergFilt2[[]], groupx = 'broadannot', groupfill = 'orig.ident', 'BroadAnnot_Sample_Stacked_Human')

saveRDS(seurmergFilt2, 'Human_Integrated_Filtered2.RDS')


