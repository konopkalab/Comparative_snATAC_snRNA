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
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(harmony)
set.seed(1234)

####
## GENERATE SEURAT OBJECTS
####

# Human #
samps = 1:4
humseurs = list()
for(i in 1:length(samps)){

	# Read the matrix
	counts_hum <- readMM(paste0("human", i, "/pc_mat.mtx"))
	counts_hum = as(counts_hum, "dgCMatrix")
	counts_hum = t(counts_hum)
	barcs_hum = read.table(paste0("human", i, "/barcodes.tsv"))
	peaks_hum = read.table(paste0("human", i, "/peaks_notab.txt"))
	colnames(counts_hum) = barcs_hum$V1
	rownames(counts_hum) = peaks_hum$V1

	# Read metadata
	meta_hum = read.table(paste0("human", i, "/human", i, "_metadata.tsv"), header = T)
	rownames(meta_hum) = meta_hum$barcode
	meta_hum$Species = "Human"

	# Add barcode multiplets
	barcmult = read.csv(paste0("human", i, "/excluded_barcodes.csv"), header = T)
	barcmult[,1] = gsub("-1", paste0("_h",i), barcmult[,1])
	barcmult[,2] = gsub("-1", paste0("_h",i), barcmult[,2])
	meta_hum$is_multiplet = ifelse(colnames(counts_hum) %in% barcmult$Excluded.Barcode, 1, 0)

	# Binarize the counts
	tmpall = counts_hum
	starts = seq(1,nrow(tmpall),50000)
	ends = c(starts[1:(length(starts)-1)] + 50000 - 1, nrow(tmpall))
	tmps = list()

	for(j in 1:length(starts)){

		tmps[[j]] = tmpall[starts[j]:ends[j],]
		tmps[[j]][tmps[[j]] > 0] = 1
		gc()

	}

	tmpall_bin = do.call(rbind, tmps)
	counts_hum = tmpall_bin


	# Create and filter seurat objects
	hum_seur = CreateSeuratObject(
	  counts = counts_hum,
	  assay = 'peaks',
	  project = 'ATAC',
	  min.cells = 1,
	  meta.data = meta_hum
	)


	hum_seur[["Sample"]] = stri_sub(rownames(hum_seur[[]]), -2,-1)

	humseurs[[i]] = hum_seur
}

allhum_seur = Reduce(merge, humseurs)

####
## QUALITY CONTROL AND FILTERING
####

meta = allhum_seur[[]]
meta$Cell = ifelse(meta$rip > 3000 &
			meta$frip > 0.15 &
			meta$rip < 100000 &
			meta$is_multiplet == 0, "Cell", "Non-cell")


meta$Sample = factor(meta$Sample, levels = c(paste0('h', 1:4)))
pdf("Frip_plot.pdf", width = 8)
ggscatter(data = meta %>% arrange(desc(is_multiplet)),
	 x = 'rip',
	 y = 'frip',
	 xlab = "Total reads in peaks",
	 ylab = "Fraction of reads in peaks",
	 alpha = 0.5,
	 color = "Cell",
	 palette = c("blue", "red"),
	 size = 2) + theme_classic() + 
	theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		strip.text = element_text(size = 20),
		legend.text=element_text(size=20),
		legend.title = element_blank(), legend.position = 'right') +
	facet_wrap(~Sample) + rotate_x_text(angle=90)
dev.off()

allhum_seur = subset(allhum_seur, subset = frip > 0.15 & rip > 3000 & rip < 100000 & is_multiplet == 0)


####
## BATCH CORRECTION AND CLUSTERING PER SPECIES
####

sps = c("Human")

# Extract one species
tmpseur = subset(merged_seur, subset = Species == sps)

# Dimensionality reduction
tmpseur <- RunTFIDF(tmpseur, assay = 'peaks')
tmpseur <- FindTopFeatures(tmpseur, min.cutoff = 'q50', assay = 'peaks')
tmpseur <- RunSVD(
  object = tmpseur,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)


unintegrated <- RunUMAP(tmpseur, reduction = 'lsi', dims = 1:30)

pdf(paste0(sps, "_preintegration_samples.pdf"))
DimPlot(unintegrated, group.by = 'Sample', pt.size = 0.1, label.size = 7) +
theme(axis.text.x = element_text(size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20))
dev.off()

# Batch correction
all_integrated <- RunHarmony(
  object = unintegrated,
  group.by.vars = 'Sample',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

# UMAP after batch correction
all_integrated <- RunUMAP(all_integrated, dims = 1:30, reduction = 'harmony')
pdf(paste0(sps, "_harmony_samples.pdf"))
DimPlot(all_integrated, group.by = 'Sample', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
dev.off()

# Clustering
all_integrated <- FindNeighbors(object = all_integrated, reduction = 'harmony', dims = 2:30)
all_integrated <- FindClusters(object = all_integrated, verbose = FALSE, resolution = 1)

Idents(all_integrated) = 'seurat_clusters'
pdf(paste0(sps, "_harmony_clusters.pdf"))
DimPlot(all_integrated, label = T, reduction = 'umap') + NoLegend()
dev.off()

saveRDS(all_integrated, "human_integrated.RDS")










