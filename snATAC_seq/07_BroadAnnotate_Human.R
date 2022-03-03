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
source("utility_functions.R")
set.seed(1234)

####
## PREPARE BOTH ASSAYS FOR LABEL TRANSFER
####

# Read snATAC-seq object and gene activity matrix
integrated = readRDS("human_integrated.RDS")
geneact = readRDS("Gene_Activity_Matrices/human_cicero_geneact.RDS")

# Match gene activity matrix cell names and order
int_names = colnames(integrated@assays$peaks@counts)
geneact = geneact[,match(int_names, colnames(geneact))]

# Add the gene RNA matrix to the Seurat object as a new assay, and normalize it
integrated[['RNA']] <- CreateAssayObject(counts = geneact)

options(future.globals.maxSize = 4000 * 1024^2)
integrated = SCTransform(integrated, ncells = 1000)
DefaultAssay(integrated) = "SCT"

# Load RNA and extract species
hrnaseur = readRDS('Human_Integrated_Filtered2.RDS')

options(future.globals.maxSize = 4000 * 1024^2)
hrnaseur = SCTransform(hrnaseur, ncells = 1000)
Idents(hrnaseur) = hrnaseur$broadannot

####
## LABEL TRANSFER
####

# Check how many genes will be used for integration
sum(VariableFeatures(hrnaseur) %in% rownames(integrated@assays$RNA@counts))

# Find anchors
transfer.anchors <- FindTransferAnchors(reference = hrnaseur, query = integrated,
					dims = 1:20, features = VariableFeatures(hrnaseur),
					reference.assay = "SCT", query.assay = "SCT",
					reduction = "cca")

# Make cell type predictions based on anchors
pred_broad <- TransferData(anchorset = transfer.anchors,
				     refdata = Idents(hrnaseur), 
				     weight.reduction = integrated[["harmony"]])

# Add predictions to meta data
integrated$pred_broad = pred_broad$predicted.id
integrated$pred_broad_score = pred_broad$prediction.score.max
pdf("Human_pred_ATAC_RNA_BROAD.pdf")
hist(integrated$pred_broad_score)
dev.off()

Idents(integrated) = integrated$pred_broad
transferPlot(integrated, hrnaseur, 'LABEL_TRANSFER_human_Rna_Atac_BROAD')
ltpred(integrated, fn = 'HEATMAP_LABEL_TRANSFER_human_Rna_Atac')

####
## ANNOTATION
####

mapnames = setNames(c("Oligos", "Oligos", "Astro", "Oligos", "Exc",
			"Oligos", "OPC", "Inh", "Exc", "Inh",
			"Oligos", "Exc", "Inh", "Exc", "Exc",
			"Oligos", "Exc", "Exc", "Exc", "Microglia",
			"Inh", "Exc", "Oligos", "Exc", "Exc",
			"Exc", "Exc"),
		    c(0:26))

integrated[["broadannot"]] <- mapnames[integrated[["seurat_clusters"]][,1]]

pdf("human_broadannot.pdf")
DimPlot(integrated, group.by = 'broadannot', label = T)
dev.off()

saveRDS(integrated, "human_annotated.RDS")

#####
## FINAL PLOTS
#####

mergedH = readRDS("human_annotated.RDS")

library("colorspace")
q4 <- qualitative_hcl(6, palette = "Set2")

mergedH$broadannot = gsub('Astro', 'Ast', mergedH$broadannot)
mergedH$broadannot = gsub('Oligos', 'Oli', mergedH$broadannot)
mergedH$broadannot = gsub('Microglia', 'Mic', mergedH$broadannot)
mergedH$broadannot = factor(mergedH$broadannot, levels = c('OPC', 'Oli', 'Ast', 'Mic', 'Exc', 'Inh'))

pdf("Human_BroadAnnot_ATAC.pdf")
DimPlot(mergedH, label = T, cols = q4, label.size = 8, raster = T, group.by = 'broadannot') +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		legend.text=element_text(size=20)) + ggtitle('')
dev.off()

# Broaden prediction label
mergedH$predicted.id = gsub('Astrocytes', 'Ast', mergedH$predicted.id)
mergedH$predicted.id = gsub('Oligo', 'Oli', mergedH$predicted.id)
mergedH$predicted.id = gsub('Microglia', 'Mic', mergedH$predicted.id)
mergedH$predicted.id = gsub('Exc_.*', 'Exc', mergedH$predicted.id)
mergedH$predicted.id = gsub('Inh_.*', 'Inh', mergedH$predicted.id)

mergedH$predicted.id = factor(mergedH$predicted.id, levels = c('OPC', 'Oli', 'Ast', 'Mic', 'Exc', 'Inh'))

ltpred(mergedH, fn = 'BROAD_ATAC_Human_PredictionPlot', vars = c('predicted.id', 'broadannot'), heatmap = F, textsize = 30)

pdf("Human_BroadAnnot_ATAC_ReadsInPeaks.pdf")
ggboxplot(mergedH[[]], x = 'broadannot', y = 'rip', fill = 'white', color = 'broadannot', palette = q4, ylim = c(0, 100000),
outlier.shape = NA) +
ylab('Total reads in peaks') + xlab('') +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		legend.text=element_text(size=20)) + NoLegend()
dev.off()

pdf("Human_BroadAnnot_ATAC_Fraction_of_ReadsInPeaks.pdf")
ggboxplot(mergedH[[]], x = 'broadannot', y = 'frip', fill = 'white', color = 'broadannot', palette = q4, ylim = c(0, 1), outlier.shape = NA) +
ylab('Fraction of reads in peaks') + xlab('') +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		legend.text=element_text(size=20)) + NoLegend()
dev.off()


# PREDICTION PLOT -- JACCARD INDEX
meta = mergedH[[]]
meta$broadannot = gsub('Astro', 'Ast', meta$broadannot)
meta$broadannot = gsub('Microglia', 'Mic', meta$broadannot)
meta$broadannot = gsub('Oligos', 'Oli', meta$broadannot)

meta$predicted.id = gsub('Exc.*', 'Exc', meta$predicted.id)
meta$predicted.id = gsub('Inh.*', 'Inh', meta$predicted.id)
meta$predicted.id = gsub('Oligo.*', 'Oli', meta$predicted.id)
meta$predicted.id = gsub('Ast.*', 'Ast', meta$predicted.id)
meta$predicted.id = gsub('Mic.*', 'Mic', meta$predicted.id)

ctypes = names(table(meta$broadannot))
jacmatL = list()
for(i in 1:length(ctypes)){

	annotCells = meta[meta$broadannot %in% ctypes[i],] %>% rownames
	jacL = list()

	for(j in 1:length(ctypes)){
		predCells = meta[meta$predicted.id %in% ctypes[j],] %>% rownames
		jacL[[j]] = sum(annotCells %in% predCells) / length(union(annotCells,predCells))
		#print(jac)
	}
	jacmatL[[i]] = unlist(jacL)
	names(jacmatL[[i]]) = ctypes
}

jacmatDF = do.call(rbind, jacmatL)
rownames(jacmatDF) = ctypes
toplot = melt(jacmatDF)
colnames(toplot) = c('Annotation', 'Prediction', 'Jaccard_Index')

pdf('BROAD_ATAC_Human_PredictionPlot_JACCARD.pdf', width = 10, height =10)
ggscatter(toplot, x = 'Prediction', y = 'Annotation', color = 'Jaccard_Index', size = 10) +
scale_color_gradient2(low = "white", high = "blue", mid = "white", midpoint = 0, limit = c(0,1)) +
xlab('Prediction') + ylab('Annotation') +
theme(text = element_text(size=25)) +
	rotate_x_text(90) + grids(linetype = "dashed", color = 'grey')
dev.off()





