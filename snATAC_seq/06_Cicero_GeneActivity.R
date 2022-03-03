rm(list = ls())
library(cicero)
library(Matrix)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(Signac)
source("custom_functions.R")

####
## FIND CONNECTIONS BETWEEN CREs
####

# Read the merged matrix from each species
humanseur = readRDS("human_integrated.RDS")

## HUMAN ##

meta = humanseur[[]]
mat = humanseur@assays$peaks@counts

# Remove mitochondria and sex chromosomes peaks
mat = mat[!grepl("chrM|chrY|chrX", rownames(mat)),]

# Remove zero or low accessible peaks
humpeaks = which(rowSums(mat) / ncol(mat) > 0.01) %>% names
length(humpeaks)

# Subset the peaks
hummat = mat[humpeaks,]

# Umap coordinates of given species
umap_coords = humanseur@reductions$umap@cell.embeddings

# Prepare cell barcodes in cicero format
cells = colnames(hummat) %>% as.data.frame()
rownames(cells) = cells[,1]
names(cells) = "cells"

# Prepare peaks in cicero format
peakinfo = rownames(hummat) %>% gsub(":|-", " ", .) %>%
		strsplit(., " ") %>% do.call(rbind, .) %>% as.data.frame()
names(peakinfo) = c("chr", "bp1", "bp2")
peakinfo$site_name = paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) = peakinfo$site_name

# Create CellDataSet (cds)
rownames(hummat) = rownames(peakinfo)
colnames(hummat) = rownames(cells)
fd = methods::new("AnnotatedDataFrame", data = peakinfo)
pd = methods::new("AnnotatedDataFrame", data = cells)

input_cds = newCellDataSet(hummat,
		phenoData = pd,
		featureData = fd,
		expressionFamily=VGAM::binomialff(),
		lowerDetectionLimit=0)

input_cds@expressionFamily@vfamily = "binomialff"
input_cds = monocle::detectGenes(input_cds)

# Make sure there are no peaks with zero accessibility
input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

# Create cicero cds which merges similar cells together to compute cis co-accessibility
human_cds = make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Calculate connections
hg38_sizes = read.table("Chrom_sizes/hg38.chrom.sizes")
human_conns = run_cicero(human_cds, hg38_sizes)
human_conns = human_conns[!is.na(human_conns$coaccess),]

saveRDS(human_conns, "Gene_Activity_Matrices/human_conns.RDS")
saveRDS(human_cds, "Gene_Activity_Matrices/human_cicero_cds.RDS")
saveRDS(input_cds, "Gene_Activity_Matrices/human_input_cds.RDS")


####
## COMPUTE GENE ACTIVITY MATRIX
####

# Reads cds and conns
human_conns = readRDS("Gene_Activity_Matrices/human_conns.RDS")
human_inputcds = readRDS("Gene_Activity_Matrices/human_input_cds.RDS")

# Keep only the transcripts of protein coding genes
hum_gtf = read.table("Hsa_GRCh38.gtf", sep = "\t", skip = 5, stringsAsFactors = F)
hum_genes = hum_gtf[(hum_gtf$V3 == "transcript") & grepl("protein_coding", hum_gtf$V9),]

# Reshape to expand metadata
tmp = hum_genes$V9 %>% gsub(";", "", .) %>% strsplit(., " ") %>% do.call(rbind, .)
hum_genes$V9 = NULL
hum_genes = cbind(hum_genes, tmp)

pos <- subset(hum_genes, V7 == "+")
pos <- pos[order(pos$V4),]
pos <- pos[!duplicated(pos[,14]),] # remove all but the first exons per transcript
pos$V5 <- pos$V4 + 1 # make a 1 base pair marker of the TSS

neg <- subset(hum_genes, V7 == "-")
neg <- neg[order(neg$V4, decreasing = TRUE),]
neg <- neg[!duplicated(neg[,14]),] # remove all but the first exons per transcript
neg$V4 <- neg$V5 - 1

gene_annotation_sub <- rbind(pos, neg)
hum_genes = gene_annotation_sub[, c(1,4,5,18)]


# Get orthologous genes #
orthologs = read.table("hum_chimp_mac_homologs_biomart_April2021.txt",
			stringsAsFactors = F, sep = "\t", header = T)
orthologs[,2] = NULL
colnames(orthologs) = c("Hum_ID", "Chimp_ID", "Mac_ID", "Gene_name")

# Find genes orthologous across all three species
orthologs = orthologs[orthologs[,1] != "" & orthologs[,2] != "" & orthologs[,3] != "",]

# Keep only orthologous genes
hum_genes = hum_genes[hum_genes[,4] %in% orthologs$Gene_name, ]
colnames(hum_genes) = c("chromosome", "start", "end", "gene")

saveRDS(hum_genes, 'Gene_Activity_Matrices/Genes_For_Cicero.RDS')

# Add to cds
human_inputcds <- annotate_cds_by_site(human_inputcds, hum_genes)

# Generate unnormalized gene activity matrix
human_unnorm_ga <- build_gene_activity_matrix(human_inputcds, human_conns)

# Remove any rows/columns with all zeroes
human_unnorm_ga <- human_unnorm_ga[!Matrix::rowSums(human_unnorm_ga) == 0, !Matrix::colSums(human_unnorm_ga) == 0]

saveRDS(human_unnorm_ga, "Gene_Activity_Matrices/human_cicero_geneact.RDS")

