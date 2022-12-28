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
library(lmtest)
library("pheatmap")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Ptroglodytes.UCSC.panTro5")
library("BSgenome.Mmulatta.UCSC.rheMac10")
source("utility_functions.R")


####
## CREATE MOTIF OBJECTS IN SIGNAC
####

humanseur = readRDS('human_annotated.RDS')
chimpseur = readRDS('chimp_annotated.RDS')
macaqueseur = readRDS('macaque_annotated.RDS')

DefaultAssay(humanseur) = 'peaks'
DefaultAssay(chimpseur) = 'peaks'
DefaultAssay(macaqueseur) = 'peaks'

# HUMAN
counts = GetAssayData(humanseur, slot = "counts")
chrobj = CreateChromatinAssay(counts = counts, genome = "hg38", sep = c(":", "-"))

seurchrobj = CreateSeuratObject(counts = chrobj, assay = "peaks", meta.data = humanseur[[]])
pfm = getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))

peaks = rownames(seurchrobj) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame
peaksgr = makeGRangesFromDataFrame(df = peaks, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1')

seurchrobj = AddMotifs(object = seurchrobj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
saveRDS(seurchrobj, 'motif_obj_human_all.RDS')

# CHIMP - CONVERT CREs TO CHIMP COORDINATES
counts = GetAssayData(chimpseur, slot = "counts")
chimp_coords = read.table('merged_lifted_chimp.bed')
chimp_coords$V5 = sub('-', ':', chimp_coords$V5)
chimp_coords$V1 = sub('-', ':', chimp_coords$V1)
chimprows = chimp_coords[match(rownames(counts), chimp_coords$V1), 'V5']

chimp_counts = counts
rownames(chimp_counts) = chimprows
chrobj = CreateChromatinAssay(counts = chimp_counts, sep = c(":", "-"))

seurchrobj = CreateSeuratObject(counts = chrobj, assay = "peaks", meta.data = chimpseur[[]])
pfm = getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))

peaks = rownames(seurchrobj) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame
peaksgr = makeGRangesFromDataFrame(df = peaks, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1')

seurchrobj = AddMotifs(object = seurchrobj, genome = BSgenome.Ptroglodytes.UCSC.panTro5, pfm = pfm)
saveRDS(seurchrobj, 'motif_obj_chimp_all.RDS')


# MACAQUE - CONVERT CREs TO MACAQUE COORDINATES
counts = GetAssayData(macaqueseur, slot = "counts")
macaque_coords = read.table('merged_lifted_macaque.bed')
macaque_coords$V5 = sub('-', ':', macaque_coords$V5)
macaque_coords$V1 = sub('-', ':', macaque_coords$V1)
macaquerows = macaque_coords[match(rownames(counts), macaque_coords$V1), 'V5']

macaque_counts = counts
rownames(macaque_counts) = macaquerows
chrobj = CreateChromatinAssay(counts = macaque_counts, sep = c(":", "-"))

seurchrobj = CreateSeuratObject(counts = chrobj, assay = "peaks", meta.data = macaqueseur[[]])
pfm = getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))

peaks = rownames(seurchrobj) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame
peaksgr = makeGRangesFromDataFrame(df = peaks, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1')

seurchrobj = AddMotifs(object = seurchrobj, genome = BSgenome.Mmulatta.UCSC.rheMac10, pfm = pfm)
saveRDS(seurchrobj, 'motif_obj_macaque_all.RDS')


####
## PREPARE FOR MOTIF ENRICHMENT
####

# Load motif objects
humanobj = readRDS('motif_obj_human_all.RDS')
Idents(humanobj) = humanobj$broadannot

# Extract motif matrix
humanMotMat = humanobj@assays$peaks@motifs@data
colnames(humanMotMat) = humanobj@assays$peaks@motifs@motif.names

# Load background peaks
finalResDF = readRDS('PSEUDOBULK_DARs_ALL.RDS')
dars = finalResDF
dars$peak = dars$Gene

####
## RUN ENRICHMENT
####

# Calculate length
allDARs = dars$peak %>% strsplit(., ':|-') %>% do.call(rbind, .) %>% as.data.frame
dars$Length = as.numeric(allDARs[,3]) - as.numeric(allDARs[,2])

# Match peak names
dars$peak2 = gsub(':', '-', dars$peak)
ctypes = unique(dars$CellType)

## HUMAN-SPECIFIC-UP ##
resdf_hupL = list()
resdfSign_hup = list()
for(i in 1:length(ctypes)){

	darsSub = dars[dars$CellType == ctypes[i],]

	# Remove non-divergent CREs
	darsSubHS = darsSub
	darsSubHS = darsSubHS[darsSubHS$Evolution != 'NS',]

	# Divide into two: specific and non-specific
	darsSubHS$Regulation = gsub('Chimp_UP|Chimp_DOWN|Macaque_UP|Macaque_DOWN|Human_DOWN', 'NS', darsSubHS$Regulation)
	darsSubHS$Regulation = factor(darsSubHS$Regulation, levels = c('NS', 'Human_UP'))

	# Fit two models and find whether motif information results in a better fit
	motMatSub = humanMotMat[darsSubHS$peak2,]
	pvalsL = list()
	coefsL = list()
	fit0L = list()
	for(j in 1:ncol(motMatSub)){

		fit0 = glm(Regulation ~ Length + motMatSub[,j], data = darsSubHS, family = "binomial")
		fit1 = glm(Regulation ~ Length, data = darsSubHS, family = "binomial")
		pval = lrtest(fit0, fit1)[2,5]

		fit0L[[j]] = fit0
		pvalsL[[j]] = pval
		coefsL[[j]] = coef(summary(fit0))[3,1]

		if(j%%50 == 0){print(j)}
	}

	# Combine p values
	hs_pvals = unlist(pvalsL)
	hs_coefs = unlist(coefsL)
	
	# Calculate odds ratio from the matrix
	hsPeaks = darsSubHS[darsSubHS$Regulation == 'Human_UP', 'peak2']
	oddsRatio = colMeans(motMatSub[hsPeaks,]) / colMeans(motMatSub[darsSubHS$peak2,])

	# results data frame
	resdf = data.frame(motif = colnames(motMatSub), CellType = ctypes[i],
				pvalue = hs_pvals, oddsRatio  = oddsRatio,
				coefficient = hs_coefs,
				percentObserved = colMeans(motMatSub[hsPeaks,])*100,
				percentBackground = colMeans(motMatSub[darsSubHS$peak2,])*100)

	resdf$FDR = p.adjust(resdf$pval, method = 'BH')
	resdf$SignificanceType = 'NS'
	resdf[resdf$FDR < 0.05 & resdf$oddsRatio > 1, 'SignificanceType'] = 'Enrichment'
	resdf[resdf$FDR < 0.05 & resdf$oddsRatio < 1, 'SignificanceType'] = 'Depleted'

	resdf_hupL[[i]] = resdf

	# Save results per cell type
	saveRDS(resdf_hupL, paste0('PSEUDOBULK_HUMAN_UP_', ctypes[i], '.RDS'))

	print(i)
}

resDF_hup = do.call(rbind, resdf_hupL)
saveRDS(resDF_hup, 'PSEUDOBULK_HUMAN_UP_MOTIF_ALLDIVERGENT.RDS')
resDF_hup = resDF_hup[grepl('^L[0-9]', resDF_hup$CellType),]


####
## Plot all significant enrichments
####

signMots = resDF_hup[resDF_hup$SignificanceType == 'Enrichment', 'motif'] %>% unique

toplot = resDF_hup[resDF_hup$motif %in% signMots,]
toplot$FDRSc = formatC(toplot$FDR, format = "e", digits = 0)
toplot$log10FDR = -log10(toplot$FDR)

# Retain only excitatory or inhibitory neurons
toplot$CellType= factor(toplot$CellType, levels = c('L2-3_1', 'L2-3_2', 'L2-3_3', 'L3-5_RORB_2', 'L3-5_RORB_1', 'L3-5_RORB_3', 'L4-5_RORB_1', 'L4-6_RORB_3', 'L4-6_RORB_2', 'L5-6_THEMIS_1', 'L5-6_THEMIS_2', 'L5-6_FEZF2_1', 'L5-6_FEZF2_2-3', 'L4-6_FEZF2'))

toplotP = toplot[, c('motif', 'CellType', 'log10FDR', 'coefficient')]
tmp = dcast(data = toplotP, formula = motif~CellType, value.var = "coefficient")
rownames(tmp) = tmp[,1]
tmp[,1] = NULL

newrown = lapply(rownames(tmp), function(x) {bquote(bold(.(x)))} )
newcoln = lapply(colnames(tmp), function(x) {bquote(bold(.(x)))} )

pdf('Motif_Coefficients_Pheatmap_HumanUP.pdf', width = 8, height = 16)
pheatmap(tmp, cluster_cols = F, fontsize_row = 10, fontsize_col = 14,
labels_row = as.expression(newrown), labels_col = as.expression(newcoln))
dev.off()


####
## ACCESSIBILITY DOTPLOT FOR FOX AND FOS/JUN TFs
####

# Read the ATAC-seq data
excATAC = readRDS('human_annotated.RDS')
DefaultAssay(excATAC) = 'RNA'
excATAC = NormalizeData(excATAC)

# All members of the family (only with enriched motifs)
fosjuns = c('FOSL2::JUND', 'FOSL2::JUNB', 'FOSB::JUNB', 'FOS::JUN', 'NFE2',
		'JUN(var.2)', 'FOSL1::JUND', 'FOSL2::JUN', 'FOS::JUND', 'FOSL1::JUN',
		'JUND', 'FOSL1::JUNB', 'FOS', 'JUNB', 'FOSL1', 'FOSL2')

foxs = c('FOXK2', 'FOXK1', 'FOXP2', 'FOXI1', 'FOXD2', 'FOXL1', 'FOXO4', 'FOXO6',
	'FOXP3', 'FOXP1', 'FOXO3', 'FOXD1', 'FOXG1')


# FOX
toplotFOX = DotPlot(excATAC, features = intersect(rownames(excATAC), foxs), group.by = 'newannot')
toplotFOX = toplotFOX$data
rownames(toplotFOX) = paste0(toplotFOX$id, '_', toplotFOX$features.plot)

pdf('EXC_HUMAN_UP_FOX_ACCESSIBILITY.pdf', width = 8, height = 6)
ggscatter(toplotFOX, y = 'id', x = 'features.plot', color = 'avg.exp', size = 'pct.exp') +
  labs(x="",y="") + scale_size_continuous(range = c(0.2,10)) +
  scale_color_gradient2(midpoint = 1.5, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(90)
dev.off()


# FOS
toplotFOS = DotPlot(excATAC, features = intersect(rownames(excATAC), fosjuns), group.by = 'newannot')
toplotFOS = toplotFOS$data
rownames(toplotFOS) = paste0(toplotFOS$id, '_', toplotFOS$features.plot)

pdf('EXC_HUMAN_UP_FOS_ACCESSIBILITY.pdf', width = 8, height = 6)
ggscatter(toplotFOS, y = 'id', x = 'features.plot', color = 'avg.exp', size = 'pct.exp') +
  labs(x="",y="") + scale_size_continuous(range = c(0.2,10)) +
  scale_color_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(90)
dev.off()


####
## SCATTER PLOTS FOR FOX AND FOS/JUN TFs
####

# FOS JUN TFs (enriched and accessible). #
fosjunsAcc = c('FOSL2::JUND', 'FOSL2::JUNB', 'FOS::JUND',
		'JUND', 'FOS', 'JUNB', 'FOSL2')


toplot = resDF_hup[resDF_hup$motif %in% fosjunsAcc,]
toplot$log10FDR = -log10(toplot$FDR)

# Combine p values for similar TFs
library(survcomp)
ctypes = unique(toplot$CellType)
combP = list()
for(i in 1:length(ctypes)){
	pvals = toplot[toplot$CellType == ctypes[i], 'pvalue']
	sizes = toplot[toplot$CellType == ctypes[i], 'percentObserved']
	combP[[i]] = combine.test(pvals, weight = sizes, method = 'fisher')
}

fosFDRs = p.adjust(unlist(combP), method = 'fdr')
names(fosFDRs) = ctypes
fosFDRs = as.data.frame(fosFDRs)
fosFDRs$CellType = rownames(fosFDRs)

toplot = toplot %>% group_by(CellType) %>%
	summarize(log10FDR = median(log10FDR), coefficient = median(coefficient)) %>%
	as.data.frame

# Add the combined p-value
toplot = merge(toplot, fosFDRs, by = 'CellType')
toplot$log10FDR_Fisher = -log10(toplot$fosFDRs)

# Add the number of hs detected per cell type
hsN = table(dars[dars$Regulation == 'Human_UP', 'CellType']) %>% as.data.frame
colnames(hsN) = c('CellType', 'Size')
toplot = merge(toplot, hsN, by = 'CellType')
toplot$CellType2 = ifelse(toplot$log10FDR_Fisher > 1.3 & toplot$coefficient > 0, toplot$CellType, '')

# Scatter plot
pdf('FOS_JUN_MOTIF_ENRICHMENT_SCATTERPLOT.pdf', width = 9)
ggscatter(toplot, y = 'coefficient', x = 'log10FDR_Fisher', color = 'coefficient', label = 'CellType2', size = 'Size', font.label = c(20, 'bold', 'black'), repel = T, ylim = c(-0.3, 0.3)) +
  labs(x="-log10(FDR)",y="Coefficient") + scale_size_continuous(range = c(2,10)) +
  scale_color_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold'),
	axis.text.x = element_text(size=20, face = 'bold'),
	axis.text.y = element_text(size=20, face = 'bold')) +
  geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
  rotate_x_text(90)
dev.off()

# Heatmap (no p-value or coefficient combining)
toplot = resDF_hup[resDF_hup$motif %in% fosjunsAcc,]
toplot$log10FDR = -log10(toplot$FDR)
toplot$coefficient2 = round(toplot$coefficient, digits = 2)

toplot$CellType= factor(toplot$CellType, levels = c('L2-3_1', 'L2-3_2', 'L2-3_3', 'L3-5_RORB_2', 'L3-5_RORB_1', 'L3-5_RORB_3', 'L4-5_RORB_1', 'L4-6_RORB_3', 'L4-6_RORB_2', 'L5-6_THEMIS_1', 'L5-6_THEMIS_2', 'L5-6_FEZF2_1', 'L5-6_FEZF2_2-3', 'L4-6_FEZF2'))

pdf('FOS_JUN_MOTIF_ENRICHMENT_HEATMAP.pdf', width = 14)
ggplot(toplot,aes(x = CellType, y = motif, fill = log10FDR)) +
  geom_tile(colour="white", size=0.2) +
  geom_text(aes(label=coefficient2), fontface = 'bold', size = 7) +
  labs(x="",y="") +
  scale_fill_gradient2(midpoint = 1.3, low = 'blue', high = 'red') +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(90)
dev.off()


# FOX TFs (enriched and accessible) #
foxsAcc = c('FOXP1', 'FOXP2', 'FOXG1')

toplot = resDF_hup[resDF_hup$motif %in% foxsAcc,]
toplot$log10FDR = -log10(toplot$FDR)

# Combine p values for similar TFs
library(survcomp)
ctypes = unique(toplot$CellType)
combP = list()
for(i in 1:length(ctypes)){
	pvals = toplot[toplot$CellType == ctypes[i], 'pvalue']
	sizes = toplot[toplot$CellType == ctypes[i], 'percentObserved']
	combP[[i]] = combine.test(pvals, weight = sizes, method = 'fisher')
}

foxFDRs = p.adjust(unlist(combP), method = 'fdr')
names(foxFDRs) = ctypes
foxFDRs = as.data.frame(foxFDRs)
foxFDRs$CellType = rownames(foxFDRs)

toplot = toplot %>% group_by(CellType) %>%
	summarize(log10FDR = median(log10FDR), coefficient = median(coefficient)) %>%
	as.data.frame

# Add the combined p-value
toplot = merge(toplot, foxFDRs, by = 'CellType')
toplot$log10FDR_Fisher = -log10(toplot$foxFDRs)

# Add the number of HS detected per cell type
hsN = table(dars[dars$Regulation == 'Human_UP', 'CellType']) %>% as.data.frame
colnames(hsN) = c('CellType', 'Size')
toplot = merge(toplot, hsN, by = 'CellType')
toplot$CellType2 = ifelse(toplot$log10FDR_Fisher > 1.3 & toplot$coefficient > 0, toplot$CellType, '')

# Scatter plot
pdf('FOX_TFs_MOTIF_ENRICHMENT.pdf', width = 9)
ggscatter(toplot, y = 'coefficient', x = 'log10FDR_Fisher', color = 'coefficient', label = 'CellType2', size = 'Size', font.label = c(20, 'bold', 'black'), repel = T) +
  labs(x="-log10(FDR)",y="Coefficient") + scale_size_continuous(range = c(2,10)) +
  scale_color_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold'),
	axis.text.x = element_text(size=20, face = 'bold'),
	axis.text.y = element_text(size=20, face = 'bold')) +
  geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
  rotate_x_text(90)
dev.off()


# Heatmap
toplot = resDF_hup[resDF_hup$motif %in% foxsAcc,]
toplot$log10FDR = -log10(toplot$FDR)
toplot$coefficient2 = round(toplot$coefficient, digits = 2)

toplot$CellType= factor(toplot$CellType, levels = c('L2-3_1', 'L2-3_2', 'L2-3_3', 'L3-5_RORB_2', 'L3-5_RORB_1', 'L3-5_RORB_3', 'L4-5_RORB_1', 'L4-6_RORB_3', 'L4-6_RORB_2', 'L5-6_THEMIS_1', 'L5-6_THEMIS_2', 'L5-6_FEZF2_1', 'L5-6_FEZF2_2-3', 'L4-6_FEZF2'))

pdf('FOX_TFs_MOTIF_ENRICHMENT_HEATMAP.pdf', width = 14)
ggplot(toplot,aes(x = CellType, y = motif, fill = log10FDR)) +
  geom_tile(colour="white", size=0.2) +
  geom_text(aes(label=coefficient2), fontface = 'bold', size = 7) +
  labs(x="",y="") +
  scale_fill_gradient2(midpoint = 1.3, low = 'blue', high = 'red', limits = c(0,4)) +
  theme_classic() +
  theme(text = element_text(size=20, face = 'bold')) +
  rotate_x_text(90)
dev.off()


