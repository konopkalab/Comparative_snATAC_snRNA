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
counts <- GetAssayData(humanseur, slot = "counts")
chrobj <- CreateChromatinAssay(counts = counts, genome = "hg38", sep = c(":", "-"))

seurchrobj <- CreateSeuratObject(counts = chrobj, assay = "peaks", meta.data = humanseur[[]])
pfm <- getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))

peaks = rownames(seurchrobj) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame
peaksgr = makeGRangesFromDataFrame(df = peaks, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1')

seurchrobj <- AddMotifs(object = seurchrobj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
saveRDS(seurchrobj, 'motif_obj_human_all.RDS')

# CHIMP - CONVERT CREs TO CHIMP COORDINATES
counts <- GetAssayData(chimpseur, slot = "counts")
chimp_coords = read.table('merged_lifted_chimp.bed')
chimp_coords$V5 = sub('-', ':', chimp_coords$V5)
chimp_coords$V1 = sub('-', ':', chimp_coords$V1)
chimprows = chimp_coords[match(rownames(counts), chimp_coords$V1), 'V5']

chimp_counts = counts
rownames(chimp_counts) = chimprows
chrobj <- CreateChromatinAssay(counts = chimp_counts, sep = c(":", "-"))

seurchrobj <- CreateSeuratObject(counts = chrobj, assay = "peaks", meta.data = chimpseur[[]])
pfm <- getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))

peaks = rownames(seurchrobj) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame
peaksgr = makeGRangesFromDataFrame(df = peaks, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1')

seurchrobj <- AddMotifs(object = seurchrobj, genome = BSgenome.Ptroglodytes.UCSC.panTro5, pfm = pfm)
saveRDS(seurchrobj, 'motif_obj_chimp_all.RDS')


# MACAQUE - CONVERT CREs TO MACAQUE COORDINATES
counts <- GetAssayData(macaqueseur, slot = "counts")
macaque_coords = read.table('merged_lifted_macaque.bed')
macaque_coords$V5 = sub('-', ':', macaque_coords$V5)
macaque_coords$V1 = sub('-', ':', macaque_coords$V1)
macaquerows = macaque_coords[match(rownames(counts), macaque_coords$V1), 'V5']

macaque_counts = counts
rownames(macaque_counts) = macaquerows
chrobj <- CreateChromatinAssay(counts = macaque_counts, sep = c(":", "-"))

seurchrobj <- CreateSeuratObject(counts = chrobj, assay = "peaks", meta.data = macaqueseur[[]])
pfm <- getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))

peaks = rownames(seurchrobj) %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame
peaksgr = makeGRangesFromDataFrame(df = peaks, start.field = 'V2', end.field = 'V3', seqnames.field = 'V1')

seurchrobj <- AddMotifs(object = seurchrobj, genome = BSgenome.Mmulatta.UCSC.rheMac10, pfm = pfm)
saveRDS(seurchrobj, 'motif_obj_macaque_all.RDS')


####
## PREPARE FOR MOTIF ENRICHMENT
####

# Load motif objects
humanobj = readRDS('motif_obj_human_all.RDS')
chimpobj = readRDS('motif_obj_chimp_all.RDS')
macaqueobj = readRDS('motif_obj_macaque_all.RDS')
Idents(humanobj) = humanobj$broadannot
Idents(chimpobj) = chimpobj$broadannot
Idents(macaqueobj) = macaqueobj$broadannot

# Load background peaks
bcg = fread('COMBINED_Allpeaks_DARs.txt') %>% as.data.frame()
dars = fread('COMBINED_Signpeaks_DARs.txt') %>% as.data.frame()

# Read each species own coordinates
chimp_coords = read.table('merged_lifted_chimp.bed')
macaque_coords = read.table('merged_lifted_macaque.bed')

# Create background CRE set in each species own coordinates
bcgs_human = bcg$peak %>% unique %>% gsub(':', '-', .)
bcgs_chimp = chimp_coords[match(bcgs_human, chimp_coords$V1), 'V5'] %>% as.character()
bcgs_macaque = macaque_coords[match(bcgs_human, macaque_coords$V1), 'V5'] %>% as.character()

####
## ENRICHMENT FOR OPEN REGIONs PER SPECIES
####

## Get top 10k accessible peaks per celltype

# Human
humanTop10k = bcg %>% group_by(cluster) %>% arrange(desc(acc1_hc)) %>% slice_head(n=10000) %>% as.data.frame
accH = humanTop10k$peak %>% unique %>% gsub(':', '-', .)
hTop = FindMotifs(object = humanobj, features = accH, background = bcgs_human)
rio::export(hTop, 'hTop_AccessibleEnrich.xlsx')

# Chimp
chimpTop10k = bcg %>% group_by(cluster) %>% arrange(desc(acc2_hc)) %>% slice_head(n=10000) %>% as.data.frame
accC = chimpTop10k$peak %>% unique %>% gsub(':', '-', .)
accCConv = chimp_coords[match(accC, chimp_coords$V1), 'V5']
cTop = FindMotifs(object = chimpobj, features = accCConv, background = bcgs_chimp)
rio::export(cTop, 'cTop_AccessibleEnrich.xlsx')

# Macaque
macaqueTop10k = bcg %>% group_by(cluster) %>% arrange(desc(acc2_hm)) %>% slice_head(n=10000) %>% as.data.frame
accM = macaqueTop10k$peak %>% unique %>% gsub(':', '-', .)
accMConv = macaque_coords[match(accM, macaque_coords$V1), 'V5']
mTop = FindMotifs(object = macaqueobj, features = accMConv, background = bcgs_macaque)
rio::export(mTop, 'mTop_AccessibleEnrich.xlsx')

# Common in all species
cmnEnr = Reduce(intersect, list(hTop[hTop$fold.enrichment > 1.5, 'motif.name'], cTop[cTop$fold.enrichment > 1.5, 'motif.name'], mTop[mTop$fold.enrichment > 1.5, 'motif.name']))

saveRDS(cmnEnr, 'OpenRegionTfsCommon.RDS')

# VENN DIAGRAM
set1 <- import('hTop_AccessibleEnrich.xlsx') %>% .[.$fold.enrichment > 1.5, 'motif.name']
set2 <- import('cTop_AccessibleEnrich.xlsx') %>% .[.$fold.enrichment > 1.5, 'motif.name']
set3 <- import('mTop_AccessibleEnrich.xlsx') %>% .[.$fold.enrichment > 1.5, 'motif.name']

setlist = list(set1, set2, set3)
cols = c('blue', 'orange', 'darkgreen')
nms = c('Human', 'Chimpanzee', 'Macaque')
names(setlist) = nms

library(ggvenn)
pdf('Accessible_Motif_VENN.pdf')
ggvenn(setlist, fill_color = cols, stroke_size = 1, set_name_size = 5,text_size = 15, digits = 0, stroke_alpha = 0.5, show_percentage = F)
dev.off()

####
## ENRICHMENT FOR SPECIES-SPECIFICALLY OPEN REGIONs PER SPECIES
####

# Motif enrichment per celltype. Background is all tested peaks for given celltype
ctypes = names(table(dars$cluster))
hres = list()
cres = list()
mres = list()
for(i in 1:length(ctypes)){

	# Find DARs per cell type
	htopda = dars[dars$Regulation == 'Human_UP' & dars$cluster == ctypes[i], 'peak'] %>% gsub(':', '-', .)
	ctopda = dars[dars$Regulation == 'Chimp_UP' & dars$cluster == ctypes[i], 'peak'] %>% gsub(':', '-', .)
	mtopda = dars[dars$MvsLCA %in% c('Macaque_UP') & dars$cluster == ctypes[i] & dars$Regulation == 'NS', 'peak'] %>% gsub(':', '-', .)

	# Find the coordinates in chimpanzee and macaque
	ctopda = chimp_coords[match(ctopda, chimp_coords$V1), 'V5']
	mtopda = macaque_coords[match(mtopda, macaque_coords$V1), 'V5']

	# Find background peaks in each species own coordinates
	bcgs_human = bcg[bcg$cluster == ctypes[i], 'peak'] %>% unique %>% gsub(':', '-', .)
	bcgs_chimp = chimp_coords[match(bcgs_human, chimp_coords$V1), 'V5'] %>% as.character()
	bcgs_macaque = macaque_coords[match(bcgs_human, macaque_coords$V1), 'V5'] %>% as.character()

	# Test motif enrichments
	hres[[i]] <- FindMotifs(object = humanobj, features = htopda, background = bcgs_human)
	cres[[i]] <- FindMotifs(object = chimpobj, features = ctopda, background = bcgs_chimp)
	mres[[i]] <- FindMotifs(object = macaqueobj, features = mtopda, background = bcgs_macaque)

	hres[[i]]$FDR = p.adjust(hres[[i]]$pvalue, method = 'BH')
	cres[[i]]$FDR = p.adjust(cres[[i]]$pvalue, method = 'BH')
	mres[[i]]$FDR = p.adjust(mres[[i]]$pvalue, method = 'BH')
}
names(hres) = ctypes
names(cres) = ctypes
names(mres) = ctypes

# Export raw results
hresdf = do.call(rbind, hres)
cresdf = do.call(rbind, cres)
mresdf = do.call(rbind, mres)

hresdf$cluster = gsub('\\..*', '', rownames(hresdf))
cresdf$cluster = gsub('\\..*', '', rownames(cresdf))
mresdf$cluster = gsub('\\..*', '', rownames(mresdf))

rio::export(hresdf, 'All_Human_UP_Motif_Enrich.xlsx')
rio::export(cresdf, 'All_Chimp_UP_Motif_Enrich.xlsx')
rio::export(mresdf, 'All_Macaque_UP_Motif_Enrich.xlsx')

# Filter based on significance and effect size.
hresdf$log10FDR = -log10(hresdf$FDR)
hresdf$log10FDR = -log10(hresdf$FDR)
#fdrFilt = quantile(hresdf$log10FDR, 0.85) # 3.8759
#fcFilt = quantile(hresdf$fold.enrichment, 0.85) # 1.283969


hsign = hresdf[hresdf$log10FDR > 4 & hresdf$fold.enrichment > 1.28 & hresdf$percent.observed > 25,]

cresdf$log10FDR = -log10(cresdf$FDR)
csign = cresdf[cresdf$log10FDR > 4 & cresdf$fold.enrichment > 1.28 & cresdf$percent.observed > 25,]

mresdf$log10FDR = -log10(mresdf$FDR)
msign = mresdf[mresdf$log10FDR > 4 & mresdf$fold.enrichment > 1.28 & mresdf$percent.observed > 25,]

rio::export(hsign, 'Sign_Human_UP_Motif_Enrich.xlsx')
rio::export(csign, 'Sign_Chimp_UP_Motif_Enrich.xlsx')
rio::export(msign, 'Sign_Macaque_UP_Motif_Enrich.xlsx')

####
## ENRICHMENT FOR SPECIES-SPECIFICALLY CLOSED REGIONs PER SPECIES
####

# Motif enrichment per celltype. Background is all tested peaks for given celltype
ctypes = names(table(dars$cluster))
hres = list()
cres = list()
mres = list()
for(i in 1:length(ctypes)){

	# Find DARs per cell type
	htopda = dars[dars$Regulation == 'Human_DOWN' & dars$cluster == ctypes[i], 'peak'] %>% gsub(':', '-', .)
	ctopda = dars[dars$Regulation == 'Chimp_DOWN' & dars$cluster == ctypes[i], 'peak'] %>% gsub(':', '-', .)
	mtopda = dars[dars$MvsLCA %in% c('Macaque_DOWN') & dars$cluster == ctypes[i] & dars$Regulation == 'NS', 'peak'] %>% gsub(':', '-', .)

	# Find the coordinates in chimpanzee and macaque
	ctopda = chimp_coords[match(ctopda, chimp_coords$V1), 'V5']
	mtopda = macaque_coords[match(mtopda, macaque_coords$V1), 'V5']

	# Find background peaks in each species own coordinates
	bcgs_human = bcg[bcg$cluster == ctypes[i], 'peak'] %>% unique %>% gsub(':', '-', .)
	bcgs_chimp = chimp_coords[match(bcgs_human, chimp_coords$V1), 'V5'] %>% as.character()
	bcgs_macaque = macaque_coords[match(bcgs_human, macaque_coords$V1), 'V5'] %>% as.character()

	# Test motif enrichments
	hres[[i]] <- FindMotifs(object = humanobj, features = htopda, background = bcgs_human)
	cres[[i]] <- FindMotifs(object = chimpobj, features = ctopda, background = bcgs_chimp)
	mres[[i]] <- FindMotifs(object = macaqueobj, features = mtopda, background = bcgs_macaque)

	hres[[i]]$FDR = p.adjust(hres[[i]]$pvalue, method = 'BH')
	cres[[i]]$FDR = p.adjust(cres[[i]]$pvalue, method = 'BH')
	mres[[i]]$FDR = p.adjust(mres[[i]]$pvalue, method = 'BH')
}
names(hres) = ctypes
names(cres) = ctypes
names(mres) = ctypes

# Export raw results
hresdf = do.call(rbind, hres)
cresdf = do.call(rbind, cres)
mresdf = do.call(rbind, mres)

hresdf$cluster = gsub('\\..*', '', rownames(hresdf))
cresdf$cluster = gsub('\\..*', '', rownames(cresdf))
mresdf$cluster = gsub('\\..*', '', rownames(mresdf))

rio::export(hresdf, 'All_Human_DOWN_Motif_Enrich.xlsx')
rio::export(cresdf, 'All_Chimp_DOWN_Motif_Enrich.xlsx')
rio::export(mresdf, 'All_Macaque_DOWN_Motif_Enrich.xlsx')

# Filter based on significance and effect size.
hresdf$log10FDR = -log10(hresdf$FDR)
#fdrFilt = quantile(hresdf$log10FDR, 0.85) # 3.8759
#fcFilt = quantile(hresdf$fold.enrichment, 0.85) # 1.283969

hsign = hresdf[hresdf$log10FDR > 4 & hresdf$fold.enrichment > 1.28 & hresdf$percent.observed > 25,]

cresdf$log10FDR = -log10(cresdf$FDR)
csign = cresdf[cresdf$log10FDR > 4 & cresdf$fold.enrichment > 1.28 & cresdf$percent.observed > 25,]

mresdf$log10FDR = -log10(mresdf$FDR)
msign = mresdf[mresdf$log10FDR > 4 & mresdf$fold.enrichment > 1.28 & mresdf$percent.observed > 25,]

rio::export(hsign, 'Sign_Human_DOWN_Motif_Enrich.xlsx')
rio::export(csign, 'Sign_Chimp_DOWN_Motif_Enrich.xlsx')
rio::export(msign, 'Sign_Macaque_DOWN_Motif_Enrich.xlsx')


####
## CHECK EXPRESSION / ACCESSIBILITY OF THE ENRICHED TFs AND PLOT
####

# Load ATAC and RNA-seq and check whether TFs are accessible/transcribed
hsign = rio::import('Sign_Human_UP_Motif_Enrich.xlsx')

hatac = readRDS('Human_Merged_Filtered_Annotated_Final.RDS')
DefaultAssay(hatac) = 'RNA'
hatac = NormalizeData(hatac)

hrna = readRDS('Human_Merged_Final_Filtered_Final_Annotated.RDS')
DefaultAssay(hrna) = 'RNA'
hrna = NormalizeData(hrna)


# Annotate corresponding ATAC-seq cell-type
hrna$atacannot = hrna$newannot
mapnames = setNames(c("Astrocyte", "Astrocyte", "Astrocyte",
			rep("Oligodendrocyte", 8), rep("OPC", 2),
			"SST_CALB1+",rep("SST_CALB1-", 6),
			rep("VIP", 9),
			"Upper_Layer", "Upper_Layer", "Upper_Layer", "Upper_Layer",
			rep("PVALB_Basket", 4), "PVALB_Chandelier"),
		      c("Ast_Interlaminar", "Ast_Protoplasmic", "Ast_Mixed",
			"NFOL1", "NFOL2", "NFOL3", "MOL1", "MOL2", "MOL3", "MOL4", "MOL5", "OPC1", "OPC2",
			"SST_L1-3", "SST_L3-5", "SST_NPY_L3-6", "SST_L4-5", "SST_L4-6", "SST_L5-6_1", "SST_L5-6_2",
			"VIP_L1-2", "VIP_L1-3", "VIP_L1-3_1", "VIP_L1-3_2", "VIP_L1-4", "VIP_L2-4", "VIP_L2-5", "VIP_L2-6", "VIP_L3-5",
			"ADARB2_SYT6_L1-3", "ADARB2_L1-2", "LAMP5_NMBR_L1", "PAX6_L1-2",
			"PVALB_L2-4", "LHX6_SULF1_L2-4", "LHX6_MEPE_L4-6", "LHX6_FILIP1_L5-6", "PVALB_L2-5"))

tmp = as.character(hrna$atacannot)
names(tmp) = tmp
tmp[which(tmp %in% names(mapnames))] = mapnames[tmp[which(tmp %in% names(mapnames))]] %>% as.character
hrna$atacannot = factor(tmp)

# Check if annotations match
all( names(table(hatac$newannot)) %in% names(table(hrna$atacannot)) )

# Rename cell types for consistency
hsign$cluster = gsub('Astro', 'Astrocyte', hsign$cluster)
hsign$cluster = gsub('Micro', 'Microglia', hsign$cluster)
hsign$cluster = gsub('Oligodendrocytes', 'Oligodendrocyte', hsign$cluster)

# Find the percentage cutoff of ATAC-seq that corresponds to 10% pct.exp in RNA-seq
hsign$RNA_Perc = 'Unk'
hsign$ATAC_Perc = 'Unk'
rnaDot = DotPlot(hrna, features = rownames(hrna), group.by = 'Species')
rnaPerc = ecdf(rnaDot$data$pct.exp)
theQuant = rnaPerc(10)

atacDot = DotPlot(hatac, features = rownames(hatac), group.by = 'Species')
atacCutoff = quantile(atacDot$data$pct.exp, theQuant)

# Find the percentage accessibility / expression of each enriched TF
ctypes = names(table(hsign$cluster))
for(j in 1:length(ctypes)){

	hsubA = subset(hatac, subset = newannot == ctypes[j])
	hsubR = subset(hrna, subset = atacannot %in% ctypes[j] & Species == 'human')

	motList = hsign[hsign$cluster == ctypes[j], 'motif.name'] %>% strsplit(., '::')
	aPerc = list()
	for(i in 1:length(motList)){
		if(sum(rownames(hsubA) %in% motList[[i]]) == 0){aPerc[[i]] = 0; next}
		adat = DotPlot(hsubA, features = motList[[i]], group.by = 'broadannot')
		aPerc[[i]] = max(adat$data$pct.exp)
	}


	rPerc = list()
	for(i in 1:length(motList)){

		if(sum(rownames(hsubR) %in% motList[[i]]) == 0){rPerc[[i]] = 0; next}
		rdat = DotPlot(hsubR, features = motList[[i]], group.by = 'newannot')
		rPerc[[i]] = max(rdat$data$pct.exp)
	}

	hsign[hsign$cluster == ctypes[j], 'ATAC_Perc'] = unlist(aPerc)
	hsign[hsign$cluster == ctypes[j], 'RNA_Perc'] = unlist(rPerc)

	print(j)
}

hsign$RNA_Perc = as.numeric(hsign$RNA_Perc)
hsign$ATAC_Perc = as.numeric(hsign$ATAC_Perc)

# Assign labels to TF accessibility / expression
hsign$is_atac = ifelse(hsign$ATAC_Perc > atacCutoff, 1, 0)
hsign$is_rna = ifelse(hsign$RNA_Perc > 10, 1, 0)
hsign$which_present = ifelse(hsign$is_rna == 1 & hsign$is_atac == 1, 'RNA+ATAC',
				ifelse(hsign$is_rna == 1 & hsign$is_atac == 0, 'RNA',
				ifelse(hsign$is_rna == 0 & hsign$is_atac == 1, 'ONLY_ATAC', 'NEITHER')))

# Assign major cell types
library(ggrepel)
hsignSub = hsign[hsign$which_present %in% c('RNA+ATAC', 'ONLY_ATAC'),]
hsignSub$which_present = factor(hsignSub$which_present, levels = c('RNA+ATAC', 'ONLY_ATAC'))
hsignSub$Major = 'Glia'
hsignSub[grepl('^L[0-9]', hsignSub$cluster), 'Major'] = 'Excitatory'
hsignSub[grepl('SST|PVALB|VIP|LAMP5', hsignSub$cluster), 'Major'] = 'Inhibitory'


# Plot per major cell type

pdf('EXC_MotifEnrich_humanUP.pdf', width = 20, height = 5)
ggscatter(hsignSub[hsignSub$Major == 'Excitatory', ], x = 'log10FDR', y = 'fold.enrichment', color = 'which_present') +
geom_text_repel(aes(log10FDR, fold.enrichment, label = motif.name), fontface = 'bold', size = 5) +
facet_grid(which_present ~ cluster) +
theme_bw() +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90)
dev.off()


pdf('INH_MotifEnrich_humanUP.pdf', width = 10, height = 5)
ggscatter(hsignSub[hsignSub$Major == 'Inhibitory', ], x = 'log10FDR', y = 'fold.enrichment', color = 'which_present') +
geom_text_repel(aes(log10FDR, fold.enrichment, label = motif.name), fontface = 'bold', size = 5) +
facet_grid(which_present ~ cluster) +
theme_bw() +
theme(text = element_text(size=20)) +
rotate_x_text(90)
dev.off()


pdf('GLIA_MotifEnrich_humanUP.pdf', width = 12, height = 5)
ggscatter(hsignSub[hsignSub$Major == 'Glia', ], x = 'log10FDR', y = 'fold.enrichment', color = 'which_present') +
geom_text_repel(aes(log10FDR, fold.enrichment, label = motif.name), fontface = 'bold', size = 5) +
facet_grid(which_present ~ cluster) +
theme_bw() +
theme(text = element_text(size=20)) +
rotate_x_text(90)
dev.off()


rio::export(hsign, 'Sign_Human_UP_Motif_Enrich_PERCENTAGES.xlsx')



####
## COMPARE IEG-TFs BETWEN HUMAN-CHIMP-MACAQUE -- UP
####

hALL = rio::import('All_Human_UP_Motif_Enrich.xlsx')
cALL = rio::import('All_Chimp_UP_Motif_Enrich.xlsx')
mALL = rio::import('All_Macaque_DOWN_Motif_Enrich.xlsx')

# Keep IEGs in mouse single-cell
hrvatin = rio::import('hrvatin_2018/ADG_IEG.xlsx')

hrvatinEXC = hrvatin[, grepl('Name|Exc', colnames(hrvatin))]
hrvatinINH = hrvatin[, grepl('Name|Int', colnames(hrvatin))]
hrvatinGLIA = hrvatin[, grepl('Name|Olig|OPC|Astro|Micro', colnames(hrvatin))]
hrvatinL = list(hrvatinEXC, hrvatinINH, hrvatinGLIA)

iegctypes = c('Exc', 'Inh', 'Glia')
iegL = list()
for(i in 1:3){

	hrvatinL[[i]][is.na(hrvatinL[[i]])] = NA
	iegs = list()
	for(j in 1:nrow(hrvatinL[[i]])){

		tmp = hrvatinL[[i]][j,2:(ncol(hrvatinL[[i]])-3)]
		tmp[is.na(tmp)] = 'NotFound'
		
		if(sum(tmp == 'a') > 0){iegs[[j]] = hrvatinL[[i]][j,]}
	}

	iegL[[i]] = do.call(rbind, iegs)
}

# Convert mouse to human
iegsEXC_H = convertMouseGeneList(unique(iegL[[1]][,1]))
iegsINH_H = convertMouseGeneList(unique(iegL[[2]][,1]))
iegsGLIA_H = convertMouseGeneList(unique(iegL[[3]][,1]))

saveRDS(iegsEXC_H, 'hrvatin_2018/IEGs_EXC_HUMAN.RDS')
saveRDS(iegsINH_H, 'hrvatin_2018/IEGs_INH_HUMAN.RDS')
saveRDS(iegsGLIA_H, 'hrvatin_2018/IEGs_GLIA_HUMAN.RDS')

# Plot EXC
iegsEXC_H = readRDS('hrvatin_2018/IEGs_EXC_HUMAN.RDS')
hALL_CONST = hALL[grepl('MEF2|CREB|SRF', hALL$motif.name) ,] %>% add_column(tftype = 'CONSTITUTIVE_TFs')
cALL_CONST = cALL[grepl('MEF2|CREB|SRF', cALL$motif.name) ,] %>% add_column(tftype = 'CONSTITUTIVE_TFs')

hALL_IEG = hALL[grepl(paste(iegsEXC_H, collapse = '|'), hALL$motif.name),] %>% add_column(tftype = 'IEG_TFs')
cALL_IEG = cALL[grepl(paste(iegsEXC_H, collapse = '|'), cALL$motif.name),] %>% add_column(tftype = 'IEG_TFs')
toplotH = rbind(hALL_CONST, hALL_IEG)
toplotC = rbind(cALL_CONST, cALL_IEG)
toplotH = toplotH[grepl('^L[0-9]', toplotH$cluster),]
toplotC = toplotC[grepl('^L[0-9]', toplotC$cluster),]
toplotH$type = 'Human'
toplotC$type = 'Chimp'

toplot = rbind(toplotH, toplotC)
comps = list(c('Human', 'Chimp'))

pdf('HumanVsChimp_UP_IEG_EXC.pdf')
ggboxplot(toplot, x = 'type', y = 'fold.enrichment', color = 'type', outlier.shape = NA,  palette = c('blue', 'orange', 'darkgreen'), ylim = c(0.5,2)) +
geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = 1.7) +
xlab('') + ylab('Fold Enrichment') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
NoLegend()
dev.off()

toplot$log10FDR = -log10(toplot$FDR)
toplot$type = factor(toplot$type, levels = c('Human', 'Chimp'))

toplotIEG = toplot[toplot$tftype != 'CONSTITUTIVE_TFs',]
toplotCONST = toplot[toplot$tftype == 'CONSTITUTIVE_TFs',]

pdf('EXC_HS_CS_UP_CRE_ENRICH_CONST_AND_IEG.pdf', width = 12, height = 5)
ggscatter(toplot, x = 'log10FDR', y = 'fold.enrichment', size = 5, color = 'grey',  xlim = c(0,30), ylim = c(0,2), alpha = 0.3) +
theme(text=element_text(size=30, face = 'bold')) + ylab('Fold Enrichment') +
xlab(expression('-log'[10]*'(FDR)')) +
geom_vline(xintercept = -log10(0.0001), linetype = 'dashed', color = 'red') + NoLegend() +
geom_text_repel(data = toplotCONST[order(toplotCONST$log10FDR,decreasing=T),][1:4,], aes(label = motif.name), nudge_y = 0.02, nudge_x = -0.05, fontface = 'bold', size = 5) +
geom_point(data = toplotCONST, aes(label = motif.name, color = type), size = 3) +
scale_colour_manual(values = c('blue', 'orange')) +
facet_wrap(~type)
dev.off()


# Plot INH
iegsINH_H = readRDS('hrvatin_2018/IEGs_INH_HUMAN.RDS')
hALL_CONST = hALL[grepl('MEF2|CREB|SRF', hALL$motif.name) ,] %>% add_column(tftype = 'CONSTITUTIVE_TFs')
cALL_CONST = cALL[grepl('MEF2|CREB|SRF', cALL$motif.name) ,] %>% add_column(tftype = 'CONSTITUTIVE_TFs')

hALL_IEG = hALL[grepl(paste(iegsINH_H, collapse = '|'), hALL$motif.name),] %>% add_column(tftype = 'IEG_TFs')
cALL_IEG = cALL[grepl(paste(iegsINH_H, collapse = '|'), cALL$motif.name),] %>% add_column(tftype = 'IEG_TFs')
toplotH = rbind(hALL_CONST, hALL_IEG)
toplotC = rbind(cALL_CONST, cALL_IEG)
toplotH = toplotH[grepl('Upper|PVALB|SST|LAMP5|VIP', toplotH$cluster),]
toplotC = toplotC[grepl('Upper|PVALB|SST|LAMP5|VIP', toplotC$cluster),]
toplotH$type = 'Human'
toplotC$type = 'Chimp'

toplot = rbind(toplotH, toplotC)
comps = list(c('Human', 'Chimp'))

pdf('HumanVsChimp_UP_IEG_INH.pdf')
ggboxplot(toplot, x = 'type', y = 'fold.enrichment', color = 'type', outlier.shape = NA,  palette = c('blue', 'orange', 'darkgreen'), ylim = c(0.5,2)) +
geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = 1.7) +
xlab('') + ylab('Fold Enrichment') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
NoLegend()
dev.off()


# Plot Constitutive vs IEG
toplot$log10FDR = -log10(toplot$FDR)
toplot$type = factor(toplot$type, levels = c('Human', 'Chimp'))

toplotIEG = toplot[toplot$tftype != 'CONSTITUTIVE_TFs',]
toplotCONST = toplot[toplot$tftype == 'CONSTITUTIVE_TFs',]

pdf('INH_HS_CS_UP_CRE_ENRICH_CONST_AND_IEG.pdf', width = 12, height = 5)
ggscatter(toplot, x = 'log10FDR', y = 'fold.enrichment', size = 5, color = 'grey',  xlim = c(0,30), ylim = c(0,2), alpha = 0.3) +
theme(text=element_text(size=30, face = 'bold')) + ylab('Fold Enrichment') +
xlab(expression('-log'[10]*'(FDR)')) +
geom_vline(xintercept = -log10(0.0001), linetype = 'dashed', color = 'red') + NoLegend() +
geom_text_repel(data = toplotCONST[order(toplotCONST$log10FDR,decreasing=T),][1:2,], aes(label = motif.name), nudge_y = 0.02, nudge_x = -0.05, fontface = 'bold', size = 5) +
geom_point(data = toplotCONST, aes(label = motif.name, color = type), size = 3) +
scale_colour_manual(values = c('blue', 'orange')) +
facet_wrap(~type)
dev.off()

####
## COMPARE CONSERVED TFs BETWEN HUMAN-CHIMP -- UP
####

cmns = readRDS('OpenRegionTfsCommon.RDS')

hCons = hALL[hALL$motif.name %in% cmns, ]
cCons = cALL[cALL$motif.name %in% cmns, ]

toplotH = hCons
toplotC = cCons

toplotH$type = 'Human'
toplotC$type = 'Chimp'

toplot = rbind(toplotH, toplotC)
comps = list(c('Human', 'Chimp'))

pdf('HumanVsChimp_UP_CONS_ALL.pdf')
ggboxplot(toplot, x = 'type', y = 'fold.enrichment', color = 'type', outlier.shape = NA, palette = c('blue', 'orange', 'darkgreen'), ylim = c(0.3,2)) +
geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = 1.7) +
xlab('') + ylab('Fold Enrichment') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
facet_wrap(~cluster, ncol = 8)
dev.off()

# EXC CONSERVED
toplotH = hCons[grepl('^L[0-9]', hCons$cluster),]
toplotC = cCons[grepl('^L[0-9]', cCons$cluster),]
toplotH$type = 'Human'
toplotC$type = 'Chimp'

toplot = rbind(toplotH, toplotC)
comps = list(c('Human', 'Chimp'))
pdf('HumanVsChimp_UP_CONS_EXC.pdf')
ggboxplot(toplot, x = 'type', y = 'fold.enrichment', color = 'type', outlier.shape = NA,  palette = c('blue', 'orange', 'darkgreen'), ylim = c(0.5,2)) +
geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = 1.7) +
xlab('') + ylab('Fold Enrichment') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
NoLegend()
dev.off()

# INH CONSERVED
toplotH = hCons[grepl('Upper|PVALB|SST|LAMP5|VIP', hCons$cluster),]
toplotC = cCons[grepl('Upper|PVALB|SST|LAMP5|VIP', cCons$cluster),]
toplotH$type = 'Human'
toplotC$type = 'Chimp'

toplot = rbind(toplotH, toplotC)
comps = list(c('Human', 'Chimp'))
pdf('HumanVsChimp_UP_CONS_INH.pdf')
ggboxplot(toplot, x = 'type', y = 'fold.enrichment', color = 'type', outlier.shape = NA,  palette = c('blue', 'orange', 'darkgreen'), ylim = c(0.5,2)) +
geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = 1.7) +
xlab('') + ylab('Fold Enrichment') +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
NoLegend()
dev.off()

####
## PLOT ENRICHED TF FAMILIES - HS
####

hsign = rio::import('Sign_Human_UP_Motif_Enrich_PERCENTAGES.xlsx')
pfm <- getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))


# Assign motif family
dflist = list()
for(i in 1:length(pfm)){
	motif.name = pfm[[i]]@name
	fm = pfm[[i]]@tags$family
	if(length(fm) == 0){fm = motif.name}
	dflist[[i]] = cbind(fm,motif.name) %>% as.data.frame()
}
motdf = do.call(rbind, dflist)
motdf = motdf[!(duplicated(motdf$motif.name)),]
hsign2 = merge(hsign, motdf, by = 'motif.name')
hsign2$log10FDR = -log10(hsign2$FDR)
keepFamilies = table(hsign2$fm) %>% sort %>% tail(5) %>% names
hsign2Top = hsign2[hsign2$fm %in% keepFamilies,]
hsign2Top$fm = factor(hsign2Top$fm, levels = keepFamilies)

hsign2Top$major = 'Other'
hsign2Top[grepl('^L[0-9]', hsign2Top$cluster), 'major'] = 'Excitatory'
hsign2Top[grepl('PVALB|SST|VIP|Upper|LAMP5', hsign2Top$cluster), 'major'] = 'Inhibitory'
hsign2Top[grepl('Ast', hsign2Top$cluster), 'major'] = 'Other'
hsign2Top[grepl('Mic', hsign2Top$cluster), 'major'] = 'Other'
hsign2Top[grepl('OPC', hsign2Top$cluster), 'major'] = 'OPC'
hsign2Top[grepl('Oli', hsign2Top$cluster), 'major'] = 'Oligodendrocyte'

hsign2Top$fm = gsub('homeo ', "homeo\n", hsign2Top$fm)
hsign2Top$fm = gsub('zinc ', "zinc\n", hsign2Top$fm)
hsign2Top$fm = gsub(' factors', "\nfactors", hsign2Top$fm)

pdf('MotifEnrich_AcrossClasses_HSIGN.pdf', width = 17, height = 8)
ggscatter(hsign2Top, y = 'fold.enrichment', x = 'log10FDR', color = 'major', alpha = 0.5, size = 4) +
facet_wrap(~fm, nrow = 2) + rotate_x_text(90) +
scale_color_manual(values = c('blue', 'red', 'darkgreen', 'green')) +
scale_x_continuous(breaks = seq(5,30,5)) +
geom_vline(xintercept = 4, linetype = 'dashed', color = 'red') +
ylab('Fold Enrichment') + xlab(expression('-log'[10]*'(FDR)')) +
theme(text = element_text(size=20, face = 'bold'))
dev.off()


# What % of enrichments are explained by these TFs?
top10Rat = nrow(hsign2Top)/nrow(hsign2)
theRest = 1-top10Rat

toplot = data.frame(vals = c(top10Rat, theRest), vars = c('Top 5', 'Others'))
toplot$vars = factor(toplot$vars, levels = c('Top 5', 'Others'))
pdf('MotifEnrich_Top5.pdf', width = 5, height = 2)
ggbarplot(toplot, x = 'vars', y = 'vals', color = 'grey', fill = 'grey') +
theme(text = element_text(size=20, face = 'bold')) +
xlab('') + ylab('Ratio') +
scale_y_continuous(breaks = seq(0,1,0.1)) +
coord_flip()
dev.off()



####
## PLOT ENRICHED TF FAMILIES - CONSERVED
####

cmnEnr = readRDS('OpenRegionTfsCommon.RDS')
hTop = rio::import('hTop_AccessibleEnrich.xlsx')
pfm <- getMatrixSet(x = JASPAR2018, opts = list(species = 9606, all_versions = FALSE))


# Assign motif family
dflist = list()
for(i in 1:length(pfm)){
	motif.name = pfm[[i]]@name
	fm = pfm[[i]]@tags$family
	if(length(fm) == 0){fm = motif.name}
	dflist[[i]] = cbind(fm,motif.name) %>% as.data.frame()
}
motdf = do.call(rbind, dflist)
motdf = motdf[!(duplicated(motdf$motif.name)),]
hTop2 = merge(hTop, motdf, by = 'motif.name')
hTop2$log10FDR = -log10(p.adjust(hTop$pvalue, method = 'BH'))
hTop2$log10FDR = ifelse(hTop2$log10FDR >30, 30, hTop2$log10FDR)
hTop2 = hTop2[hTop2$motif.name %in% cmnEnr,]
keepFamilies = table(hTop2$fm) %>% sort %>% tail(5) %>% names
hTop2Top = hTop2[hTop2$fm %in% keepFamilies,]
hTop2Top$fm = factor(hTop2Top$fm, levels = keepFamilies)

hTop2Top$fm = gsub('3 ', "3\n", hTop2Top$fm)
hTop2Top$fm = gsub('finger ', "finger\n", hTop2Top$fm)
hTop2Top$fm = gsub('related ', "related\n", hTop2Top$fm)
selectMot = hTop2Top %>% group_by(fm) %>% arrange(desc(fold.enrichment)) %>%
		slice_head(n=3) %>% dplyr::select(motif.name) %>% as.data.frame %>% .[,2]
hTop2Top$Label = ifelse(hTop2Top$motif.name %in% selectMot, 1, 0)

pdf('MotifEnrich_AcrossClasses_HTOP.pdf', width = 14, height = 8)
ggscatter(hTop2Top, y = 'fold.enrichment', x = 'fm', color = 'fm', alpha = 0.5, size = 4) +
ylab('Fold Enrichment') + xlab('') +
theme(text = element_text(size=20, face = 'bold')) +
geom_text(data = hTop2Top[hTop2Top$Label == 1,], aes(y=fold.enrichment, x=fm, label = motif.name), nudge_x = 0.1, size = 6, fontface = 'bold') +
NoLegend()
dev.off()

####
## EXPORT FINAL TABLE
####

hsign = rio::import('Sign_Human_UP_Motif_Enrich_PERCENTAGES.xlsx')
hsign = hsign %>% dplyr::rename(Activity_of_TF = 'which_present', CellType = 'cluster') %>%
		dplyr::select(-log10FDR, -ATAC_Perc, -RNA_Perc, -is_rna, -is_atac)
hsign[grepl('Astro', hsign$CellType), 'CellType'] = 'Astrocyte'
hsign[grepl('Micro', hsign$CellType), 'CellType'] = 'Microglia'
hsign[grepl('Oli', hsign$CellType), 'CellType'] = 'Oligodendrocyte'
hsign$TF_Family = motdf[match(hsign$motif.name, motdf$motif.name), 'fm'] %>% droplevels

csign = rio::import('Sign_Chimp_UP_Motif_Enrich.xlsx')
csign = csign %>% dplyr::rename(CellType = 'cluster') %>%
		dplyr::select(-log10FDR)
csign[grepl('Astro', csign$CellType), 'CellType'] = 'Astrocyte'
csign[grepl('Micro', csign$CellType), 'CellType'] = 'Microglia'
csign[grepl('Oli', csign$CellType), 'CellType'] = 'Oligodendrocyte'
csign$TF_Family = motdf[match(csign$motif.name, motdf$motif.name), 'fm'] %>% droplevels

msign = rio::import('Sign_Macaque_UP_Motif_Enrich.xlsx')
msign = msign %>% dplyr::rename(CellType = 'cluster') %>%
		dplyr::select(-log10FDR)
msign[grepl('Astro', msign$CellType), 'CellType'] = 'Astrocyte'
msign[grepl('Micro', msign$CellType), 'CellType'] = 'Microglia'
msign[grepl('Oli', msign$CellType), 'CellType'] = 'Oligodendrocyte'
msign$TF_Family = motdf[match(msign$motif.name, motdf$motif.name), 'fm'] %>% droplevels

hAll = rio::import('All_Human_UP_Motif_Enrich.xlsx')
hAll = hAll %>% dplyr::rename(CellType = 'cluster')
hAll[grepl('Astro', hAll$CellType), 'CellType'] = 'Astrocyte'
hAll[grepl('Micro', hAll$CellType), 'CellType'] = 'Microglia'
hAll[grepl('Oli', hAll$CellType), 'CellType'] = 'Oligodendrocyte'
hAll$TF_Family = motdf[match(hAll$motif.name, motdf$motif.name), 'fm'] %>% droplevels

cAll = rio::import('All_Chimp_UP_Motif_Enrich.xlsx')
cAll = cAll %>% dplyr::rename(CellType = 'cluster')
cAll[grepl('Astro', cAll$CellType), 'CellType'] = 'Astrocyte'
cAll[grepl('Micro', cAll$CellType), 'CellType'] = 'Microglia'
cAll[grepl('Oli', cAll$CellType), 'CellType'] = 'Oligodendrocyte'
cAll$TF_Family = motdf[match(cAll$motif.name, motdf$motif.name), 'fm'] %>% droplevels

mAll = rio::import('All_Macaque_UP_Motif_Enrich.xlsx')
mAll = mAll %>% dplyr::rename(CellType = 'cluster')
mAll[grepl('Astro', mAll$CellType), 'CellType'] = 'Astrocyte'
mAll[grepl('Micro', mAll$CellType), 'CellType'] = 'Microglia'
mAll[grepl('Oli', mAll$CellType), 'CellType'] = 'Oligodendrocyte'
mAll$TF_Family = motdf[match(mAll$motif.name, motdf$motif.name), 'fm'] %>% droplevels

toExport = list(hsign, csign, msign, hAll, cAll, mAll)
names(toExport) = c('HS_Open_Sign_Motif_Enrichment', 'CS_Open_Sign_Motif_Enrichment',
			'MvsHC_Open_Sign_Motif_Enrichment', 'HS_Open_All_Motif_Enrichment',
			'CS_Open_All_Motif_Enrichment', 'MvsHC_Open_All_Motif_Enrichment')

rio::export(toExport, 'STable_MotifEnrich.xlsx')

