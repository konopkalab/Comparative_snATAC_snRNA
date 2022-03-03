rm(list = ls())
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(Matrix)
library(ggplot2)
library(plyr)
library(Matrix.utils)
library(ggpubr)
library(biomaRt)
library(brglm)
library(lmtest)
library(tidyr)
library(rio)
library(dplyr)
library(stringr)
library('bedr')
library(data.table)
source("~/onlybiohpc/pr3/OUR_DATA/FINAL_SCRIPTS_ATAC/custom_functions.R")
set.seed(1234)

####
## LOAD DATA
####

# Load DAR results
dars = fread('COMBINED_Allpeaks_DARs.txt') %>% as.data.frame
dars$cluster = gsub('Astro', 'Astrocyte', dars$cluster)
dars$cluster = gsub('Micro', 'Microglia', dars$cluster)
dars$cluster = gsub('Oligodendrocytes', 'Oligodendrocyte', dars$cluster)

# Create peak - gene overlaps
genecoord = readRDS('Human_GENEBODY_3kbUP_All.RDS')
darsUniq = unique(dars$peak)
tmp = gsub(':|-', '_', darsUniq) %>% strsplit(., '_') %>% do.call(rbind,.) %>% as.data.frame
tmp$V2 = as.numeric(tmp$V2)
tmp$V3 = as.numeric(tmp$V3)
darsUniq = cbind(tmp, darsUniq)
darsUniq = bedr.sort.region(darsUniq, verbose = F)
genecoord = bedr.sort.region(genecoord, verbose = F)

peakOv = bedr(input = list(a = darsUniq, b = genecoord), method = "intersect", params = "-loj", verbose = F)
export(peakOv, 'peakOv.xlsx')

# Load peak-gene overlap and degs
peakOv = import('peakOv.xlsx')
degs = import('COMBINED_AllGenes_DEGs.xlsx')

# Annotate corresponding ATAC-seq cell-type
degs$atacannot = degs$cluster
mapnames = setNames(c("Astrocyte", "Astrocyte", "Astrocyte",
			rep("Oligodendrocyte", 8), rep("OPC", 3),
			"SST_CALB1+",rep("SST_CALB1-", 6),
			rep("VIP", 9),
			"Upper_Layer", "Upper_Layer", "Upper_Layer", "Upper_Layer",
			rep("PVALB_Basket", 4), "PVALB_Chandelier"),
		      c("Ast_Interlaminar", "Ast_Protoplasmic", "Ast_Mixed",
			"NFOL1", "NFOL2", "NFOL3", "MOL1", "MOL2", "MOL3", "MOL4", "MOL5", "OPC1", "OPC2", "COP",
			"SST_L1-3", "SST_L3-5", "SST_NPY_L3-6", "SST_L4-5", "SST_L4-6", "SST_L5-6_1", "SST_L5-6_2",
			"VIP_L1-2", "VIP_L1-3", "VIP_L1-3_1", "VIP_L1-3_2", "VIP_L1-4", "VIP_L2-4", "VIP_L2-5", "VIP_L2-6", "VIP_L3-6",
			"ADARB2_L1-2", "ADARB2_SYT6_L1-3", "PAX6_L1-2", "LAMP5_NMBR_L1",
			"PVALB_L2-4", "LHX6_MEPE_L4-6", "LHX6_SULF1_L2-4", "LHX6_FILIP1_L5-6", "PVALB_L2-5"))

tmp = as.character(degs$atacannot)
names(tmp) = tmp
tmp[which(tmp %in% names(mapnames))] = mapnames[tmp[which(tmp %in% names(mapnames))]] %>% as.character
degs$atacannot = factor(tmp)

# Check if all match
all( names(table(dars$cluster)) %in% names(table(degs$atacannot)) )
ctypes = names(table(dars$cluster))

###
## TEST HUMAN UP GENE-PEAK OVERLAP WITH BACKGROUND
####

peakOvGene = peakOv[peakOv$gene_name != '.',]

obsH = list()
obsC = list()
obsM = list()

pvalsH = list()
pvalsC = list()
pvalsM = list()

ratioH = list()
ratioC = list()
ratioM = list()

randomsH = list()
randomsC = list()
randomsM = list()

hTotL = list()
cTotL = list()
mTotL = list()

hOvL = list()
cOvL = list()
mOvL = list()
for(i in 1:length(ctypes)){

	peakOvGeneSUB = peakOvGene[peakOvGene$gene_name %in% degs[degs$atacannot == ctypes[i], 'Gene'],]
	darsSUB = dars[dars$cluster == ctypes[i] & dars$peak %in% peakOvGeneSUB$peak,]
	darsSUB$gene_name = peakOvGeneSUB[match(darsSUB$peak, peakOvGeneSUB$peak), 'gene_name']
	degsSUB = degs[degs$atacannot == ctypes[i],]

	darsSUB_HS_UP = darsSUB[darsSUB$Regulation == 'Human_UP',]
	degsSUB_HS_UP = degsSUB[degsSUB$Regulation == 'Human_UP',]
	
	darsSUB_CS_UP = darsSUB[darsSUB$Regulation == 'Chimp_UP',]
	degsSUB_CS_UP = degsSUB[degsSUB$Regulation == 'Chimp_UP',]

	darsSUB_MS_UP = darsSUB[darsSUB$MvsLC == 'Macaque_UP',]
	degsSUB_MS_UP = degsSUB[degsSUB$MvsLC == 'Macaque_UP',]

	hTotL[[i]] = nrow(darsSUB_HS_UP)
	cTotL[[i]] = nrow(darsSUB_CS_UP)
	mTotL[[i]] = nrow(darsSUB_MS_UP)

	hOvL[[i]] = sum(darsSUB_HS_UP$gene_name %in% degsSUB_HS_UP$Gene)
	cOvL[[i]] = sum(darsSUB_CS_UP$gene_name %in% degsSUB_CS_UP$Gene)
	mOvL[[i]] = sum(darsSUB_MS_UP$gene_name %in% degsSUB_MS_UP$Gene)

	obsH[[i]] = sum(darsSUB_HS_UP$gene_name %in% degsSUB_HS_UP$Gene) / nrow(darsSUB_HS_UP)
	obsC[[i]] = sum(darsSUB_CS_UP$gene_name %in% degsSUB_CS_UP$Gene) / nrow(darsSUB_CS_UP)
	obsM[[i]] = sum(darsSUB_MS_UP$gene_name %in% degsSUB_MS_UP$Gene) / nrow(darsSUB_MS_UP)

	# Randomize 1000 times
	randsH = list()
	randsC = list()
	randsM = list()
	for(j in 1:1000){
	
		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_HS_UP)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_HS_UP)),]
		randsH[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_CS_UP)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_CS_UP)),]
		randsC[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_MS_UP)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_MS_UP)),]
		randsM[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)
	}

	randomsH[[i]] = unlist(randsH)
	randomsC[[i]] = unlist(randsC)
	randomsM[[i]] = unlist(randsM)

	# P-value
	pvalsH[[i]] = sum(obsH[[i]] <= unlist(randsH)) / length(unlist(randsH))
	pvalsC[[i]] = sum(obsC[[i]] <= unlist(randsC)) / length(unlist(randsC))
	pvalsM[[i]] = sum(obsM[[i]] <= unlist(randsM)) / length(unlist(randsM))

	ratioH[[i]] = obsH[[i]] / mean(unlist(randsH))
	ratioC[[i]] = obsC[[i]] / mean(unlist(randsC))
	ratioM[[i]] = obsM[[i]] / mean(unlist(randsM))

	print(i)
}

pvalHUP = unlist(pvalsH)
pvalCUP = unlist(pvalsC)
pvalMUP = unlist(pvalsM)

ratioHUP = unlist(ratioH)
ratioCUP = unlist(ratioC)
ratioMUP = unlist(ratioM)

hTot = unlist(hTotL)
cTot = unlist(cTotL)
mTot = unlist(mTotL)

hOv = unlist(hOvL)
cOv = unlist(cOvL)
mOv = unlist(mOvL)

# Calculate pval with respect to background
pvalBCG = data.frame(CellType = ctypes, Pval_Human = pvalHUP, hsOpen_hsUp = hOv, hsOpen = hTot,
					Pval_Chimp = pvalCUP, csOpen_csUp = cOv, csOpen = cTot,
					Pval_Macaque = pvalMUP, msOpen_msUp = mOv, msOpen = mTot)

export(pvalBCG, 'pvalBCG.xlsx')


# Calculate pval for comparative species overlap with respect to background
fdrHC = list()
fdrHM = list()
fdrCM = list()
for(i in 1:length(hOv)){
	fdrHC[[i]] = prop.test(c(hTot[i] - hOv[i], cTot[i] - cOv[i]), c(hTot[i], cTot[i]), alternative = 'greater')$p.value
	fdrHM[[i]] = prop.test(c(hTot[i] - hOv[i], mTot[i] - mOv[i]), c(hTot[i], mTot[i]), alternative = 'greater')$p.value
	fdrCM[[i]] = prop.test(c(cTot[i] - cOv[i], mTot[i] - mOv[i]), c(cTot[i], mTot[i]), alternative = 'greater')$p.value

}
fdrHC = p.adjust(unlist(fdrHC), method = 'BH')
fdrHM = p.adjust(unlist(fdrHM), method = 'BH')
fdrCM = p.adjust(unlist(fdrCM), method = 'BH')

statsHC = data.frame(CellType = ctypes, FDR_Open_Up = fdrHC,
			hsOpen_hsUp = hOv, hsOpen = hTot, csOpen_csUp = cOv, csOpen = cTot,
			msOpen_msUp = mOv, msOpen = mTot, Comparison = 'Human-Chimp')

statsHM = data.frame(CellType = ctypes, FDR_Open_Up = fdrHM,
			hsOpen_hsUp = hOv, hsOpen = hTot, csOpen_csUp = cOv, csOpen = cTot,
			msOpen_msUp = mOv, msOpen = mTot, Comparison = 'Human-Macaque')

statsCM = data.frame(CellType = ctypes, FDR_Open_Up = fdrCM,
			hsOpen_hsUp = hOv, hsOpen = hTot, csOpen_csUp = cOv, csOpen = cTot,
			msOpen_msUp = mOv, msOpen = mTot, Comparison = 'Chimp-Macaque')

statSaveUP = rbind(statsHC, statsHM, statsCM)
export(statSaveUP, 'statSaveUP.xlsx')



## PLOT OBSERVED VS EXPECTED ##
randHDF = do.call(cbind, randomsH)
randHDF = rbind(randHDF, unlist(obsH))
colnames(randHDF) = ctypes
rownames(randHDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randHDF = melt(randHDF)
randHDF$Var1 = gsub('_.*', '', randHDF$Var1)
randHDF$Species = 'Human'

randCDF = do.call(cbind, randomsC)
randCDF = rbind(randCDF, unlist(obsC))
colnames(randCDF) = ctypes
rownames(randCDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randCDF = melt(randCDF)
randCDF$Var1 = gsub('_.*', '', randCDF$Var1)
randCDF$Species = 'Chimp'

randMDF = do.call(cbind, randomsM)
randMDF = rbind(randMDF, unlist(obsM))
colnames(randMDF) = ctypes
rownames(randMDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randMDF = melt(randMDF)
randMDF$Var1 = gsub('_.*', '', randMDF$Var1)
randMDF$Species = 'Macaque'

randAll = rbind(randHDF, randCDF, randMDF)
randAll$Species = factor(randAll$Species, levels = c('Human', 'Chimp', 'Macaque'))
randAll$Major = 'Glia'
randAll[grepl('^L[0-9]', randAll$Var2), 'Major'] = 'Excitatory'
randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Major'] = 'Inhibitory'
randAll$Var2 = factor(randAll$Var2, levels = c(as.character(unique(randAll[grepl('^L[0-9]', randAll$Var2), 'Var2'])),
				as.character(unique(randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Var2'])),
						c('Astro', 'Microglia', 'OPC', 'Oligodendrocyte')))

randEXC = randAll[randAll$Major == 'Excitatory',]
pdf('DARUP_DEGUP_EXC.pdf', width = 15, height = 8)
ggplot(randEXC, aes(x = Species, y = value)) + ylim(0,0.4) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--UP') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randEXC[randEXC$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randINH = randAll[randAll$Major == 'Inhibitory',]
pdf('DARUP_DEGUP_INH.pdf', width = 10, height = 7)
ggplot(randINH, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--UP') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randINH[randINH$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randGLI = randAll[randAll$Major == 'Glia',]
pdf('DARUP_DEGUP_GLI.pdf', width = 7, height = 7)
ggplot(randGLI, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--UP') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randGLI[randGLI$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()




## TEST SPECIES DAR-DEG OVERLAP ##

# Create data frame for all celltypes
statSaveUP$log10FDR = round(-log10(statSaveUP$FDR), digits = 2)
statSaveUP$labelFDR = ifelse(statSaveUP$log10FDR > 10, 10, statSaveUP$log10FDR )
statSaveUP$Sign = ifelse(statSaveUP$FDR_Open_Up < 0.05, 'FDR < 0.05', 'FDR > 0.05')
statSaveUP$Major = 'Glia'
statSaveUP[grepl('^L[0-9]', statSaveUP$CellType), 'Major'] = 'Excitatory'
statSaveUP[grepl('Upper|VIP|SST|PVALB|LAMP5', statSaveUP$CellType), 'Major'] = 'Inhibitory'
statSaveUP$CellType = factor(statSaveUP$CellType, levels = c(as.character(unique(statSaveUP[grepl('^L[0-9]', statSaveUP$CellType), 'CellType'])),
				as.character(unique(statSaveUP[grepl('Upper|VIP|SST|PVALB|LAMP5', statSaveUP$CellType), 'CellType'])),
						c('Astro', 'Microglia', 'OPC', 'Oligodendrocyte')))


# Plot Excitatory
statSaveUP_EXC = statSaveUP[statSaveUP$Major == 'Excitatory',]
statSaveUP_EXC$Comparison = factor(statSaveUP_EXC$Comparison, levels = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))
statSaveUP_EXC$CellType = factor(statSaveUP_EXC$CellType, levels = rev(levels(statSaveUP_EXC$CellType)))

pdf('DARUP_DEGUP_EXC_FDRs.pdf', width = 7, height = 15)
ggscatter(statSaveUP_EXC, x = 'labelFDR', y = 'CellType', color = 'Sign', palette = c('red', 'grey'), ylab = '', xlab = expression('-log'[10]*'(FDR)'), size = 5, jitter = T) +
theme(text = element_text(size=20)) +
facet_wrap(~Comparison, ncol = 1) +
rotate_x_text(90) +
geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed')
dev.off()


# Plot Inhibitory
statSaveUP_INH = statSaveUP[statSaveUP$Major == 'Inhibitory',]
statSaveUP_INH$Comparison = factor(statSaveUP_INH$Comparison, levels = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))
statSaveUP_INH$CellType = factor(statSaveUP_INH$CellType, levels = rev(levels(statSaveUP_INH$CellType)))

pdf('DARUP_DEGUP_INH_FDRs.pdf', width = 7, height = 10)
ggscatter(statSaveUP_INH, x = 'labelFDR', y = 'CellType', color = 'Sign', palette = c('red', 'grey'), ylab = '', xlab = expression('-log'[10]*'(FDR)'), size = 5, jitter = T) +
theme(text = element_text(size=20)) +
facet_wrap(~Comparison, ncol = 1) +
rotate_x_text(90) +
geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed')
dev.off()


# Plot Glia
statSaveUP_GLI = statSaveUP[statSaveUP$Major == 'Glia',]
statSaveUP_GLI$Comparison = factor(statSaveUP_GLI$Comparison, levels = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))
statSaveUP_GLI$CellType = factor(statSaveUP_GLI$CellType, levels = rev(levels(statSaveUP_GLI$CellType)))

pdf('DARUP_DEGUP_GLI_FDRs.pdf', width = 7, height = 7)
ggscatter(statSaveUP_GLI, x = 'labelFDR', y = 'CellType', color = 'Sign', palette = c('red', 'grey'), ylab = '', xlab = expression('-log'[10]*'(FDR)'), size = 5, jitter = T) +
theme(text = element_text(size=20)) +
facet_wrap(~Comparison, ncol = 1) +
rotate_x_text(90) +
geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed')
dev.off()


###
## TEST HUMAN DOWN GENE-PEAK OVERLAP WITH BACKGROUND
####

peakOvGene = peakOv[peakOv$gene_name != '.',]

obsH = list()
obsC = list()
obsM = list()

pvalsH = list()
pvalsC = list()
pvalsM = list()

ratioH = list()
ratioC = list()
ratioM = list()

randomsH = list()
randomsC = list()
randomsM = list()

hTotL = list()
cTotL = list()
mTotL = list()

hOvL = list()
cOvL = list()
mOvL = list()
for(i in 1:length(ctypes)){

	peakOvGeneSUB = peakOvGene[peakOvGene$gene_name %in% degs[degs$atacannot == ctypes[i], 'Gene'],]
	darsSUB = dars[dars$cluster == ctypes[i] & dars$peak %in% peakOvGeneSUB$peak,]
	darsSUB$gene_name = peakOvGeneSUB[match(darsSUB$peak, peakOvGeneSUB$peak), 'gene_name']
	degsSUB = degs[degs$atacannot == ctypes[i],]

	darsSUB_HS_DOWN = darsSUB[darsSUB$Regulation == 'Human_DOWN',]
	degsSUB_HS_DOWN = degsSUB[degsSUB$Regulation == 'Human_DOWN',]
	
	darsSUB_CS_DOWN = darsSUB[darsSUB$Regulation == 'Chimp_DOWN',]
	degsSUB_CS_DOWN = degsSUB[degsSUB$Regulation == 'Chimp_DOWN',]

	darsSUB_MS_DOWN = darsSUB[darsSUB$MvsLC == 'Macaque_DOWN',]
	degsSUB_MS_DOWN = degsSUB[degsSUB$MvsLC == 'Macaque_DOWN',]

	hTotL[[i]] = nrow(darsSUB_HS_DOWN)
	cTotL[[i]] = nrow(darsSUB_CS_DOWN)
	mTotL[[i]] = nrow(darsSUB_MS_DOWN)

	hOvL[[i]] = sum(darsSUB_HS_DOWN$gene_name %in% degsSUB_HS_DOWN$Gene)
	cOvL[[i]] = sum(darsSUB_CS_DOWN$gene_name %in% degsSUB_CS_DOWN$Gene)
	mOvL[[i]] = sum(darsSUB_MS_DOWN$gene_name %in% degsSUB_MS_DOWN$Gene)

	obsH[[i]] = sum(darsSUB_HS_DOWN$gene_name %in% degsSUB_HS_DOWN$Gene) / nrow(darsSUB_HS_DOWN)
	obsC[[i]] = sum(darsSUB_CS_DOWN$gene_name %in% degsSUB_CS_DOWN$Gene) / nrow(darsSUB_CS_DOWN)
	obsM[[i]] = sum(darsSUB_MS_DOWN$gene_name %in% degsSUB_MS_DOWN$Gene) / nrow(darsSUB_MS_DOWN)

	# Randomize 1000 times
	randsH = list()
	randsC = list()
	randsM = list()
	for(j in 1:1000){
	
		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_HS_DOWN)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_HS_DOWN)),]
		randsH[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_CS_DOWN)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_CS_DOWN)),]
		randsC[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_MS_DOWN)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_MS_DOWN)),]
		randsM[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)
	}

	randomsH[[i]] = unlist(randsH)
	randomsC[[i]] = unlist(randsC)
	randomsM[[i]] = unlist(randsM)

	# P-value
	pvalsH[[i]] = sum(obsH[[i]] <= unlist(randsH)) / length(unlist(randsH))
	pvalsC[[i]] = sum(obsC[[i]] <= unlist(randsC)) / length(unlist(randsC))
	pvalsM[[i]] = sum(obsM[[i]] <= unlist(randsM)) / length(unlist(randsM))

	ratioH[[i]] = obsH[[i]] / mean(unlist(randsH))
	ratioC[[i]] = obsC[[i]] / mean(unlist(randsC))
	ratioM[[i]] = obsM[[i]] / mean(unlist(randsM))

	print(i)
}

pvalHDOWN = unlist(pvalsH)
pvalCDOWN = unlist(pvalsC)
pvalMDOWN = unlist(pvalsM)

hTot = unlist(hTotL)
cTot = unlist(cTotL)
mTot = unlist(mTotL)

hOv = unlist(hOvL)
cOv = unlist(cOvL)
mOv = unlist(mOvL)

# Calculate pval with respect to background
pvalBCG = data.frame(CellType = ctypes, Pval_Human = pvalHDOWN, hsOpen_hsDown = hOv, hsOpen = hTot,
					Pval_Chimp = pvalCDOWN, csOpen_csDown = cOv, csOpen = cTot,
					Pval_Macaque = pvalMDOWN, msOpen_msDown = mOv, msOpen = mTot)

export(pvalBCG, 'pvalBCG_DOWN.xlsx')


# Calculate pval for comparative species overlap with respect to background
fdrHC = list()
fdrHM = list()
fdrCM = list()
#logoddsHC = list()
for(i in 1:length(hOv)){
	fdrHC[[i]] = prop.test(c(hTot[i] - hOv[i], cTot[i] - cOv[i]), c(hTot[i], cTot[i]), alternative = 'greater')$p.value
	fdrHM[[i]] = prop.test(c(hTot[i] - hOv[i], mTot[i] - mOv[i]), c(hTot[i], mTot[i]), alternative = 'greater')$p.value
	fdrCM[[i]] = prop.test(c(cTot[i] - cOv[i], mTot[i] - mOv[i]), c(cTot[i], mTot[i]), alternative = 'greater')$p.value

	#logoddsHC[[i]] = log2(((hOv[i]/hTot[i]) / (cOv[i]/cTot[i])))

}
fdrHC = p.adjust(unlist(fdrHC), method = 'BH')
fdrHM = p.adjust(unlist(fdrHM), method = 'BH')
fdrCM = p.adjust(unlist(fdrCM), method = 'BH')

statsHC = data.frame(CellType = ctypes, FDR_Open_Down = fdrHC,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Human-Chimp')

statsHM = data.frame(CellType = ctypes, FDR_Open_Down = fdrHM,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Human-Macaque')

statsCM = data.frame(CellType = ctypes, FDR_Open_Down = fdrCM,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Chimp-Macaque')

statSaveDOWN = rbind(statsHC, statsHM, statsCM)
export(statSaveDOWN, 'statSaveDOWN.xlsx')



## PLOT OBSERVED VS EXPECTED ##
randHDF = do.call(cbind, randomsH)
randHDF = rbind(randHDF, unlist(obsH))
colnames(randHDF) = ctypes
rownames(randHDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randHDF = melt(randHDF)
randHDF$Var1 = gsub('_.*', '', randHDF$Var1)
randHDF$Species = 'Human'

randCDF = do.call(cbind, randomsC)
randCDF = rbind(randCDF, unlist(obsC))
colnames(randCDF) = ctypes
rownames(randCDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randCDF = melt(randCDF)
randCDF$Var1 = gsub('_.*', '', randCDF$Var1)
randCDF$Species = 'Chimp'

randMDF = do.call(cbind, randomsM)
randMDF = rbind(randMDF, unlist(obsM))
colnames(randMDF) = ctypes
rownames(randMDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randMDF = melt(randMDF)
randMDF$Var1 = gsub('_.*', '', randMDF$Var1)
randMDF$Species = 'Macaque'

randAll = rbind(randHDF, randCDF, randMDF)
randAll$Species = factor(randAll$Species, levels = c('Human', 'Chimp', 'Macaque'))
randAll$Major = 'Glia'
randAll[grepl('^L[0-9]', randAll$Var2), 'Major'] = 'Excitatory'
randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Major'] = 'Inhibitory'
randAll$Var2 = factor(randAll$Var2, levels = c(as.character(unique(randAll[grepl('^L[0-9]', randAll$Var2), 'Var2'])),
				as.character(unique(randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Var2'])),
						c('Astro', 'Microglia', 'OPC', 'Oligodendrocyte')))

randEXC = randAll[randAll$Major == 'Excitatory',]
pdf('DARDOWN_DEGDOWN_EXC.pdf', width = 15, height = 8)
ggplot(randEXC, aes(x = Species, y = value)) + ylim(0,0.4) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randEXC[randEXC$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randINH = randAll[randAll$Major == 'Inhibitory',]
pdf('DARDOWN_DEGDOWN_INH.pdf', width = 10, height = 7)
ggplot(randINH, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randINH[randINH$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randGLI = randAll[randAll$Major == 'Glia',]
pdf('DARDOWN_DEGDOWN_GLI.pdf', width = 7, height = 7)
ggplot(randGLI, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randGLI[randGLI$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()




## TEST SPECIES DAR-DEG OVERLAP ##

# Create data frame for all celltypes
statSaveDOWN$log10FDR = round(-log10(statSaveDOWN$FDR), digits = 2)
statSaveDOWN$labelFDR = ifelse(statSaveDOWN$log10FDR > 10, 10, statSaveDOWN$log10FDR )
statSaveDOWN$Sign = ifelse(statSaveDOWN$FDR_Open_Down < 0.05, 'FDR < 0.05', 'FDR > 0.05')
statSaveDOWN$Major = 'Glia'
statSaveDOWN[grepl('^L[0-9]', statSaveDOWN$CellType), 'Major'] = 'Excitatory'
statSaveDOWN[grepl('Upper|VIP|SST|PVALB|LAMP5', statSaveDOWN$CellType), 'Major'] = 'Inhibitory'
statSaveDOWN$CellType = factor(statSaveDOWN$CellType, levels = c(as.character(unique(statSaveDOWN[grepl('^L[0-9]', statSaveDOWN$CellType), 'CellType'])),
				as.character(unique(statSaveDOWN[grepl('Upper|VIP|SST|PVALB|LAMP5', statSaveDOWN$CellType), 'CellType'])),
						c('Astro', 'Microglia', 'OPC', 'Oligodendrocyte')))


# Plot Excitatory
statSaveDOWN_EXC = statSaveDOWN[statSaveDOWN$Major == 'Excitatory',]
statSaveDOWN_EXC$Comparison = factor(statSaveDOWN_EXC$Comparison, levels = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))
statSaveDOWN_EXC$CellType = factor(statSaveDOWN_EXC$CellType, levels = rev(levels(statSaveDOWN_EXC$CellType)))

pdf('DARDOWN_DEGDOWN_EXC_FDRs.pdf', width = 7, height = 15)
ggscatter(statSaveDOWN_EXC, x = 'labelFDR', y = 'CellType', color = 'Sign', palette = c('red', 'grey'), ylab = '', xlab = expression('-log'[10]*'(FDR)'), size = 5, jitter = T) +
theme(text = element_text(size=20)) +
facet_wrap(~Comparison, ncol = 1) +
rotate_x_text(90) +
geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed')
dev.off()


# Plot Inhibitory
statSaveDOWN_INH = statSaveDOWN[statSaveDOWN$Major == 'Inhibitory',]
statSaveDOWN_INH$Comparison = factor(statSaveDOWN_INH$Comparison, levels = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))
statSaveDOWN_INH$CellType = factor(statSaveDOWN_INH$CellType, levels = rev(levels(statSaveDOWN_INH$CellType)))

pdf('DARDOWN_DEGDOWN_INH_FDRs.pdf', width = 7, height = 10)
ggscatter(statSaveDOWN_INH, x = 'labelFDR', y = 'CellType', color = 'Sign', palette = c('red', 'grey'), ylab = '', xlab = expression('-log'[10]*'(FDR)'), size = 5, jitter = T) +
theme(text = element_text(size=20)) +
facet_wrap(~Comparison, ncol = 1) +
rotate_x_text(90) +
geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed')
dev.off()


# Plot Glia
statSaveDOWN_GLI = statSaveDOWN[statSaveDOWN$Major == 'Glia',]
statSaveDOWN_GLI$Comparison = factor(statSaveDOWN_GLI$Comparison, levels = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))
statSaveDOWN_GLI$CellType = factor(statSaveDOWN_GLI$CellType, levels = rev(levels(statSaveDOWN_GLI$CellType)))

pdf('DARDOWN_DEGDOWN_GLI_FDRs.pdf', width = 7, height = 7)
ggscatter(statSaveDOWN_GLI, x = 'labelFDR', y = 'CellType', color = 'Sign', palette = c('red', 'grey'), ylab = '', xlab = expression('-log'[10]*'(FDR)'), size = 5, jitter = T) +
theme(text = element_text(size=20)) +
facet_wrap(~Comparison, ncol = 1) +
rotate_x_text(90) +
geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed')
dev.off()



###
## TEST HUMAN DOWN GENE - UP PEAK OVERLAP WITH BACKGROUND
####

peakOvGene = peakOv[peakOv$gene_name != '.',]

obsH = list()
obsC = list()
obsM = list()

pvalsH = list()
pvalsC = list()
pvalsM = list()

ratioH = list()
ratioC = list()
ratioM = list()

randomsH = list()
randomsC = list()
randomsM = list()

hTotL = list()
cTotL = list()
mTotL = list()

hOvL = list()
cOvL = list()
mOvL = list()
for(i in 1:length(ctypes)){

	peakOvGeneSUB = peakOvGene[peakOvGene$gene_name %in% degs[degs$atacannot == ctypes[i], 'Gene'],]
	darsSUB = dars[dars$cluster == ctypes[i] & dars$peak %in% peakOvGeneSUB$peak,]
	darsSUB$gene_name = peakOvGeneSUB[match(darsSUB$peak, peakOvGeneSUB$peak), 'gene_name']
	degsSUB = degs[degs$atacannot == ctypes[i],]

	darsSUB_HS_UP = darsSUB[darsSUB$Regulation == 'Human_UP',]
	degsSUB_HS_DOWN = degsSUB[degsSUB$Regulation == 'Human_DOWN',]
	
	darsSUB_CS_UP = darsSUB[darsSUB$Regulation == 'Chimp_UP',]
	degsSUB_CS_DOWN = degsSUB[degsSUB$Regulation == 'Chimp_DOWN',]

	darsSUB_MS_UP = darsSUB[darsSUB$MvsLC == 'Macaque_UP',]
	degsSUB_MS_DOWN = degsSUB[degsSUB$MvsLC == 'Macaque_DOWN',]

	hTotL[[i]] = nrow(darsSUB_HS_UP)
	cTotL[[i]] = nrow(darsSUB_CS_UP)
	mTotL[[i]] = nrow(darsSUB_MS_UP)

	hOvL[[i]] = sum(darsSUB_HS_UP$gene_name %in% degsSUB_HS_DOWN$Gene)
	cOvL[[i]] = sum(darsSUB_CS_UP$gene_name %in% degsSUB_CS_DOWN$Gene)
	mOvL[[i]] = sum(darsSUB_MS_UP$gene_name %in% degsSUB_MS_DOWN$Gene)

	obsH[[i]] = sum(darsSUB_HS_UP$gene_name %in% degsSUB_HS_DOWN$Gene) / nrow(darsSUB_HS_UP)
	obsC[[i]] = sum(darsSUB_CS_UP$gene_name %in% degsSUB_CS_DOWN$Gene) / nrow(darsSUB_CS_UP)
	obsM[[i]] = sum(darsSUB_MS_UP$gene_name %in% degsSUB_MS_DOWN$Gene) / nrow(darsSUB_MS_UP)

	# Randomize 1000 times
	randsH = list()
	randsC = list()
	randsM = list()
	for(j in 1:1000){
	
		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_HS_UP)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_HS_DOWN)),]
		randsH[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_CS_UP)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_CS_DOWN)),]
		randsC[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_MS_UP)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_MS_DOWN)),]
		randsM[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)
	}

	randomsH[[i]] = unlist(randsH)
	randomsC[[i]] = unlist(randsC)
	randomsM[[i]] = unlist(randsM)

	# P-value
	pvalsH[[i]] = sum(obsH[[i]] <= unlist(randsH)) / length(unlist(randsH))
	pvalsC[[i]] = sum(obsC[[i]] <= unlist(randsC)) / length(unlist(randsC))
	pvalsM[[i]] = sum(obsM[[i]] <= unlist(randsM)) / length(unlist(randsM))

	ratioH[[i]] = obsH[[i]] / mean(unlist(randsH))
	ratioC[[i]] = obsC[[i]] / mean(unlist(randsC))
	ratioM[[i]] = obsM[[i]] / mean(unlist(randsM))

	print(i)
}

pvalH_UPDOWN = unlist(pvalsH)
pvalC_UPDOWN = unlist(pvalsC)
pvalM_UPDOWN = unlist(pvalsM)


hTot = unlist(hTotL)
cTot = unlist(cTotL)
mTot = unlist(mTotL)

hOv = unlist(hOvL)
cOv = unlist(cOvL)
mOv = unlist(mOvL)

# Calculate pval with respect to background
pvalBCG = data.frame(CellType = ctypes, Pval_Human = pvalH_UPDOWN, hsOpen_hsDown = hOv, hsOpen = hTot,
					Pval_Chimp = pvalC_UPDOWN, csOpen_csDown = cOv, csOpen = cTot,
					Pval_Macaque = pvalM_UPDOWN, msOpen_msDown = mOv, msOpen = mTot)

export(pvalBCG, 'pvalBCG_UPDOWN.xlsx')


# Calculate pval for comparative species overlap with respect to background
fdrHC = list()
fdrHM = list()
fdrCM = list()
#logoddsHC = list()
for(i in 1:length(hOv)){
	fdrHC[[i]] = prop.test(c(hOv[i], cOv[i]), c(hTot[i], cTot[i]), alternative = 'less')$p.value
	fdrHM[[i]] = prop.test(c(hOv[i], mOv[i]), c(hTot[i], mTot[i]), alternative = 'less')$p.value
	fdrCM[[i]] = prop.test(c(cOv[i], mOv[i]), c(cTot[i], mTot[i]), alternative = 'less')$p.value

	#logoddsHC[[i]] = log2(((hOv[i]/hTot[i]) / (cOv[i]/cTot[i])))

}
fdrHC = p.adjust(unlist(fdrHC), method = 'BH')
fdrHM = p.adjust(unlist(fdrHM), method = 'BH')
fdrCM = p.adjust(unlist(fdrCM), method = 'BH')

statsHC = data.frame(CellType = ctypes, FDR_Open_Down = fdrHC,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Human-Chimp')

statsHM = data.frame(CellType = ctypes, FDR_Open_Down = fdrHM,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Human-Macaque')

statsCM = data.frame(CellType = ctypes, FDR_Open_Down = fdrCM,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Chimp-Macaque')

statSaveDOWN = rbind(statsHC, statsHM, statsCM)
export(statSaveDOWN, 'statSave_UPDOWN.xlsx')



## PLOT OBSERVED VS EXPECTED ##
randHDF = do.call(cbind, randomsH)
randHDF = rbind(randHDF, unlist(obsH))
colnames(randHDF) = ctypes
rownames(randHDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randHDF = melt(randHDF)
randHDF$Var1 = gsub('_.*', '', randHDF$Var1)
randHDF$Species = 'Human'

randCDF = do.call(cbind, randomsC)
randCDF = rbind(randCDF, unlist(obsC))
colnames(randCDF) = ctypes
rownames(randCDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randCDF = melt(randCDF)
randCDF$Var1 = gsub('_.*', '', randCDF$Var1)
randCDF$Species = 'Chimp'

randMDF = do.call(cbind, randomsM)
randMDF = rbind(randMDF, unlist(obsM))
colnames(randMDF) = ctypes
rownames(randMDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randMDF = melt(randMDF)
randMDF$Var1 = gsub('_.*', '', randMDF$Var1)
randMDF$Species = 'Macaque'

randAll = rbind(randHDF, randCDF, randMDF)
randAll$Species = factor(randAll$Species, levels = c('Human', 'Chimp', 'Macaque'))
randAll$Major = 'Glia'
randAll[grepl('^L[0-9]', randAll$Var2), 'Major'] = 'Excitatory'
randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Major'] = 'Inhibitory'
randAll$Var2 = factor(randAll$Var2, levels = c(as.character(unique(randAll[grepl('^L[0-9]', randAll$Var2), 'Var2'])),
				as.character(unique(randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Var2'])),
						c('Astro', 'Microglia', 'OPC', 'Oligodendrocyte')))

randEXC = randAll[randAll$Major == 'Excitatory',]
pdf('DARUP_DEGDOWN_EXC.pdf', width = 15, height = 8)
ggplot(randEXC, aes(x = Species, y = value)) + ylim(0,0.4) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randEXC[randEXC$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randINH = randAll[randAll$Major == 'Inhibitory',]
pdf('DARUP_DEGDOWN_INH.pdf', width = 10, height = 7)
ggplot(randINH, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randINH[randINH$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randGLI = randAll[randAll$Major == 'Glia',]
pdf('DARUP_DEGDOWN_GLI.pdf', width = 7, height = 7)
ggplot(randGLI, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randGLI[randGLI$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

###
## TEST HUMAN UP GENE - DOWN PEAK OVERLAP WITH BACKGROUND
####

peakOvGene = peakOv[peakOv$gene_name != '.',]

obsH = list()
obsC = list()
obsM = list()

pvalsH = list()
pvalsC = list()
pvalsM = list()

ratioH = list()
ratioC = list()
ratioM = list()

randomsH = list()
randomsC = list()
randomsM = list()

hTotL = list()
cTotL = list()
mTotL = list()

hOvL = list()
cOvL = list()
mOvL = list()
for(i in 1:length(ctypes)){

	peakOvGeneSUB = peakOvGene[peakOvGene$gene_name %in% degs[degs$atacannot == ctypes[i], 'Gene'],]
	darsSUB = dars[dars$cluster == ctypes[i] & dars$peak %in% peakOvGeneSUB$peak,]
	darsSUB$gene_name = peakOvGeneSUB[match(darsSUB$peak, peakOvGeneSUB$peak), 'gene_name']
	degsSUB = degs[degs$atacannot == ctypes[i],]

	darsSUB_HS_DOWN = darsSUB[darsSUB$Regulation == 'Human_DOWN',]
	degsSUB_HS_UP = degsSUB[degsSUB$Regulation == 'Human_UP',]
	
	darsSUB_CS_DOWN = darsSUB[darsSUB$Regulation == 'Chimp_DOWN',]
	degsSUB_CS_UP = degsSUB[degsSUB$Regulation == 'Chimp_UP',]

	darsSUB_MS_DOWN = darsSUB[darsSUB$MvsLC == 'Macaque_DOWN',]
	degsSUB_MS_UP = degsSUB[degsSUB$MvsLC == 'Macaque_UP',]

	hTotL[[i]] = nrow(darsSUB_HS_DOWN)
	cTotL[[i]] = nrow(darsSUB_CS_DOWN)
	mTotL[[i]] = nrow(darsSUB_MS_DOWN)

	hOvL[[i]] = sum(darsSUB_HS_DOWN$gene_name %in% degsSUB_HS_UP$Gene)
	cOvL[[i]] = sum(darsSUB_CS_DOWN$gene_name %in% degsSUB_CS_UP$Gene)
	mOvL[[i]] = sum(darsSUB_MS_DOWN$gene_name %in% degsSUB_MS_UP$Gene)

	obsH[[i]] = sum(darsSUB_HS_DOWN$gene_name %in% degsSUB_HS_UP$Gene) / nrow(darsSUB_HS_DOWN)
	obsC[[i]] = sum(darsSUB_CS_DOWN$gene_name %in% degsSUB_CS_UP$Gene) / nrow(darsSUB_CS_DOWN)
	obsM[[i]] = sum(darsSUB_MS_DOWN$gene_name %in% degsSUB_MS_UP$Gene) / nrow(darsSUB_MS_DOWN)

	# Randomize 1000 times
	randsH = list()
	randsC = list()
	randsM = list()
	for(j in 1:1000){
	
		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_HS_DOWN)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_HS_UP)),]
		randsH[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_CS_DOWN)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_CS_UP)),]
		randsC[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)

		darsSUB_RAND = darsSUB[sample(1:nrow(darsSUB), nrow(darsSUB_MS_DOWN)),]
		degsSUB_RAND = degsSUB[sample(1:nrow(degsSUB), nrow(degsSUB_MS_UP)),]
		randsM[[j]] = sum(darsSUB_RAND$gene_name %in% degsSUB_RAND$Gene) / nrow(darsSUB_RAND)
	}

	randomsH[[i]] = unlist(randsH)
	randomsC[[i]] = unlist(randsC)
	randomsM[[i]] = unlist(randsM)

	# P-value
	pvalsH[[i]] = sum(obsH[[i]] <= unlist(randsH)) / length(unlist(randsH))
	pvalsC[[i]] = sum(obsC[[i]] <= unlist(randsC)) / length(unlist(randsC))
	pvalsM[[i]] = sum(obsM[[i]] <= unlist(randsM)) / length(unlist(randsM))

	ratioH[[i]] = obsH[[i]] / mean(unlist(randsH))
	ratioC[[i]] = obsC[[i]] / mean(unlist(randsC))
	ratioM[[i]] = obsM[[i]] / mean(unlist(randsM))

	print(i)
}

pvalH_UPDOWN = unlist(pvalsH)
pvalC_UPDOWN = unlist(pvalsC)
pvalM_UPDOWN = unlist(pvalsM)


hTot = unlist(hTotL)
cTot = unlist(cTotL)
mTot = unlist(mTotL)

hOv = unlist(hOvL)
cOv = unlist(cOvL)
mOv = unlist(mOvL)

# Calculate pval with respect to background
pvalBCG = data.frame(CellType = ctypes, Pval_Human = pvalH_UPDOWN, hsOpen_hsDown = hOv, hsOpen = hTot,
					Pval_Chimp = pvalC_UPDOWN, csOpen_csDown = cOv, csOpen = cTot,
					Pval_Macaque = pvalM_UPDOWN, msOpen_msDown = mOv, msOpen = mTot)

export(pvalBCG, 'pvalBCG_DOWNUP.xlsx')


# Calculate pval for comparative species overlap with respect to background
fdrHC = list()
fdrHM = list()
fdrCM = list()
for(i in 1:length(hOv)){
	fdrHC[[i]] = prop.test(c(hOv[i], cOv[i]), c(hTot[i], cTot[i]), alternative = 'less')$p.value
	fdrHM[[i]] = prop.test(c(hOv[i], mOv[i]), c(hTot[i], mTot[i]), alternative = 'less')$p.value
	fdrCM[[i]] = prop.test(c(cOv[i], mOv[i]), c(cTot[i], mTot[i]), alternative = 'less')$p.value

}
fdrHC = p.adjust(unlist(fdrHC), method = 'BH')
fdrHM = p.adjust(unlist(fdrHM), method = 'BH')
fdrCM = p.adjust(unlist(fdrCM), method = 'BH')

statsHC = data.frame(CellType = ctypes, FDR_Open_Down = fdrHC,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Human-Chimp')

statsHM = data.frame(CellType = ctypes, FDR_Open_Down = fdrHM,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Human-Macaque')

statsCM = data.frame(CellType = ctypes, FDR_Open_Down = fdrCM,
			hsOpen_hsDown = hOv, hsOpen = hTot, csOpen_csDown = cOv, csOpen = cTot,
			msOpen_msDown = mOv, msOpen = mTot, Comparison = 'Chimp-Macaque')

statSaveDOWN = rbind(statsHC, statsHM, statsCM)
export(statSaveDOWN, 'statSave_DOWNUP.xlsx')

## PLOT OBSERVED VS EXPECTED ##
randHDF = do.call(cbind, randomsH)
randHDF = rbind(randHDF, unlist(obsH))
colnames(randHDF) = ctypes
rownames(randHDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randHDF = melt(randHDF)
randHDF$Var1 = gsub('_.*', '', randHDF$Var1)
randHDF$Species = 'Human'

randCDF = do.call(cbind, randomsC)
randCDF = rbind(randCDF, unlist(obsC))
colnames(randCDF) = ctypes
rownames(randCDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randCDF = melt(randCDF)
randCDF$Var1 = gsub('_.*', '', randCDF$Var1)
randCDF$Species = 'Chimp'

randMDF = do.call(cbind, randomsM)
randMDF = rbind(randMDF, unlist(obsM))
colnames(randMDF) = ctypes
rownames(randMDF) = c(paste0('Randomized_', 1:1000), 'Observed')
randMDF = melt(randMDF)
randMDF$Var1 = gsub('_.*', '', randMDF$Var1)
randMDF$Species = 'Macaque'

randAll = rbind(randHDF, randCDF, randMDF)
randAll$Species = factor(randAll$Species, levels = c('Human', 'Chimp', 'Macaque'))
randAll$Major = 'Glia'
randAll[grepl('^L[0-9]', randAll$Var2), 'Major'] = 'Excitatory'
randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Major'] = 'Inhibitory'
randAll$Var2 = factor(randAll$Var2, levels = c(as.character(unique(randAll[grepl('^L[0-9]', randAll$Var2), 'Var2'])),
				as.character(unique(randAll[grepl('Upper|VIP|SST|PVALB|LAMP5', randAll$Var2), 'Var2'])),
						c('Astro', 'Microglia', 'OPC', 'Oligodendrocyte')))

randEXC = randAll[randAll$Major == 'Excitatory',]
pdf('DARDOWN_DEGUP_EXC.pdf', width = 15, height = 8)
ggplot(randEXC, aes(x = Species, y = value)) + ylim(0,0.4) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randEXC[randEXC$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randINH = randAll[randAll$Major == 'Inhibitory',]
pdf('DARDOWN_DEGUP_INH.pdf', width = 10, height = 7)
ggplot(randINH, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randINH[randINH$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()

randGLI = randAll[randAll$Major == 'Glia',]
pdf('DARDOWN_DEGUP_GLI.pdf', width = 7, height = 7)
ggplot(randGLI, aes(x = Species, y = value)) + ylim(0,1) +
geom_boxplot(outlier.shape = NA, aes(color = Species)) +
scale_color_manual(values = c('blue', 'orange', 'darkgreen')) +
theme_classic() +
theme(text = element_text(size=20)) +
ylab('(DAR in DEG) / (DAR)') + xlab('') +
rotate_x_text(90) + ggtitle('Accessible--DOWN') +
facet_wrap(~Var2, nrow = 1) +
geom_point(data = randGLI[randGLI$Var1 == 'Observed',], color = 'red') +
theme(strip.text.x = element_text(angle = 90))
dev.off()


