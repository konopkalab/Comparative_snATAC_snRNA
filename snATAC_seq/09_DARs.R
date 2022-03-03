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
library(tidyverse)
library(rio)
library(dplyr)
library(data.table)
library(stringr)
source("FINAL_SCRIPTS_ATAC/custom_functions.R")
set.seed(1234)

####
## PREPARE DATA
####

# Load objects
human = readRDS("Exc_human_annotated.RDS")
chimp = readRDS("Exc_chimp_annotated.RDS")
macaque = readRDS("Exc_macaque_annotated.RDS")

# Merge and add total activity per cell
merged = merge(human, merge(chimp, macaque))
merged$tot_acc_cicero = log2(merged$nCount_RNA)

# Extract raw matrix and metadata
mat = merged@assays$peaks@counts
meta = merged[[]]

# Add sample metadata
smeta = import('ba23_atac_metadata.xlsx')
smeta$Species = NULL
meta = cbind(meta, smeta[match(meta$Sample, smeta$Sample_id), ])
meta$lib_batch = factor(meta$lib_batch)

# Order species so that it is Human, Chimp, Macaque
meta$Species = factor(meta$Species, levels = c("Human", "Chimp", "Macaque"))

####
## DETERMINE PEAKS TO BE TESTED
####

# Note that following loop creates lists that can be directly...
#... used in the da_peaksPar function for pairwise DAR analysis between species

## Human vs Chimp ##
meta_split = split(meta, meta$newannot)
reshcs = list()
reshms = list()
rescms = list()
for(i in 1:length(meta_split)){

	print(paste("Cluster", i))

	# Take only cluster cells
	clustermeta = meta_split[[i]]
	clustermat = mat[, rownames(meta_split[[i]])]

	# Check order of species. Should be Human, Chimp, Macaque
	print(names(table(clustermeta$Species)))

	# Extract cell barcodes
	cells1 = na.omit(clustermeta[clustermeta$Species == names(table(clustermeta$Species))[1], "barcode"])
	cells2 = na.omit(clustermeta[clustermeta$Species == names(table(clustermeta$Species))[2], "barcode"])
	cells3 = na.omit(clustermeta[clustermeta$Species == names(table(clustermeta$Species))[3], "barcode"])

	# How many peaks are above 85% quantile in human?
	q = rowSums(clustermat[,cells1]) / length(cells1)
	hq = quantile(q, 0.85)
	h_pkstot = sum(rowSums(clustermat[,cells1]) / length(cells1) > hq)

	# Get same number of peaks from all species
	hpeaks = (rowSums(clustermat[,cells1]) / length(cells1)) %>% sort(decreasing=T) %>% head(h_pkstot) %>% names
	cpeaks = (rowSums(clustermat[,cells2]) / length(cells2)) %>% sort(decreasing=T) %>% head(h_pkstot) %>% names
	mpeaks = (rowSums(clustermat[,cells3]) / length(cells3)) %>% sort(decreasing=T) %>% head(h_pkstot) %>% names

	# Take union of peaks from all species
	finalpks = Reduce(union, list(hpeaks,cpeaks,mpeaks))
	print(length(finalpks))

	# Filter for final peak set
	clustermat = clustermat[finalpks,]

	# Take common peaks between all species
	cmnpeaks = Reduce(intersect, list(hpeaks, cpeaks, mpeaks))

	# Mean accessibility of common peaks
	fpk1 = (rowSums(clustermat[cmnpeaks,cells1]) / length(cells1))
	fpk2 = (rowSums(clustermat[cmnpeaks,cells2]) / length(cells2))
	fpk3 = (rowSums(clustermat[cmnpeaks,cells3]) / length(cells3))

	# Ratio of mean accessibility is the scaling factor
	sf_hc = mean(fpk1) / mean(fpk2)
	sf_hm = mean(fpk1) / mean(fpk3)
	sf_cm = mean(fpk2) / mean(fpk3)

	# Create pairwise cluster matrix
	hc_clustermat = clustermat[, c(cells1, cells2)]
	hc_clustermeta = clustermeta[c(cells1, cells2),]
	hc_clustermeta$Species = droplevels(hc_clustermeta$Species)

	hm_clustermat = clustermat[, c(cells1, cells3)]
	hm_clustermeta = clustermeta[c(cells1, cells3),]
	hm_clustermeta$Species = droplevels(hm_clustermeta$Species)

	cm_clustermat = clustermat[, c(cells2, cells3)]
	cm_clustermeta = clustermeta[c(cells2, cells3),]
	cm_clustermeta$Species = droplevels(cm_clustermeta$Species)

	comp = paste(names(table(hc_clustermeta$Species)), collapse = '_')

	print(paste0(nrow(clustermat), ' peaks will be tested'))
	print(paste(length(cells1), length(cells2), ' cells for ', comp))

	# Determine covariates. The first one will be tested in log likelihood ratio test.
	covars = c("Species", "tot_acc_cicero", "sex", "human_age")

	# Save in a list to run DAR
	hc_clustermats = split(as.data.frame(hc_clustermat), rep(1:20, each = floor(nrow(hc_clustermat)/20)))
	hc_clustermats = lapply(hc_clustermats, FUN=as.matrix)
	da_package = lapply(hc_clustermats, function(x){list(x, hc_clustermeta, sf_hc, covars)})

	fn = paste0("Excitatory/DA_PACKAGES/", comp, "_", i, '.RDS')
	saveRDS(da_package, fn)

}


####
## TEST FOR DIFFERENTIAL PEAKS
####

# Parallelize the following function from bash

da_peaksPar = function(allcomb){

	# This function fits two models: with and without 'Species' covariate.
	# Then uses likelihood ratio test to compute the p-value.

	# Get variables from the list
	require(Matrix)
	clustermat = as(allcomb[[1]], 'dgCMatrix')

	clustermeta = allcomb[[2]]
	sf = allcomb[[3]]
	covars = allcomb[[4]]

	rm(allcomb)
	gc()

	# Get cells
	cells1 = clustermeta[clustermeta$Species == names(table(clustermeta$Species))[1], "barcode"]
	cells2 = clustermeta[clustermeta$Species == names(table(clustermeta$Species))[2], "barcode"]

	print(paste0("Testing: ", nrow(clustermat), " peaks"))

	# Check sample size before running
	print(table(clustermeta$Species))
	print(table(clustermeta$Sample))

	# Fit two models
	glmfit_h0 = list()
	glmfit_h1 = list()
	lr_res = list()
	acc1s = list()
	acc2s = list()
	deltas = list()

	fcs = list()
	for(i in 1:nrow(clustermat)){

		reg = rownames(clustermat)[i]

		acc1s[[i]] = mean(clustermat[i, cells1])
		acc2s[[i]] = mean(clustermat[i, cells2])

		# Calculate delta and fold change. Always normalize second group to first group.
		deltas[[i]] = mean(clustermat[i, cells1]) - (mean(clustermat[i, cells2]) * sf)
		fcs[[i]] = mean(clustermat[i, cells1]) / (mean(clustermat[i, cells2]) * sf)

		# Only total accessibility as covariate
		glmfit_h0 = glm(as.formula(paste("clustermat[reg, ] ~ ", paste(covars[2:length(covars)], collapse = '+'))), data = clustermeta, family="binomial")

		# Only species as binomial predictor in alternative model
		glmfit_h1 = glm(as.formula(paste("clustermat[reg, ] ~ ", paste(covars, collapse = '+'))), data = clustermeta, family="binomial")

		# Does adding the species result in significantly better fit?
		# Note that h0 will always give equal or worse likelihood than h1, therefore we can use lr test.
		lr_res[[i]] = lmtest::lrtest(glmfit_h0, glmfit_h1)

		if((i %% 10) == 0){print(i)}
	}

	# Create results matrix
	pval = lapply(lr_res, function(x){x[2,5]}) # Likelihood test p-value

	# Create results. Convert fold change to log2 scale for easier interpretation. Adjust p-value with BH.
	res = data.frame(pval = unlist(pval), fc = unlist(fcs),
			 acc1 = unlist(acc1s), acc2 = unlist(acc2s), delta = unlist(deltas))

	res$adj_pval = p.adjust(res$pval, method = 'BH', n = length(res$pval))
	rownames(res) = rownames(clustermat)[1:nrow(res)]

	res$condition = ifelse(res$delta > 0, "UP", "DOWN") # >0 is up in human

	rm(glmfit_h1, glmfit_h0); gc()

	return(res)
}


# Combine runs per cell type
thedirs = list.dirs(path = "Excitatory/Human_Chimp/", full = T)
thedirs = thedirs[-1]
ctypes = names(table(merged$newannot))
writedir = "Excitatory/Human_Chimp/COMBINED/"

for(i in 1:length(thedirs)){

	compfn = list.files(thedirs[i], pattern = '.RDS')
	comps = lapply(1:length(compfn), function(x){readRDS(paste0(thedirs[i], '/', compfn[x]))})
	comps = do.call(rbind, comps)
		
	comps$is_sign = ifelse(comps$adj_pval < 0.05 & (abs(comps$delta) > sd(comps$delta)*1.5), 'Sign', 'NS')

	pdf(paste0(thedirs[i], '/ScatterPlot.pdf'))
	print(ggscatter(comps, x = 'acc1', y = 'acc2', color = 'is_sign', alpha = 0.5, palette = c('grey', 'red'))+
	xlab('Chimp') + ylab('Macaque'))
	dev.off()

	comps$log10 = -log10(comps$adj_pval)
	pdf(paste0(thedirs[i], '/VolcanoPlot.pdf'))
	print(ggscatter(comps, x = 'delta', y = 'log10', color = 'is_sign', palette = c('grey', 'red'), alpha = 0.3) +
	geom_vline(xintercept = c(sd(comps$delta)*1.5, -sd(comps$delta)*1.5), linetype = 'dashed', color = 'lightblue') +
	geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'lightblue') +
	xlab('Chimp-Macaque') +
	ylab('-log10'))
	dev.off()

	# Numbers are according to the celltype order in seurat object
	ctypes_ind = gsub('.*_', '', thedirs[i]) %>% as.numeric()
	saveRDS(comps, paste0(writedir, 'Human_Chimp_', ctypes[ctypes_ind], '.RDS'))
}



# Combine runs across all cell types per pairwise comparison

# Human-Chimp
hcpath = "Excitatory/Human_Chimp/COMBINED"
hcfn = list.files(hcpath, pattern = '.RDS', full = T)

hcs = list()
for(i in 1:length(hcfn)){
	cname = sub('.*Chimp_', '', hcfn[i]) %>% sub('.RDS', '', .)
	tmp = readRDS(hcfn[i])
	tmp$cluster = cname
	tmp$peak = rownames(tmp)
	hcs[[i]] = tmp
}

reshcs = do.call(rbind, hcs)
colnames(reshcs)[1:9] = paste(colnames(reshcs)[1:9], 'hc', sep = '_')
reshcs$adj_pval_hc = p.adjust(reshcs$pval_hc, method = 'fdr')

# Human-Macaque
hmpath = "Excitatory/Human_Macaque/COMBINED"
hmfn = list.files(hmpath, pattern = '.RDS', full = T)

hms = list()
for(i in 1:length(hmfn)){
	cname = sub('.*Macaque_', '', hmfn[i]) %>% sub('.RDS', '', .)
	tmp = readRDS(hmfn[i])
	tmp$cluster = cname
	tmp$peak = rownames(tmp)
	hms[[i]] = tmp
}

reshms = do.call(rbind, hms)
colnames(reshms)[1:9] = paste(colnames(reshms)[1:9], 'hm', sep = '_')
reshms$adj_pval_hm = p.adjust(reshms$pval_hm, method = 'fdr')

# Chimp-Macaque
cmpath = "Excitatory/Chimp_Macaque/COMBINED"
cmfn = list.files(cmpath, pattern = '.RDS', full = T)

cms = list()
for(i in 1:length(cmfn)){
	cname = sub('.*Macaque_', '', cmfn[i]) %>% sub('.RDS', '', .)
	tmp = readRDS(cmfn[i])
	tmp$cluster = cname
	tmp$peak = rownames(tmp)
	cms[[i]] = tmp
}

rescms = do.call(rbind, cms)
colnames(rescms)[1:9] = paste(colnames(rescms)[1:9], 'cm', sep = '_')
rescms$adj_pval_cm = p.adjust(rescms$pval_cm, method = 'fdr')

####
## FIND SPECIES SPECIFIC PEAKS
####

# Loop over each cell type to find species specific peaks
peaklist_list = list()
ctypes = names(table(reshcs$cluster))
for(i in 1:length(ctypes)){


	# Load results
	hc = reshcs[reshcs$cluster == ctypes[i],]
	hm = reshms[reshms$cluster == ctypes[i],]
	cm = rescms[rescms$cluster == ctypes[i],]

	# Create deg result table
	alldegs = Reduce(function(x,y){
			merge(x,y, by='peak')}, list(hc,hm,cm))
	alldegs$Evolution = 'NS'
	alldegs$Regulation = 'NS'

	# Find and label human specific peaks
	hsp_up = alldegs[(alldegs$adj_pval_hc < 0.05 & alldegs$delta_hc > sd(alldegs$delta_hc)*1.5) &
			  (alldegs$adj_pval_hm < 0.05 & alldegs$delta_hm > sd(alldegs$delta_hm)*1.5), ]

	hsp_down = alldegs[(alldegs$adj_pval_hc < 0.05 & alldegs$delta_hc < sd(alldegs$delta_hc)*-1.5) &
			    (alldegs$adj_pval_hm < 0.05 & alldegs$delta_hm < sd(alldegs$delta_hm)*-1.5), ]


	# Detect ambiguous DARs
	hup_cdown = alldegs[(alldegs$adj_pval_hm < 0.05 & alldegs$delta_hm > sd(alldegs$delta_hm)*1.5) &
			    (alldegs$adj_pval_hc < 0.05 & alldegs$delta_hc > sd(alldegs$delta_hc)*1.5) &
			    (alldegs$adj_pval_cm < 0.05 & alldegs$delta_cm < sd(alldegs$delta_cm)*-1.5), ]

	hdown_cup = alldegs[(alldegs$adj_pval_hc < 0.05 & alldegs$delta_hc < sd(alldegs$delta_hc)*-1.5) &
			    (alldegs$adj_pval_hm < 0.05 & alldegs$delta_hm < sd(alldegs$delta_hm)*-1.5) &
			    (alldegs$adj_pval_cm < 0.05 & alldegs$delta_cm > sd(alldegs$delta_cm)*1.5), ]


	# Remove ambiguous DARs from human specific peaks
	hsp_up = hsp_up[!(hsp_up$peak %in% hup_cdown$peak),]
	hsp_down = hsp_down[!(hsp_down$peak %in% hdown_cup$peak),]

	# Assign decision
	alldegs[match(hsp_up$peak, alldegs$peak), 'Evolution'] = 'Human_Specific'
	alldegs[match(hsp_up$peak, alldegs$peak), 'Regulation'] = 'Human_UP'
	alldegs[match(hsp_down$peak, alldegs$peak), 'Evolution'] = 'Human_Specific'
	alldegs[match(hsp_down$peak, alldegs$peak), 'Regulation'] = 'Human_DOWN'
	alldegs[match(hup_cdown$peak, alldegs$peak), 'Evolution'] = 'Unspecific'
	alldegs[match(hup_cdown$peak, alldegs$peak), 'Regulation'] = 'H>M>C'
	alldegs[match(hdown_cup$peak, alldegs$peak), 'Evolution'] = 'Unspecific'
	alldegs[match(hdown_cup$peak, alldegs$peak), 'Regulation'] = 'H<M<C'

	# Find and label Chimp specific peaks and their mode of evolution
	csp_up = alldegs[(alldegs$adj_pval_hc < 0.05 & alldegs$delta_hc < sd(alldegs$delta_hc)*-1.5) &
			  (alldegs$adj_pval_cm < 0.05 & alldegs$delta_cm > sd(alldegs$delta_cm)*1.5), ]

	csp_down = alldegs[(alldegs$adj_pval_hc < 0.05 & alldegs$delta_hc > sd(alldegs$delta_hc)*1.5) &
			    (alldegs$adj_pval_cm < 0.05 & alldegs$delta_cm < sd(alldegs$delta_cm)*-1.5), ]

	# Remove ambiguous DARs from chimp specific peaks
	csp_up = csp_up[!(csp_up$peak %in% hdown_cup$peak),]
	csp_down = csp_down[!(csp_down$peak %in% hup_cdown$peak),]

	alldegs[match(csp_up$peak, alldegs$peak), 'Evolution'] = 'Chimp_Specific'
	alldegs[match(csp_up$peak, alldegs$peak), 'Regulation'] = 'Chimp_UP'
	alldegs[match(csp_down$peak, alldegs$peak), 'Evolution'] = 'Chimp_Specific'
	alldegs[match(csp_down$peak, alldegs$peak), 'Regulation'] = 'Chimp_DOWN'


	# Find and label DARs between Macaque and Human-Chimp
	alldegs$MvsHC = 'NS'

	m_up = alldegs[(alldegs$adj_pval_hm < 0.05 & alldegs$delta_hm < sd(alldegs$delta_hm)*-1.5) &
			(alldegs$adj_pval_cm < 0.05 & alldegs$delta_cm < sd(alldegs$delta_cm)*-1.5), ]

	m_down = alldegs[(alldegs$adj_pval_hm < 0.05 & alldegs$delta_hm > sd(alldegs$delta_hm)*1.5) &
			(alldegs$adj_pval_cm < 0.05 & alldegs$delta_cm > sd(alldegs$delta_cm)*1.5), ]

	alldegs[match(m_up$peak, alldegs$peak), 'MvsHC'] = 'Macaque_UP'
	alldegs[match(m_down$peak, alldegs$peak), 'MvsHC'] = 'Macaque_DOWN'

	peaklist_list[[i]] = alldegs
}

# Save all peak results
toplot = do.call(rbind, peaklist_list)
toplot_all = toplot
fwrite(toplot_all, 'Excitatory/Excitatory_Allpeaks_DARs.txt')

# Save only species-specific peak results
toplot = toplot[!(toplot$Evolution %in% c('NS', 'Unspecific') & toplot$MvsHC == 'NS'),]
export(toplot, 'Excitatory/Excitatory_Signpeaks_DARs.xlsx')

##
## Test human vs chimpanzee specificity with randomized background
##

toplot_all = fread('Excitatory/Excitatory_Allpeaks_DARs.txt') %>% as.data.frame
toplot = import('Excitatory/Excitatory_Signpeaks_DARs.xlsx')
ctypes = names(table(toplot$cluster))
pvalHC = list()
pvalCH = list()
pvalHM = list()
pvalMH = list()
dfs = list()
for(i in 1:length(ctypes)){

	# Extract raw results per cell type
	q = toplot_all[toplot_all$cluster == ctypes[i],]

	# Shuffle delta accessibility for background
	q$delta_hc = sample(q$delta_hc)
	q$delta_hm = sample(q$delta_hm)
	q$delta_cm = sample(q$delta_cm)

	# Number of features found per species-specific group
	hs = toplot[toplot$cluster == ctypes[i] & toplot$Evolution %in% c('Human_Specific'),] %>% nrow
	cs = toplot[toplot$cluster == ctypes[i] & toplot$Evolution %in% c('Chimp_Specific'),] %>% nrow
	ms = toplot[toplot$cluster == ctypes[i] & toplot$MvsHC != 'NS',] %>% nrow

	# Randomly select 500 (shuffled) features 10000 times...
	# ...and find species-specificity only using sign
	hsrands = list()
	csrands = list()
	msrands = list()
	hsrats = list()
	hmrats = list()
	for(j in 1:1000){
		q2 = q[sample(1:nrow(q), 500), ]
		hsrands[[j]] = sum((q2$delta_hc > 0 & q2$delta_hm > 0) | (q2$delta_hc < 0 & q2$delta_hm < 0))
		csrands[[j]] = sum((q2$delta_hc < 0 & q2$delta_cm > 0) | (q2$delta_hc > 0 & q2$delta_cm < 0))
		msrands[[j]] = sum((q2$delta_hm < 0 & q2$delta_cm < 0) | (q2$delta_hm > 0 & q2$delta_cm > 0))
		hsrats[[j]] = hsrands[[j]] / csrands[[j]]
		hmrats[[j]] = hsrands[[j]] / msrands[[j]]
	}

	# Find empirical p-value. Correct MvsHC changes (ms) by branch length based on phylogeny.
	pvalHC[[i]] = 1- (sum(hs/cs > unlist(hsrats)) / length(hsrats))
	pvalCH[[i]] = (sum(hs/cs > unlist(hsrats)) / length(hsrats))
	pvalHM[[i]] = 1- (sum(hs/(ms*6/45) > unlist(hsrats)) / length(hsrats))
	pvalMH[[i]] = (sum(hs/(ms*6/45) > unlist(hmrats)) / length(hmrats))

	# Save all results in a list
	dfRAND = data.frame(hsrat = unlist(hsrats), hmrat = unlist(hmrats), CellType = ctypes[i], type = 'RAND')
	dfOBS = data.frame(hsrat = hs/cs, hmrat = hs/(ms*6/45), CellType = ctypes[i], type = 'OBS')
	dfs[[i]] = rbind(dfRAND, dfOBS)
}

randDF = do.call(rbind, dfs)
pvalDF = data.frame(pvalHC = unlist(pvalHC), pvalCH = unlist(pvalCH), pvalHM = unlist(pvalHM), pvalMH = unlist(pvalMH), CellType = ctypes)

# Adjust p-value using FDR
pvalDF$pvalHC = p.adjust(pvalDF$pvalHC, method = 'BH')
pvalDF$pvalCH = p.adjust(pvalDF$pvalCH, method = 'BH')
pvalDF$pvalHM = p.adjust(pvalDF$pvalHM, method = 'BH')
pvalDF$pvalMH = p.adjust(pvalDF$pvalMH, method = 'BH')
finalDF = merge(randDF, pvalDF, by = 'CellType')

saveRDS(finalDF, 'Excitatory/HC_HM_RATIOS.RDS')


# Plot results #

# Human-Chimp
hcPlot = finalDF
hcPlot$is_sign = 'NS'
hcPlot[hcPlot$pvalHC < 0.05, 'is_sign'] = 'H>C'
hcPlot[hcPlot$pvalCH < 0.05, 'is_sign'] = 'H<C'

pdf('HC_Ratio_ATAC_EXC.pdf', width = 7, height = 8)
ggplot(hcPlot, aes(y = CellType, x = hsrat)) + xlim(0.5,2) +
geom_boxplot(outlier.shape = NA, color = 'black') +
geom_text(data=hcPlot[hcPlot$type == 'OBS',], aes(label=is_sign, color = is_sign), nudge_x = 0.1, fontface = 'bold', size = 7) +
scale_colour_manual(values=c('orange', 'blue', 'black')) +
theme_classic() +
theme(text = element_text(size=25, face = 'bold')) +
ylab('') + xlab('') +
rotate_x_text(90) +
geom_point(data = hcPlot[hcPlot$type == 'OBS',], color = 'red') +
theme(strip.text.x = element_text(angle = 90)) +
grids(linetype = "dashed", color = 'grey')
dev.off()


# Human-Macaque
hmPlot = finalDF
hmPlot$is_sign = 'NS'
hmPlot[hmPlot$pvalHM < 0.05, 'is_sign'] = 'H>C'
hmPlot[hmPlot$pvalMH < 0.05, 'is_sign'] = 'H<C'

pdf('HM_Ratio_ATAC_EXC.pdf', width = 7, height = 8)
ggplot(hmPlot, aes(y = CellType, x = hmrat)) + xlim(0.5,4) +
geom_boxplot(outlier.shape = NA, color = 'black') +
geom_text(data=hmPlot[hmPlot$type == 'OBS',], aes(label=is_sign, color = is_sign), nudge_x = 0.1, fontface = 'bold', size = 7) +
scale_colour_manual(values=c('orange', 'blue', 'black')) +
theme_classic() +
theme(text = element_text(size=25, face = 'bold')) +
ylab('') + xlab('') +
rotate_x_text(90) +
geom_point(data = hmPlot[hmPlot$type == 'OBS',], color = 'red') +
theme(strip.text.x = element_text(angle = 90)) +
grids(linetype = "dashed", color = 'grey')
dev.off()

####
## PLOTS
####

# Plot total number of species specific features
htot = length(unique(toplot[toplot$Evolution == 'Human_Specific', 'peak'])) # Total hspecific
ctot = length(unique(toplot[toplot$Evolution == 'Chimp_Specific', 'peak'])) # Total cspecific
mtot = length(unique(toplot[toplot$MvsHC != 'NS', 'peak'])) # Total mspecific
df = data.frame(vals = c(htot,ctot,mtot), vars = c('Human', 'Chimp', 'Macaque'))

pdf('EXC_Unique_DARs.pdf')
ggbarplot(df, x = 'vars', y = 'vals', fill = 'vars', palette = c('blue','orange','darkgreen')) +
rotate_x_text(90) +
ylab('Total Unique DARs') +
NoLegend()
dev.off()

# Dot plot of species specific peaks as validation of analysis
library(gridExtra)
randhsdf = toplot %>% filter(Evolution == 'Human_Specific') %>%
			group_by(cluster) %>% slice_sample(n=1) %>% as.data.frame
reg = randhsdf$Regulation
randhs = randhsdf$peak
cls = randhsdf$cluster

Idents(merged) = merged$Species
ggs = list()
for(i in 1:length(randhs)){
	ggs[[i]] = DotPlot(subset(merged, subset = newannot == cls[i]),
				randhs[i], assay = 'peaks', group.by = 'Species',
				scale.min = 0, scale.max = 30, col.min = 0, col.max = 0,
				cols = c("orange", "orange"), dot.scale = 20) +
				ggtitle(paste(cls[i], randhs[i], reg[i], sep = '\n'))
}

pdf("EXC_hspecific_dar_per_cluster.pdf", width = 15, height = 20)
grid.arrange(grobs = ggs, nrow = 5, ncol = 3)
dev.off()

# Plot up/down peaks per species
library(reshape2)
toplotm = melt(id.vars = c('peak', 'cluster'), measure.vars = c('MvsHC','Regulation') ,toplot)
toplotm = toplotm[!(toplotm$value == 'NS'),]

toplotm$peak = 1
toplotm = aggregate(peak~value+cluster, toplotm, FUN = sum)
toplotm$value = factor(toplotm$value,
			  levels = c('Human_UP', 'Human_DOWN', 'Chimp_UP',
					'Chimp_DOWN', 'Macaque_UP', 'Macaque_DOWN'))
toplotm$updown = factor(gsub('.*_', '', toplotm$value))

pdf('EXC_DARs_SpeciesSpecific_UPDOWN.pdf', width = 18, height = 10)
ggbarplot(toplotm, x = 'value', y = 'peak', fill = 'updown', palette = c('red', 'blue')) + facet_wrap(~cluster) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		strip.text = element_text(size = 20),
		legend.text=element_text(size=20),
		legend.title=element_text(size=20)) +
xlab('') + ylab('# of peaks') + rotate_x_text(90)
dev.off()

toplotm$peak = ifelse(toplotm$updown == 'UP', toplotm$peak, -toplotm$peak)
toplotm$Species = gsub('_.*', '', toplotm$value)
toplotm$Species = factor(toplotm$Species, levels = rev(c('Human', 'Chimp', 'Macaque')))

toplotm$cluster2 = toplotm$cluster
toplotm$cluster2 = gsub('_THEMIS', "_\nTHEMIS", toplotm$cluster2)
toplotm$cluster2 = gsub('_FEZF2', "_\nFEZF2", toplotm$cluster2)
toplotm$cluster2 = gsub('_RORB', "_\nRORB", toplotm$cluster2)
toplotm$updown = gsub('UP', "GAIN", toplotm$updown)
toplotm$updown = gsub('DOWN', "LOSS", toplotm$updown)
toplotm$updown = factor(toplotm$updown, levels = c('LOSS', 'GAIN'))

pdf('EXC_DARs_SpeciesSpecific_UPDOWN_VERTICAL.pdf', width = 10, height = 12)
ggbarplot(toplotm, x = 'Species', y = 'peak', fill = 'updown', palette = c('red', 'blue')) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		strip.text = element_text(size = 15, face = 'bold'),
		legend.text=element_text(size=20),
		legend.title=element_text(size=20)) +
xlab('') + ylab('# of peaks') + rotate_x_text(90) +
facet_wrap(~cluster2, ncol = 2, strip.position = 'right') +
coord_flip()
dev.off()



####
## DISTRIBUTION OF SUBTYPE SPECIFIC FEATURES AFTER EQUALIZING #OF FEATURES
####

toplot = import('Excitatory/Excitatory_Signpeaks_DARs.xlsx')
toplot = toplot %>% rename(peak = 'peak')
toplot_up = toplot[toplot$Regulation %in% c('Human_UP', 'Chimp_UP') |
			toplot$MvsLCA == 'Macaque_UP',]
toplot_down = toplot[toplot$Regulation %in% c('Human_DOWN', 'Chimp_DOWN') |
			toplot$MvsLCA == 'Macaque_DOWN',]

ln = length(table(toplot$cluster))

aggL = list()
sn = 500 # Number of features to subset. Results are robust to variations of this number.
for(i in 1:100){

	# Randomly subset open chromatin changes
	rand_freq_up = toplot %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Random', type = 'UP')

	# Randomly subset closed chromatin changes
	rand_freq_down = toplot %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Random', type = 'DOWN')

	colnames(rand_freq_up)[1] = 'FoundIn'
	colnames(rand_freq_down)[1] = 'FoundIn'

	rand_freq_up$peak_ratios = prop.table(rand_freq_up$Freq)
	rand_freq_down$peak_ratios = prop.table(rand_freq_down$Freq)

	rand_freq_up$change_ratios = prop.table( (rand_freq_up$Freq) * as.numeric(rand_freq_up$FoundIn) )
	rand_freq_down$change_ratios = prop.table( (rand_freq_down$Freq) * as.numeric(rand_freq_down$FoundIn) )

	# HUMAN
	hs_freq_up = toplot_up %>% filter(Regulation == 'Human_UP') %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Human', type = 'UP')

	hs_freq_down = toplot_down %>% filter(Regulation == 'Human_DOWN') %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Human', type = 'DOWN')

	colnames(hs_freq_up)[1] = 'FoundIn'
	colnames(hs_freq_down)[1] = 'FoundIn'

	hs_freq_up$peak_ratios = prop.table(hs_freq_up$Freq)
	hs_freq_down$peak_ratios = prop.table(hs_freq_down$Freq)

	hs_freq_up$change_ratios = prop.table( (hs_freq_up$Freq) * as.numeric(hs_freq_up$FoundIn) )
	hs_freq_down$change_ratios = prop.table( (hs_freq_down$Freq) * as.numeric(hs_freq_down$FoundIn) )

	# CHIMP
	cs_freq_up = toplot_up %>% filter(Regulation == 'Chimp_UP') %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Chimp', type = 'UP')

	cs_freq_down = toplot_down %>% filter(Regulation == 'Chimp_DOWN') %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Chimp', type = 'DOWN')

	colnames(cs_freq_up)[1] = 'FoundIn'
	colnames(cs_freq_down)[1] = 'FoundIn'

	cs_freq_up$peak_ratios = prop.table(cs_freq_up$Freq)
	cs_freq_down$peak_ratios = prop.table(cs_freq_down$Freq)

	cs_freq_up$change_ratios = prop.table( (cs_freq_up$Freq) * as.numeric(cs_freq_up$FoundIn) )
	cs_freq_down$change_ratios = prop.table( (cs_freq_down$Freq) * as.numeric(cs_freq_down$FoundIn) )

	# MACAQUE
	ms_freq_up = toplot_up %>% filter(MvsLCA == 'Macaque_UP') %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Macaque', type = 'UP')

	ms_freq_down = toplot_down %>% filter(MvsLCA == 'Macaque_DOWN') %>% sample_n(sn, replace = F) %>% group_by(peak) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Macaque', type = 'DOWN')

	colnames(ms_freq_up)[1] = 'FoundIn'
	colnames(ms_freq_down)[1] = 'FoundIn'

	ms_freq_up$peak_ratios = prop.table(ms_freq_up$Freq)
	ms_freq_down$peak_ratios = prop.table(ms_freq_down$Freq)

	ms_freq_up$change_ratios = prop.table( (ms_freq_up$Freq) * as.numeric(ms_freq_up$FoundIn) )
	ms_freq_down$change_ratios = prop.table( (ms_freq_down$Freq) * as.numeric(ms_freq_down$FoundIn) )

	# Compute percentage of subtype specificity
	library(RColorBrewer)
	combFreqs = rbind(rand_freq_up, rand_freq_down, hs_freq_up, hs_freq_down,
				cs_freq_up, cs_freq_down, ms_freq_up, ms_freq_down)
	combFreqs$FoundIn = as.numeric(combFreqs$FoundIn)
	combFreqs$presence = 'Unk'
	combFreqs[combFreqs$FoundIn == 1, 'presence'] = 'Specific'
	combFreqs[combFreqs$FoundIn > 1, 'presence'] = 'Shared'
	comFreqsAgg = aggregate(peak_ratios ~ presence + Species + type, combFreqs, FUN = sum)
	comFreqsAgg$presence = factor(comFreqsAgg$presence, levels = rev(c('Specific', 'Shared')))
	comFreqsAgg$Species = factor(comFreqsAgg$Species, levels = c('Random', 'Human', 'Chimp', 'Macaque'))
	comFreqsAgg$type = factor(comFreqsAgg$type, levels = c('UP', 'DOWN'))

	comFreqsAgg$Rands = i

	aggL[[i]] = comFreqsAgg
}

aggDF = do.call(rbind, aggL)

# UP
aggDF_UP = aggDF[aggDF$type == 'UP' & aggDF$presence == 'Specific',]
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'), c('Chimp', 'Macaque'))
tmpdf = aggDF_UP[aggDF_UP$Species %in% c('Human', 'Chimp', 'Macaque', 'Random'),]

pdf('EXC_ATAC_SubtypeSpecificity_UP_SPECIES.pdf')
ggboxplot(tmpdf, y = 'peak_ratios', x = 'Species', outlier.shape = NA, color = 'Species', palette = c('grey', 'blue', 'orange', 'darkgreen'), ylim = c(0.7,1.1)) + ylab('Ratio of subtype specific peaks') + xlab('') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = c(1.02,1.06,1.1)) +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

# Extract p-values
hRats = tmpdf[tmpdf$Species == 'Human', 'peak_ratios']
cRats = tmpdf[tmpdf$Species == 'Chimp', 'peak_ratios']
mRats = tmpdf[tmpdf$Species == 'Macaque', 'peak_ratios']

hcComp = wilcox.test(hRats, cRats)$p.value
hmComp = wilcox.test(hRats, mRats)$p.value
cmComp = wilcox.test(cRats, mRats)$p.value

upcomps = data.frame(pvals = c(hcComp, hmComp, cmComp), comparisons = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))

rio::export(upcomps, 'SubtypeSpecificity_Up_Comps_EXC_ATAC.xlsx')


# DOWN
aggDF_DOWN = aggDF[aggDF$type == 'DOWN' & aggDF$presence == 'Specific',]
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'), c('Chimp', 'Macaque'))
tmpdf = aggDF_DOWN[aggDF_DOWN$Species %in% c('Human', 'Chimp', 'Macaque', 'Random'),]

pdf('EXC_ATAC_SubtypeSpecificity_DOWN_SPECIES.pdf')
ggboxplot(tmpdf, y = 'peak_ratios', x = 'Species', outlier.shape = NA, color = 'Species', palette = c('grey', 'blue', 'orange', 'darkgreen'), ylim = c(0.7,1.1)) + ylab('Ratio of subtype specific peaks') + xlab('') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = c(1.02,1.06,1.1)) +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

# Extract p-values
hRats = tmpdf[tmpdf$Species == 'Human', 'peak_ratios']
cRats = tmpdf[tmpdf$Species == 'Chimp', 'peak_ratios']
mRats = tmpdf[tmpdf$Species == 'Macaque', 'peak_ratios']

hcComp = wilcox.test(hRats, cRats)$p.value
hmComp = wilcox.test(hRats, mRats)$p.value
cmComp = wilcox.test(cRats, mRats)$p.value

downcomps = data.frame(pvals = c(hcComp, hmComp, cmComp), comparisons = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))

rio::export(downcomps, 'SubtypeSpecificity_Down_Comps_EXC_ATAC.xlsx')





