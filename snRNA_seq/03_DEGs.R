rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(Matrix.utils)
library(ggpubr)
library(reshape2)
library(data.table)
library(rio)
library(scater)
library(MAST)
library(gridExtra)
library(pheatmap)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

####
## DEG ANALYSIS EXCITATORY
####

# Load RNA
allrna = readRDS('EXC_Integrated_Annotated.RDS')

metakeep = c('newannot', 'orig.ident', 'Species', 'seurat_clusters', 'nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT')
allrna@meta.data = allrna[[metakeep]]

# Add sample metadata
smeta = import('ba23_rna_metadata.xlsx')
smeta$Species = NULL
meta = allrna[[]]
meta = cbind(meta, smeta[match(meta$orig.ident, smeta$Sample_id), ])
meta$lib_batch = factor(meta$lib_batch)

# Order species so that it is Human, Chimp, Macaque
meta$Species = factor(meta$Species, levels = c("human", "chimp", "macaque"))
allrna@meta.data = meta

# Normalize and log transform RNA
allrna = NormalizeData(allrna, normalization.method = "LogNormalize", scale.factor = 10000, assay = 'RNA')

# Species
sps = c('human', 'chimp', 'macaque')

# Extract cell types
ctypes = table(allrna$newannot) %>% names
Idents(allrna) = allrna$Species

# Set comparisons
comp_hc = c('human', 'chimp')
comp_hm = c('human', 'macaque')
comp_cm = c('chimp', 'macaque')

# Find DEGs between pairwise comparison. Function uses MAST
find_degs(allrna, comp_hc, ctypes, qNorm = T)
find_degs(allrna, comp_hm, ctypes, qNorm = T)
find_degs(allrna, comp_cm, ctypes, qNorm = T)

# Load and find species specific genes
genelist_list = list()
for(i in 1:length(ctypes)){

	# Load DEG results
	hc = list.files(pattern = paste0('human_chimp_', ctypes[i], '_DEGs.RDS'),
			path = 'DEG_RESULTS/', full = T) %>% readRDS
	colnames(hc) = paste(colnames(hc), 'hc', sep = '_')
	hc$Gene = rownames(hc)

	hm = list.files(pattern = paste0('human_macaque_', ctypes[i], '_DEGs.RDS'),
			path = 'DEG_RESULTS/', full = T) %>% readRDS
	colnames(hm) = paste(colnames(hm), 'hm', sep = '_')
	hm$Gene = rownames(hm)

	cm = list.files(pattern = paste0('chimp_macaque_', ctypes[i], '_DEGs.RDS'),
			path = 'DEG_RESULTS/', full = T) %>% readRDS
	colnames(cm) = paste(colnames(cm), 'cm', sep = '_')
	cm$Gene = rownames(cm)	

	# Create deg result table
	alldegs = Reduce(function(x,y){
			merge(x,y, by='Gene')}, list(hc,hm,cm))
	alldegs$Evolution = 'NS'
	alldegs$Regulation = 'NS'

	fcCutOff = 0.25

	# Find and label human specific genes
	hsp_up = alldegs[(alldegs$adj_p_value_hc < 0.05 & alldegs$avg_logfc_hc > fcCutOff) &
			  (alldegs$adj_p_value_hm < 0.05 & alldegs$avg_logfc_hm > fcCutOff),]

	hsp_down = alldegs[(alldegs$adj_p_value_hc < 0.05 & alldegs$avg_logfc_hc < -fcCutOff) &
			    (alldegs$adj_p_value_hm < 0.05 & alldegs$avg_logfc_hm < -fcCutOff),]

	# Detect ambiguous DEGs
	hdown_cup = alldegs[(alldegs$adj_p_value_hm < 0.05 & alldegs$avg_logfc_hm < -fcCutOff) &
			    (alldegs$adj_p_value_hc < 0.05 & alldegs$avg_logfc_hc < -fcCutOff) &
			    (alldegs$adj_p_value_cm < 0.05 & alldegs$avg_logfc_cm > fcCutOff), ]

	hup_cdown = alldegs[(alldegs$adj_p_value_hm < 0.05 & alldegs$avg_logfc_hm > fcCutOff) &
			    (alldegs$adj_p_value_hc < 0.05 & alldegs$avg_logfc_hc > fcCutOff) &
			    (alldegs$adj_p_value_cm < 0.05 & alldegs$avg_logfc_cm < -fcCutOff), ]

	# Remove ambiguous DEGs from human specific genes
	hsp_up = hsp_up[!(hsp_up$Gene %in% hup_cdown$Gene),]
	hsp_down = hsp_down[!(hsp_down$Gene %in% hdown_cup$Gene),]

	# Assign decision - human
	alldegs[match(hsp_up$Gene, alldegs$Gene), 'Evolution'] = 'Human_Specific'
	alldegs[match(hsp_up$Gene, alldegs$Gene), 'Regulation'] = 'Human_UP'
	alldegs[match(hsp_down$Gene, alldegs$Gene), 'Evolution'] = 'Human_Specific'
	alldegs[match(hsp_down$Gene, alldegs$Gene), 'Regulation'] = 'Human_DOWN'
	alldegs[match(hup_cdown$Gene, alldegs$Gene), 'Evolution'] = 'Unspecific'
	alldegs[match(hup_cdown$Gene, alldegs$Gene), 'Regulation'] = 'H>M>C'
	alldegs[match(hdown_cup$Gene, alldegs$Gene), 'Evolution'] = 'Unspecific'
	alldegs[match(hdown_cup$Gene, alldegs$Gene), 'Regulation'] = 'H<M<C'

	# Find and label chimp specific genes
	csp_up = alldegs[(alldegs$adj_p_value_hc < 0.05 & alldegs$avg_logfc_hc < -fcCutOff) &
			  (alldegs$adj_p_value_cm < 0.05 & alldegs$avg_logfc_cm > fcCutOff), ]

	csp_down = alldegs[(alldegs$adj_p_value_hc < 0.05 & alldegs$avg_logfc_hc > fcCutOff) &
			    (alldegs$adj_p_value_cm < 0.05 & alldegs$avg_logfc_cm < -fcCutOff), ]

	# Remove ambiguous DEGs from human specific genes
	csp_up = csp_up[!(csp_up$Gene %in% hdown_cup$Gene),]
	csp_down = csp_down[!(csp_down$Gene %in% hup_cdown$Gene),]

	# Assign decision - chimpanzee
	alldegs[match(csp_up$Gene, alldegs$Gene), 'Evolution'] = 'Chimp_Specific'
	alldegs[match(csp_up$Gene, alldegs$Gene), 'Regulation'] = 'Chimp_UP'
	alldegs[match(csp_down$Gene, alldegs$Gene), 'Evolution'] = 'Chimp_Specific'
	alldegs[match(csp_down$Gene, alldegs$Gene), 'Regulation'] = 'Chimp_DOWN'

	# Find and label DEGs between Macaque and HC-LCA
	alldegs$MvsLCA = 'NS'

	m_up = alldegs[(alldegs$adj_p_value_hm < 0.05 & alldegs$avg_logfc_hm < -fcCutOff) &
			 (alldegs$adj_p_value_cm < 0.05 & alldegs$avg_logfc_cm < -fcCutOff), ]

	m_down = alldegs[(alldegs$adj_p_value_hm < 0.05 & alldegs$avg_logfc_hm > fcCutOff) &
			   (alldegs$adj_p_value_cm < 0.05 & alldegs$avg_logfc_cm > fcCutOff), ]

	alldegs[match(m_up$Gene, alldegs$Gene), 'MvsLCA'] = 'Macaque_UP'
	alldegs[match(m_down$Gene, alldegs$Gene), 'MvsLCA'] = 'Macaque_DOWN'

	alldegs$cluster = ctypes[i]
	genelist_list[[i]] = alldegs
}

# Save all results
toplot = do.call(rbind, genelist_list)
toplot_all = toplot
export(toplot_all, 'Excitatory_AllGenes_DEGs.xlsx')

# Keep only significant changes
toplot = toplot[(toplot$Evolution %in% c('Human_Specific', 'Chimp_Specific') | toplot$MvsLCA != 'NS'),]

# Remove ribosomal and mitochondrial genes from significant genes for downstream comparisons
toplot = toplot[!(grepl('^MT-|^RPL|^RPS|^MRPL',toplot$Gene)),]

export(toplot, 'Excitatory_SignGenes_DEGs.xlsx')

##
## Test human vs chimpanzee specificity with randomized background
##

toplot_all = import('Excitatory_AllGenes_DEGs.xlsx')
toplot = import('Excitatory_SignGenes_DEGs.xlsx')
ctypes = names(table(toplot$cluster))
pvalHC = list()
pvalCH = list()
pvalHM = list()
pvalMH = list()
dfs = list()
for(i in 1:length(ctypes)){
	q = toplot_all[toplot_all$cluster == ctypes[i],]

	q$avg_logfc_hc = sample(q$avg_logfc_hc)
	q$avg_logfc_hm = sample(q$avg_logfc_hm)
	q$avg_logfc_cm = sample(q$avg_logfc_cm)

	hs = toplot[toplot$cluster == ctypes[i] & toplot$Evolution %in% c('Human_Specific'),] %>% nrow
	cs = toplot[toplot$cluster == ctypes[i] & toplot$Evolution %in% c('Chimp_Specific'),] %>% nrow
	ms = toplot[toplot$cluster == ctypes[i] & toplot$MvsLCA != 'NS',] %>% nrow

	hsrands = list()
	csrands = list()
	msrands = list()
	hsrats = list()
	hmrats = list()
	for(j in 1:10000){
		q2 = q[sample(1:nrow(q), 500), ]
		hsrands[[j]] = sum((q2$avg_logfc_hc > 0 & q2$avg_logfc_hm > 0) | (q2$avg_logfc_hc < 0 & q2$avg_logfc_hm < 0))
		csrands[[j]] = sum((q2$avg_logfc_hc < 0 & q2$avg_logfc_cm > 0) | (q2$avg_logfc_hc > 0 & q2$avg_logfc_cm < 0))
		msrands[[j]] = sum((q2$avg_logfc_hm < 0 & q2$avg_logfc_cm < 0) | (q2$avg_logfc_hm > 0 & q2$avg_logfc_cm > 0))
		hsrats[[j]] = hsrands[[j]] / csrands[[j]]
		hmrats[[j]] = hsrands[[j]] / msrands[[j]]
	}

	pvalHC[[i]] = 1- (sum(hs/cs > unlist(hsrats)) / length(hsrats))
	pvalCH[[i]] = (sum(hs/cs > unlist(hsrats)) / length(hsrats))
	pvalHM[[i]] = 1- (sum(hs/(ms*6/45) > unlist(hsrats)) / length(hsrats))
	pvalMH[[i]] = (sum(hs/(ms*6/45) > unlist(hmrats)) / length(hmrats))

	dfRAND = data.frame(hsrat = unlist(hsrats), hmrat = unlist(hmrats), CellType = ctypes[i], type = 'RAND')
	dfOBS = data.frame(hsrat = hs/cs, hmrat = hs/(ms*6/45), CellType = ctypes[i], type = 'OBS')
	dfs[[i]] = rbind(dfRAND, dfOBS)
}

randDF = do.call(rbind, dfs)
pvalDF = data.frame(pvalHC = unlist(pvalHC), pvalCH = unlist(pvalCH), pvalHM = unlist(pvalHM), pvalMH = unlist(pvalMH), CellType = ctypes)
pvalDF$pvalHC = p.adjust(pvalDF$pvalHC, method = 'BH')
pvalDF$pvalCH = p.adjust(pvalDF$pvalCH, method = 'BH')
pvalDF$pvalHM = p.adjust(pvalDF$pvalHM, method = 'BH')
pvalDF$pvalMH = p.adjust(pvalDF$pvalMH, method = 'BH')
finalDF = merge(randDF, pvalDF, by = 'CellType')

saveRDS(finalDF, 'HC_HM_RATIOS.RDS')


# H-C
hcPlot = finalDF
hcPlot$is_sign = 'NS'
hcPlot[hcPlot$pvalHC < 0.05, 'is_sign'] = 'H>C'
hcPlot[hcPlot$pvalCH < 0.05, 'is_sign'] = 'H<C'

pdf('HC_Ratio_RNA_EXC.pdf', width = 7, height = 8)
ggplot(hcPlot, aes(y = CellType, x = hsrat)) + xlim(0.5,2) +
geom_boxplot(outlier.shape = NA, color = 'black') +
geom_text(data=hcPlot[hcPlot$type == 'OBS',], aes(label=is_sign, color = is_sign), nudge_x = 0.1, fontface = 'bold', size = 7) +
scale_colour_manual(values=c('blue', 'black')) +
theme_classic() +
theme(text = element_text(size=25, face = 'bold')) +
ylab('') + xlab('') +
rotate_x_text(90) +
geom_point(data = hcPlot[hcPlot$type == 'OBS',], color = 'red') +
theme(strip.text.x = element_text(angle = 90)) +
grids(linetype = "dashed", color = 'grey')
dev.off()


# H-M
hmPlot = finalDF
hmPlot$is_sign = 'NS'
hmPlot[hmPlot$pvalHM < 0.05, 'is_sign'] = 'H>C'
hmPlot[hmPlot$pvalMH < 0.05, 'is_sign'] = 'H<C'

pdf('HM_Ratio_RNA_EXC.pdf', width = 7, height = 8)
ggplot(hmPlot, aes(y = CellType, x = hmrat)) + xlim(0.5,3) +
geom_boxplot(outlier.shape = NA, color = 'black') +
geom_text(data=hmPlot[hmPlot$type == 'OBS',], aes(label=is_sign, color = is_sign), nudge_x = 0.1, fontface = 'bold', size = 7) +
scale_colour_manual(values=c('blue', 'black')) +
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

# Plot total unique species specific genes
toplot_all = import('Excitatory_AllGenes_DEGs.xlsx')
toplot = import('Excitatory_SignGenes_DEGs.xlsx')
htot = length(unique(toplot[toplot$Evolution == 'Human_Specific', 'Gene'])) # Total hspecific
ctot = length(unique(toplot[toplot$Evolution == 'Chimp_Specific', 'Gene'])) # Total cspecific
mtot = length(unique(toplot[toplot$MvsLCA != 'NS', 'Gene']))*6/45 # Total mspecific, normalized
df = data.frame(vals = c(htot,ctot,mtot), vars = c('Human', 'Chimp', 'Macaque'))

pdf('EXC_Unique_DEGs.pdf')
ggbarplot(df, x = 'vars', y = 'vals', fill = 'vars', palette = c('blue','orange','darkgreen')) +
rotate_x_text(90) +
ylab('Total Unique DARs') +
scale_colour_manual(values=c('blue', 'black')) +
theme_classic() +
theme(text = element_text(size=20, face = 'bold')) +
ylab('# of unique genes') + xlab('') +
rotate_x_text(90) +
NoLegend()
dev.off()


# Plot human specific genes for validation
library(gridExtra)
randhsdf = toplot %>% filter(Evolution == 'Human_Specific') %>%
			group_by(cluster) %>% slice_sample(n=1) %>% as.data.frame
reg = randhsdf$Regulation
randhs = randhsdf$Gene
cls = randhsdf$cluster

Idents(allrna) = allrna$Species
ggs = list()
for(i in 1:length(ctypes)){
	ggs[[i]] = VlnPlot(subset(allrna, subset = newannot == cls[i]),
				randhs[i], pt.size = 0, assay = 'RNA') +
				geom_jitter(alpha = 0.1) +
				ggtitle(paste(cls[i], randhs[i], reg[i], sep = '\n'))
}

pdf("EXC_hspecific_deg_per_celltype.pdf", width = 24, height = 18)
grid.arrange(grobs = ggs, nrow = 7, ncol = 10)
dev.off()

# Plot specific genes
library(gridExtra)
gns = c('FOXP2')
ctypes = names(table(toplot$cluster))

for(j in 1:length(gns)){

	Idents(allrna) = factor(allrna$Species, levels = rev(c('human', 'chimp', 'macaque')))
	ggs = list()
	for(i in 1:length(ctypes)){
		ggs[[i]] = VlnPlot(subset(allrna, subset = newannot == ctypes[i]),
					gns[j], pt.size = 0, assay = 'RNA', cols = rev(c('blue', 'orange', 'darkgreen'))) +
					ggtitle(ctypes[i]) + NoLegend() +
					theme(axis.text.x = element_text(size=20),
						axis.text.y = element_text(size=20),
						axis.title = element_text(size=20)) + coord_flip()
	}

	pdf(paste0(gns[j], "_VlnPlot.pdf"), height = 8, width = 20)
	print( grid.arrange(grobs = ggs, nrow = 3) )
	dev.off()
}


# Plot up/down genes per species
library(reshape2)
toplotm = melt(id.vars = c('Gene', 'cluster'), measure.vars = c('MvsLCA','Regulation') ,toplot)
toplotm = toplotm[!(toplotm$value == 'NS'),]

toplotm$Gene = 1
toplotm = aggregate(Gene~value+cluster, toplotm, FUN = sum)
toplotm$value = factor(toplotm$value,
			  levels = c('Human_UP', 'Human_DOWN', 'Chimp_UP',
					'Chimp_DOWN', 'Macaque_UP', 'Macaque_DOWN'))
toplotm$updown = factor(gsub('.*_', '', toplotm$value))


toplotm$Gene = ifelse(toplotm$updown == 'UP', toplotm$Gene, -toplotm$Gene)
toplotm$Species = gsub('_.*', '', toplotm$value)
toplotm$Species = factor(toplotm$Species, levels = rev(c('Human', 'Chimp', 'Macaque')))

toplotm$cluster2 = toplotm$cluster
toplotm$cluster2 = gsub('L4-5_RORB_2-L4-6_RORB_1', 'L4-6_RORB_3', toplotm$cluster2)
toplotm$cluster2 = gsub('_THEMIS', "_\nTHEMIS", toplotm$cluster2)
toplotm$cluster2 = gsub('_FEZF2', "_\nFEZF2", toplotm$cluster2)
toplotm$cluster2 = gsub('_RORB', "_\nRORB", toplotm$cluster2)
toplotm$updown = gsub('UP', "GAIN", toplotm$updown)
toplotm$updown = gsub('DOWN', "LOSS", toplotm$updown)
toplotm$updown = factor(toplotm$updown, levels = c('LOSS', 'GAIN'))


pdf('EXC_DEGs_SpeciesSpecific_UPDOWN_VERTICAL.pdf', width = 10, height = 12)
ggbarplot(toplotm, x = 'Species', y = 'Gene', fill = 'updown', palette = c('red', 'blue')) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20),
		strip.text = element_text(size = 15, face = 'bold'),
		legend.text=element_text(size=20),
		legend.title=element_text(size=20)) +
xlab('') + ylab('# of genes') + rotate_x_text(90) +
facet_wrap(~cluster2, ncol = 2, strip.position = 'right') +
coord_flip()
dev.off()

####
## DISTRIBUTION OF SUBTYPE SPECIFIC FEATURES AFTER EQUALIZING #OF FEATURES
####

toplot = import('Excitatory_SignGenes_DEGs.xlsx')
toplot = toplot %>% rename(Gene = 'Gene')
toplot_up = toplot[toplot$Regulation %in% c('Human_UP', 'Chimp_UP') |
			toplot$MvsLCA == 'Macaque_UP',]
toplot_down = toplot[toplot$Regulation %in% c('Human_DOWN', 'Chimp_DOWN') |
			toplot$MvsLCA == 'Macaque_DOWN',]

ln = length(table(toplot$cluster))

aggL = list()
sn = 500 # Number of features to subset. Results are robust to variations of this number.
for(i in 1:100){

	# Randomly subset open chromatin changes
	rand_freq_up = toplot %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Random', type = 'UP')

	# Randomly subset closed chromatin changes
	rand_freq_down = toplot %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Random', type = 'DOWN')

	colnames(rand_freq_up)[1] = 'FoundIn'
	colnames(rand_freq_down)[1] = 'FoundIn'

	rand_freq_up$Gene_ratios = prop.table(rand_freq_up$Freq)
	rand_freq_down$Gene_ratios = prop.table(rand_freq_down$Freq)

	rand_freq_up$change_ratios = prop.table( (rand_freq_up$Freq) * as.numeric(rand_freq_up$FoundIn) )
	rand_freq_down$change_ratios = prop.table( (rand_freq_down$Freq) * as.numeric(rand_freq_down$FoundIn) )

	# HUMAN
	hs_freq_up = toplot_up %>% filter(Regulation == 'Human_UP') %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Human', type = 'UP')

	hs_freq_down = toplot_down %>% filter(Regulation == 'Human_DOWN') %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Human', type = 'DOWN')

	colnames(hs_freq_up)[1] = 'FoundIn'
	colnames(hs_freq_down)[1] = 'FoundIn'

	hs_freq_up$Gene_ratios = prop.table(hs_freq_up$Freq)
	hs_freq_down$Gene_ratios = prop.table(hs_freq_down$Freq)

	hs_freq_up$change_ratios = prop.table( (hs_freq_up$Freq) * as.numeric(hs_freq_up$FoundIn) )
	hs_freq_down$change_ratios = prop.table( (hs_freq_down$Freq) * as.numeric(hs_freq_down$FoundIn) )

	# CHIMP
	cs_freq_up = toplot_up %>% filter(Regulation == 'Chimp_UP') %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Chimp', type = 'UP')

	cs_freq_down = toplot_down %>% filter(Regulation == 'Chimp_DOWN') %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Chimp', type = 'DOWN')

	colnames(cs_freq_up)[1] = 'FoundIn'
	colnames(cs_freq_down)[1] = 'FoundIn'

	cs_freq_up$Gene_ratios = prop.table(cs_freq_up$Freq)
	cs_freq_down$Gene_ratios = prop.table(cs_freq_down$Freq)

	cs_freq_up$change_ratios = prop.table( (cs_freq_up$Freq) * as.numeric(cs_freq_up$FoundIn) )
	cs_freq_down$change_ratios = prop.table( (cs_freq_down$Freq) * as.numeric(cs_freq_down$FoundIn) )

	# MACAQUE
	ms_freq_up = toplot_up %>% filter(MvsLCA == 'Macaque_UP') %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Macaque', type = 'UP')

	ms_freq_down = toplot_down %>% filter(MvsLCA == 'Macaque_DOWN') %>% sample_n(sn, replace = F) %>% group_by(Gene) %>% group_size %>% factor(levels = 1:ln) %>% table %>% as.data.frame %>% add_column(Species = 'Macaque', type = 'DOWN')

	colnames(ms_freq_up)[1] = 'FoundIn'
	colnames(ms_freq_down)[1] = 'FoundIn'

	ms_freq_up$Gene_ratios = prop.table(ms_freq_up$Freq)
	ms_freq_down$Gene_ratios = prop.table(ms_freq_down$Freq)

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
	comFreqsAgg = aggregate(Gene_ratios ~ presence + Species + type, combFreqs, FUN = sum)
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

pdf('EXC_RNA_SubtypeSpecificity_UP_SPECIES.pdf')
ggboxplot(tmpdf, y = 'Gene_ratios', x = 'Species', outlier.shape = NA, color = 'Species', palette = c('grey', 'blue', 'orange', 'darkgreen'), ylim = c(0.4,1.1)) + ylab('Ratio of subtype specific Genes') + xlab('') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = c(1.02,1.06,1.1)) +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

# Extract p-values
hRats = tmpdf[tmpdf$Species == 'Human', 'Gene_ratios']
cRats = tmpdf[tmpdf$Species == 'Chimp', 'Gene_ratios']
mRats = tmpdf[tmpdf$Species == 'Macaque', 'Gene_ratios']

hcComp = wilcox.test(hRats, cRats)$p.value
hmComp = wilcox.test(hRats, mRats)$p.value
cmComp = wilcox.test(cRats, mRats)$p.value

upcomps = data.frame(pvals = c(hcComp, hmComp, cmComp), comparisons = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))

rio::export(upcomps, 'SubtypeSpecificity_Up_Comps_INH_RNA.xlsx')

# DOWN
aggDF_DOWN = aggDF[aggDF$type == 'DOWN' & aggDF$presence == 'Specific',]
comps = list(c('Human', 'Chimp'), c('Human', 'Macaque'), c('Chimp', 'Macaque'))
tmpdf = aggDF_DOWN[aggDF_DOWN$Species %in% c('Human', 'Chimp', 'Macaque', 'Random'),]

pdf('EXC_RNA_SubtypeSpecificity_DOWN_SPECIES.pdf')
ggboxplot(tmpdf, y = 'Gene_ratios', x = 'Species', outlier.shape = NA, color = 'Species', palette = c('grey', 'blue', 'orange', 'darkgreen'), ylim = c(0.4,1.1)) + ylab('Ratio of subtype specific Genes') + xlab('') +
stat_compare_means(comparisons = comps, method = 'wilcox.test', label.y = c(1.02,1.06,1.1)) +
theme(text = element_text(size=20)) +
theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
rotate_x_text(90) +
NoLegend()
dev.off()

# Extract p-values
hRats = tmpdf[tmpdf$Species == 'Human', 'Gene_ratios']
cRats = tmpdf[tmpdf$Species == 'Chimp', 'Gene_ratios']
mRats = tmpdf[tmpdf$Species == 'Macaque', 'Gene_ratios']

hcComp = wilcox.test(hRats, cRats)$p.value
hmComp = wilcox.test(hRats, mRats)$p.value
cmComp = wilcox.test(cRats, mRats)$p.value

upcomps = data.frame(pvals = c(hcComp, hmComp, cmComp), comparisons = c('Human-Chimp', 'Human-Macaque', 'Chimp-Macaque'))

rio::export(upcomps, 'SubtypeSpecificity_Down_Comps_INH_RNA.xlsx')



