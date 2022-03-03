rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggpubr)
library(reshape2)
library(data.table)
library(Signac)
library(Seurat)
library(dplyr)
library(tidyverse)
library(rio)
library(GenomicRanges)
library(liftOver)
library('bedr')
source("utility_functions.R")

####
## MERGE HAR DATASETS
####

# Read HARs
doan = rio::import('Doan_2016_Compendium.xlsx')
gittel = rio::import('Gittelman_2015_hg19.xlsx')

doan = doan[, c(1:3)]
gittel = gittel[, c(1:3)]
colnames(gittel) = colnames(doan)

doan = makeGRangesFromDataFrame(doan)
gittel = makeGRangesFromDataFrame(gittel)

# LiftOver HARs to hg38
hg19hg38 = import.chain('hg19ToHg38.over.chain')

doan_hg38 = liftOver(doan, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))
gittel_hg38 = liftOver(gittel, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))

# Merge all HARs
doan_hg38 = bedr.sort.region(doan_hg38, verbose = F)
gittel_hg38 = bedr.sort.region(gittel_hg38, verbose = F)

allHars = rbind(doan_hg38, gittel_hg38)
allHars = bedr.merge.region(allHars) # 3267 -> 3176
rio::export(allHars %>% dplyr::rename(width = 'names'), 'Doan_Combined_HARs_hg38.xlsx')

# Read HARs
allHars = rio::import('Doan_Combined_HARs_hg38.xlsx')
allHars = allHars[, c(1:3)]

####
## OVERLAP PEAKS AND HARs
####

# DARs (Differential CREs)
dars = fread('COMBINED_Signpeaks_DARs.txt', header = T)
dars = as.data.frame(dars)
dars[dars$cluster == 'OPC', 'major'] = 'OPC'

# All CREs
allpks = fread('COMBINED_Allpeaks_DARs.txt', header = T)
allpks = as.data.frame(allpks)
allpks[allpks$cluster == 'OPC', 'major'] = 'OPC'

# Overlap with all CREs
bcg_bed = allpks$peak %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3)))
bcg_bed = cbind(bcg_bed, peak = allpks$peak %>% gsub(':|-', '_', .))
bcg_bed = bcg_bed[!duplicated(bcg_bed$peak),]
bcg_bed = bedr.sort.region(bcg_bed, verbose = F)
bcg_bed = bedr(input = list(a = bcg_bed, b = allHars), method = "intersect", params = c("-loj"), verbose = F)

bcg_bed$start.b = as.numeric(bcg_bed$start.b)
bcg_bed$end.b = as.numeric(bcg_bed$end.b)
bcg_bed$ov = pmin(bcg_bed$V3, bcg_bed$end.b) - pmax(bcg_bed$V2, bcg_bed$start.b)
bcg_bed$ov = ifelse(bcg_bed$ov < 0, 0, bcg_bed$ov)

# Assign species specific CREs
hsdars = dars[dars$Evolution == 'Human_Specific', 'peak'] %>% gsub(':|-', '_', .) %>% unique
csdars = dars[dars$Evolution == 'Chimp_Specific', 'peak'] %>% gsub(':|-', '_', .) %>% unique

bcg_bed$is_hs = ifelse(bcg_bed$peak %in% hsdars, 'HS', 'NS')
bcg_bed$is_cs = ifelse(bcg_bed$peak %in% csdars, 'CS', 'NS')
bcg_bed$is_hs = factor(bcg_bed$is_hs, levels = c('NS', 'HS'))
bcg_bed$is_cs = factor(bcg_bed$is_cs, levels = c('NS', 'CS'))

# Length of CRE
bcg_bed$length = bcg_bed$V3 - bcg_bed$V2

# Save HAR-DAR overlap
darmeta = dars[, c('peak', 'Evolution', 'Regulation', 'MvsLCA', 'cluster')]
darmeta$peak = gsub(':|-', '_', darmeta$peak)
tosave = merge(bcg_bed, darmeta, by.y = 'peak')
tosave = tosave[tosave$ov > 0,]
colnames(tosave) = c('CRE_peak', 'CRE_chr', 'CRE_start', 'CRE_end', 'CRE_chr', 'CRE_start',
			'HAR_end', 'Overlap', 'is_hs', 'is_cs', 'length', 'Evolution', 'Regulation', 'MvsHC', 'CellType')
tosave$HARCoord = paste0(tosave$HAR_chr, ':', tosave$HAR_start, '-', tosave$HAR_end)
rio::export(tosave, 'HAR_DAR_Ov.xlsx')

# Plot number of HARs overlapping human or chimp
hardf = table(tosave$Evolution) %>% as.data.frame
hardf = hardf[-nrow(hardf),]
pdf('HAR_HumanChimp_Number.pdf')
ggbarplot(hardf, x = 'Var1', y = 'Freq', fill = 'Var1', color = 'Var1') + NoLegend() +
ylab('Number of HAR-DAR Overlap') + xlab('')
dev.off()

# Plot distribution of non-zeros
pdf('Bcg_HAR_overlap_amount_dist.pdf')
hist(bcg_bed[bcg_bed$ov > 0, 'ov'], 100)
dev.off()

####
## ENRICHMENT BY REGRESSION (MAJOR CELL TYPES)
####

bcg_bed = bcg_bed[!duplicated(bcg_bed$peak),]
bcg_bed$length = log2(bcg_bed$length)

# Test enrichment on all CREs. Fit length as covariate. Collect results in a data frame
hres <- glm(ov ~ length + is_hs, data = bcg_bed)
cres <- glm(ov ~ length + is_cs, data = bcg_bed)

human_pval = coef(summary(hres))[3,4]
chimp_pval = coef(summary(cres))[3,4]

OR_human = coef(summary(hres))[3,1]
OR_chimp = coef(summary(cres))[3,1]

all_presdf = data.frame(human_pval = human_pval, chimp_pval = chimp_pval, Coef_Human = OR_human, Coef_Chimp = OR_chimp, CellType = 'All_Cells')

require(pscl)
require(MASS)
require(boot)

# Run enrichment per major cell type
ctypes = names(table(allpks$major))
hresCt = list()
cresCt = list()
for(i in 1:length(ctypes)){

	tested = allpks[allpks$major == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)
	subOv = bcg_bed[bcg_bed$peak %in% tested,]

	hsubdars = dars[dars$Evolution == 'Human_Specific' & dars$major == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)
	csubdars = dars[dars$Evolution == 'Chimp_Specific' & dars$major == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)

	subOv$is_hs = ifelse(subOv$peak %in% hsubdars, 'HS', 'NS')
	subOv$is_cs = ifelse(subOv$peak %in% csubdars, 'CS', 'NS')
	subOv$is_hs = factor(subOv$is_hs, levels = c('NS', 'HS'))
	subOv$is_cs = factor(subOv$is_cs, levels = c('NS', 'CS'))

	hresCt[[i]] <- glm(ov ~ length + is_hs, data = subOv)
	cresCt[[i]] <- glm(ov ~ length + is_cs, data = subOv)

	print(i)
}

# Collect results in a data frame
human_pval = sapply(hresCt, function(x){coef(summary(x))[3,4]})
chimp_pval = sapply(cresCt, function(x){coef(summary(x))[3,4]})

OR_human = sapply(hresCt, function(x){coef(summary(x))[3,1]})
OR_chimp = sapply(cresCt, function(x){coef(summary(x))[3,1]})

presdf = data.frame(human_pval = human_pval, chimp_pval = chimp_pval, Coef_Human = OR_human, Coef_Chimp = OR_chimp, CellType = ctypes)
presdf = rbind(all_presdf, presdf)

# Adjust for p-value using FDR
presdf$log10FDR_Human = -log10(p.adjust(presdf$human_pval, method = 'BH'))
presdf$log10FDR_Chimp = -log10(p.adjust(presdf$chimp_pval, method = 'BH'))
presdf[is.na(presdf$log10FDR_Human) | is.nan(presdf$log10FDR_Human), 'log10FDR_Human'] = 0
presdf[is.na(presdf$log10FDR_Chimp) | is.nan(presdf$log10FDR_Chimp), 'log10FDR_Chimp'] = 0

# Plot the results
Coefmelt = melt(measure.vars = c('Coef_Human', 'Coef_Chimp'), variable.name = 'Var1', value.name = 'Coef', presdf[,c(3,4,5)])
fdrmelt = melt(measure.vars = c('log10FDR_Human', 'log10FDR_Chimp'), variable.name = 'Var1', value.name = 'log10FDR', presdf[,c(6,7,5)])
toplot = cbind(Coefmelt, fdrmelt)

toplot$Species = gsub('.*_', '', toplot$Var1)
colnames(toplot) = make.unique(colnames(toplot))
toplot$condition = ifelse(toplot$Coef > 0, 'Enrichment', 'Depletion')

library(ggrepel)
pdf('HAR_DAR_Enrichment_GLM_LOG2LENGTH.pdf', width = 10, height = 5)
ggscatter(toplot, x = 'log10FDR', y = 'CellType', color = 'condition', size = 'log10FDR', palette = c('blue', 'red')) +
geom_label_repel(data = toplot[toplot$log10FDR > 1,], aes(label = CellType), nudge_y = 0.02, nudge_x = -0.05) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
facet_wrap(~Species)
dev.off()


####
## ENRICHMENT BY REGRESSION (SUBTYPE)
####

# Run enrichment per subtype
ctypes = names(table(allpks$cluster))
hresCt = list()
cresCt = list()
for(i in 1:length(ctypes)){

	tested = allpks[allpks$cluster == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)
	subOv = bcg_bed[bcg_bed$peak %in% tested,]

	hsubdars = dars[dars$Evolution == 'Human_Specific' & dars$cluster == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)
	csubdars = dars[dars$Evolution == 'Chimp_Specific' & dars$cluster == ctypes[i], 'peak'] %>% unique %>% gsub(':|-', '_', .)

	subOv$is_hs = ifelse(subOv$peak %in% hsubdars, 'HS', 'NS')
	subOv$is_cs = ifelse(subOv$peak %in% csubdars, 'CS', 'NS')
	subOv$is_hs = factor(subOv$is_hs, levels = c('NS', 'HS'))
	subOv$is_cs = factor(subOv$is_cs, levels = c('NS', 'CS'))

	hresCt[[i]] <- glm(ov ~ length + is_hs, data = subOv)
	cresCt[[i]] <- glm(ov ~ length + is_cs, data = subOv)

	print(i)
}

# Collect results in a data frame
human_pval = sapply(hresCt, function(x){coef(summary(x))[3,4]})
chimp_pval = sapply(cresCt, function(x){coef(summary(x))[3,4]})

OR_human = sapply(hresCt, function(x){coef(summary(x))[3,1]})
OR_chimp = sapply(cresCt, function(x){coef(summary(x))[3,1]})

presdf = data.frame(human_pval = human_pval, chimp_pval = chimp_pval, Coef_Human = OR_human, Coef_Chimp = OR_chimp, CellType = ctypes)
presdf = rbind(all_presdf, presdf)

# Keep only neurons for subtypes
presdf = presdf[grepl('^L|VIP|SST|PVALB|LAMP5|Upper', presdf$CellType), ]


# Adjust for p-value using FDR
presdf$log10FDR_Human = -log10(p.adjust(presdf$human_pval, method = 'BH'))
presdf$log10FDR_Chimp = -log10(p.adjust(presdf$chimp_pval, method = 'BH'))
presdf[is.na(presdf$log10FDR_Human) | is.nan(presdf$log10FDR_Human), 'log10FDR_Human'] = 0
presdf[is.na(presdf$log10FDR_Chimp) | is.nan(presdf$log10FDR_Chimp), 'log10FDR_Chimp'] = 0

# Plot the results
Coefmelt = melt(measure.vars = c('Coef_Human', 'Coef_Chimp'), variable.name = 'Var1', value.name = 'Coef', presdf[,c(3,4,5)])
fdrmelt = melt(measure.vars = c('log10FDR_Human', 'log10FDR_Chimp'), variable.name = 'Var1', value.name = 'log10FDR', presdf[,c(6,7,5)])
toplot = cbind(Coefmelt, fdrmelt)

toplot$Species = gsub('.*_', '', toplot$Var1)
colnames(toplot) = make.unique(colnames(toplot))
toplot$condition = ifelse(toplot$Coef > 0, 'Enrichment', 'Depletion')

library(ggrepel)
pdf('HAR_DAR_Enrichment_GLM_SUBTYPE.pdf', width = 10, height = 10)
ggscatter(toplot, x = 'log10FDR', y = 'CellType', color = 'condition', size = 'log10FDR', palette = c('blue', 'red')) +
geom_label_repel(data = toplot[toplot$log10FDR > 1.3,], aes(label = CellType), nudge_y = 0.02, nudge_x = -0.05) +
geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
facet_wrap(~Species)
dev.off()

# All peaks overlap with HARs
hbcg = read.table('merged_lifted_human.bed')
hbcg = cbind(hbcg, peak = paste0(hbcg$V1, '_', hbcg$V2, '_', hbcg$V3))

hbcg = bedr.sort.region(hbcg, verbose = F)
hbcgOv = bedr(input = list(a = hbcg, b = allHars), method = "intersect", params = "-loj", verbose = F)
hbcgOv = hbcgOv[!duplicated(hbcgOv$peak),]

# Pie-charts
library(scales)

# Find percentag of HARs in all peaks
harRat = sum(hbcgOv$start.b > -1) / nrow(allHars)
darRat = sum(bcg_bed$start.b > -1 & bcg_bed$is_hs == 'HS') / nrow(allHars)

df1 = data.frame(vals = c(harRat, 1 - harRat), vars = c('HAR+', 'HAR-'))
df2 = data.frame(vals = c(darRat, 1 - darRat), vars = c('Human_DAR+', 'Human_DAR-'))

df1$vals = df1$vals * 100
df1 <- df1 %>% 
  arrange(desc(vars)) %>%
  mutate(prop = vals / sum(df1$vals) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

df2$vals = df2$vals * 100
df2 <- df2 %>% 
  arrange(desc(vars)) %>%
  mutate(prop = vals / sum(df2$vals) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )


p1 = ggplot(df1, aes(x="", y=vals, fill=vars))+
	geom_bar(width = 1, stat = "identity") +
	coord_polar("y", start=0) +
	scale_fill_brewer("Blues") + blank_theme +
	  theme(axis.text.x=element_blank()) +
	  geom_text(aes(y = ypos, label = percent(vals/100)), size=4) +
	  ggtitle('Among All Peaks')

p2 = ggplot(df2, aes(x="", y=vals, fill=vars))+
	geom_bar(width = 1, stat = "identity") +
	coord_polar("y", start=0) +
	scale_fill_brewer("Blues") + blank_theme +
	  theme(axis.text.x=element_blank())+
	  geom_text(aes(y = ypos, label = percent(vals/100)), size=4) +
	  ggtitle('Among Human Peaks')

pdf('HAR_AllPeaks_Overlap.pdf', height = 5, width = 8)
p1+p2
dev.off()




