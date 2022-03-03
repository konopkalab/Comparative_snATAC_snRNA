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
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

####
## LOAD DATASETS
####

# All DARs
dars = fread('COMBINED_Signpeaks_DARs.txt', header = T)
dars = as.data.frame(dars)

darsEXC_ALL = dars[dars$major == 'Excitatory', 'peak'] %>% unique() %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsEXC_HS = dars[dars$major == 'Excitatory' & dars$Regulation == 'Human_UP','peak'] %>% unique() %>% gsub(':|-', '_', .) %>% 			strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsEXC_CS = dars[dars$major == 'Excitatory' & dars$Regulation == 'Chimp_UP','peak'] %>% unique() %>% gsub(':|-', '_', .) %>% 			strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_ALL = dars[dars$major == 'Inhibitory','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_HS = dars[dars$major == 'Inhibitory' & dars$Regulation == 'Human_UP','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_CS = dars[dars$major == 'Inhibitory' & dars$Regulation == 'Chimp_UP','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

# BOULTING CREs
bou = import_list('BOULTING_2021/41593_2020_786_MOESM5_ESM.xlsx')

bou15mInc = bou[[5]][-1,8:10] %>% dplyr::rename(chr = '...8', start = '...9', end = '...10') %>% mutate(start = as.numeric(start), end = as.numeric(end))
bou15mDec = bou[[6]][-1,8:10] %>% dplyr::rename(chr = '...8', start = '...9', end = '...10') %>% mutate(start = as.numeric(start), end = as.numeric(end))

bou2hInc = bou[[7]][-1,8:10] %>% dplyr::rename(chr = '...8', start = '...9', end = '...10') %>% mutate(start = as.numeric(start), end = as.numeric(end))
bou2hDec = bou[[8]][-1,8:10] %>% dplyr::rename(chr = '...8', start = '...9', end = '...10') %>% mutate(start = as.numeric(start), end = as.numeric(end))

# Keep only 15m or only 2h
bou15mInc$ID = paste0(bou15mInc[,1], bou15mInc[,2], bou15mInc[,3])
bou15mDec$ID = paste0(bou15mDec[,1], bou15mDec[,2], bou15mDec[,3])
bou2hInc$ID = paste0(bou2hInc[,1], bou2hInc[,2], bou2hInc[,3])
bou2hDec$ID = paste0(bou2hDec[,1], bou2hDec[,2], bou2hDec[,3])

bou15mIncOnly = bou15mInc[!(bou15mInc$ID %in% bou2hInc$ID), -4]
bou15mDecOnly = bou15mDec[!(bou15mDec$ID %in% bou2hInc$ID), -4]
bou2hIncOnly = bou2hInc[!(bou2hInc$ID %in% bou15mInc$ID), -4]
bou2hDecOnly = bou2hDec[!(bou2hDec$ID %in% bou15mDec$ID), -4]

# LiftOver to hg38
hg19hg38 = import.chain('hg19ToHg38.over.chain')

rownames(bou15mInc) = 1:nrow(bou15mInc)
tmp = makeGRangesFromDataFrame(bou15mInc)
tmp = liftOver(tmp, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))
bou15mInc = tmp[,1:3]

rownames(bou15mDec) = 1:nrow(bou15mDec)
tmp = makeGRangesFromDataFrame(bou15mDec)
tmp = liftOver(tmp, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))
bou15mDec = tmp[,1:3]

rownames(bou2hInc) = 1:nrow(bou2hInc)
tmp = makeGRangesFromDataFrame(bou2hInc)
tmp = liftOver(tmp, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))
bou2hInc = tmp[,1:3]

rownames(bou2hDec) = 1:nrow(bou2hDec)
tmp = makeGRangesFromDataFrame(bou2hDec)
tmp = liftOver(tmp, chain = hg19hg38) %>% unlist %>% as.data.frame() %>% mutate(seqnames = as.character(seqnames))
bou2hDec = tmp[,1:3]

bou15mInc = bedr.sort.region(bou15mInc, verbose = F)
bou15mDec = bedr.sort.region(bou15mDec, verbose = F)
bou2hInc = bedr.sort.region(bou2hInc, verbose = F)
bou2hDec = bedr.sort.region(bou2hDec, verbose = F)

####
## OVERLAP WITH ALL BOULTING CREs
####

## 15 min vs 2 HOURS ##
colnames(testINH) = colnames(bou2hInc)
colnames(darsEXC_HS) = colnames(bou2hInc)
colnames(darsINH_HS) = colnames(bou2hInc)
colnames(darsEXC_CS) = colnames(bou2hInc)
colnames(darsINH_CS) = colnames(bou2hInc)

# EXCITATORY - INCREASING - HUMAN
ovEXC_HS_15m = bedr(input = list(a = darsEXC_HS, b = bou15mInc), method = "intersect", params = "-loj", verbose = F)
ovEXC_HS_2h = bedr(input = list(a = darsEXC_HS, b = bou2hInc), method = "intersect", params = "-loj", verbose = F)

ovEXC_HS_15m_sum = sum(ovEXC_HS_15m[,5] != '-1')
ovEXC_HS_2h_sum = sum(ovEXC_HS_2h[,5] != '-1')

excIncHS = prop.test(c(ovEXC_HS_15m_sum, ovEXC_HS_2h_sum), c(nrow(ovEXC_HS_15m), nrow(ovEXC_HS_2h)))

# EXCITATORY - INCREASING - CHIMP
ovEXC_CS_15m = bedr(input = list(a = darsEXC_CS, b = bou15mInc), method = "intersect", params = "-loj", verbose = F)
ovEXC_CS_2h = bedr(input = list(a = darsEXC_CS, b = bou2hInc), method = "intersect", params = "-loj", verbose = F)

ovEXC_CS_15m_sum = sum(ovEXC_CS_15m[,5] != '-1')
ovEXC_CS_2h_sum = sum(ovEXC_CS_2h[,5] != '-1')

excIncCS = prop.test(c(ovEXC_CS_15m_sum, ovEXC_CS_2h_sum), c(nrow(ovEXC_CS_15m), nrow(ovEXC_CS_2h)))


# INHIBITORY - INCREASING - HUMAN
ovINH_HS_15m = bedr(input = list(a = darsINH_HS, b = bou15mInc), method = "intersect", params = "-loj", verbose = F)
ovINH_HS_2h = bedr(input = list(a = darsINH_HS, b = bou2hInc), method = "intersect", params = "-loj", verbose = F)

ovINH_HS_15m_sum = sum(ovINH_HS_15m[,5] != '-1')
ovINH_HS_2h_sum = sum(ovINH_HS_2h[,5] != '-1')

inhIncHS = prop.test(c(ovINH_HS_15m_sum, ovINH_HS_2h_sum), c(nrow(ovINH_HS_15m), nrow(ovINH_HS_2h)))


# INHIBITORY - INCREASING - CHIMP
ovINH_CS_15m = bedr(input = list(a = darsINH_CS, b = bou15mInc), method = "intersect", params = "-loj", verbose = F)
ovINH_CS_2h = bedr(input = list(a = darsINH_CS, b = bou2hInc), method = "intersect", params = "-loj", verbose = F)

ovINH_CS_15m_sum = sum(ovINH_CS_15m[,5] != '-1')
ovINH_CS_2h_sum = sum(ovINH_CS_2h[,5] != '-1')

inhIncCS = prop.test(c(ovINH_CS_15m_sum, ovINH_CS_2h_sum), c(nrow(ovINH_CS_15m), nrow(ovINH_CS_2h)))


# PLOT INCREASING
library('RColorBrewer')

pvals = c(excIncHS$p.value, inhIncHS$p.value, excIncCS$p.value, inhIncCS$p.value)
cellType = c(rep('EXC', 4), rep('INH', 4))
timeInt = c('15m', '2h', '15m', '2h', '15m', '2h', '15m', '2h')
excVals = c(ovEXC_HS_15m_sum, ovEXC_HS_2h_sum, ovEXC_CS_15m_sum, ovEXC_CS_2h_sum)
inhVals = c(ovINH_HS_15m_sum, ovINH_HS_2h_sum, ovINH_CS_15m_sum, ovINH_CS_2h_sum)

toplot = c(excVals, inhVals) %>%
		as.data.frame %>%
		dplyr::rename(Overlap = ".") %>%
		add_column(cellType = cellType) %>%
		add_column(timeInt = timeInt) %>%
		add_column(Species = c(rep('HUMAN', 2), rep('CHIMP', 2), rep('HUMAN', 2), rep('CHIMP', 2))) %>%
		add_column(FDR = rep(p.adjust(pvals, method = 'BH'), 2))


toplotM = melt(toplot, id.vars = c('FDR', 'cellType', 'Species', 'timeInt'))
toplotM$cellType = factor(toplotM$cellType, levels = c('EXC', 'INH'))
toplotM$Species = factor(toplotM$Species, levels = c('HUMAN', 'CHIMP'))


# Plot
plt = ggbarplot(toplotM, x = 'timeInt', y = 'value', color = 'Species',
	fill = 'Species', palette = c('blue', 'orange'), ylim = c(0,175)) +
	 facet_wrap(~Species + cellType, nrow = 1) +
	 ylab('Number Of Species-Specific\n And Activity Regulated CREs') + xlab('Time after stimulation') +
	 rotate_x_text(90) +
	 theme(text=element_text(size=20, face = 'bold')) + NoLegend()

dat_text <- data.frame(label = c("*", "ns"),  Species = factor(c('HUMAN', 'CHIMP'), levels = c('HUMAN', 'CHIMP')))

pdf('Boulting_Compare_15m_to_2h.pdf', width = 8)
plt + geom_text(data = dat_text, size = 6, fontface = 'bold', mapping = aes(x = 1, y = 170, label = label), nudge_x = 0.4)
dev.off()






