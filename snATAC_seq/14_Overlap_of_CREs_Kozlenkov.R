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

## ALL DARs
dars = fread('COMBINED_Signpeaks_DARs.txt', header = T)
dars = as.data.frame(dars)

darsEXC_ALL = dars[dars$major == 'Excitatory', 'peak'] %>% unique() %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsEXC_HS = dars[dars$major == 'Excitatory' & dars$Evolution == 'Human_Specific','peak'] %>% unique() %>% gsub(':|-', '_', .) %>% 			strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsEXC_CS = dars[dars$major == 'Excitatory' & dars$Evolution == 'Chimp_Specific','peak'] %>% unique() %>% gsub(':|-', '_', .) %>% 			strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsEXC_MHC = dars[dars$major == 'Excitatory' & dars$MvsLCA != 'NS','peak'] %>% unique() %>% gsub(':|-', '_', .) %>% 			strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_ALL = dars[dars$major == 'Inhibitory','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_HS = dars[dars$major == 'Inhibitory' & dars$Evolution == 'Human_Specific','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_CS = dars[dars$major == 'Inhibitory' & dars$Evolution == 'Chimp_Specific','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

darsINH_MHC = dars[dars$major == 'Inhibitory' & dars$MvsLCA != 'NS','peak'] %>% unique() %>% gsub(':|-', '_', .) %>%
		strsplit(., '_') %>% do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

# ALL TESTED
alltest = fread('COMBINED_Allpeaks_DARs.txt', header = T)
alltest = as.data.frame(alltest)

testEXC = alltest[alltest$major == 'Excitatory', 'peak'] %>% unique() %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

testINH = alltest[alltest$major == 'Inhibitory', 'peak'] %>% unique() %>% gsub(':|-', '_', .) %>% strsplit(., '_') %>%
		do.call(rbind, .) %>% as.data.frame() %>%
		mutate(V1 = as.character(V1), V2 = as.numeric(as.character(V2)), V3 = as.numeric(as.character(V3))) %>% bedr.sort.region(., verbose = F)

# ALL PEAKS
allpks = fread('merged_lifted_human_sorted.bed')
allpks = as.data.frame(allpks)
allpks = bedr.sort.region(allpks, verbose = F)

# KOZLENKOV CREs
koz = import_list('KOZLENKOV_2020/pnas.2011884117.sd06.xlsx')
kozGLU = koz[[1]][,1:3]
kozGABA = koz[[2]][,1:3]
kozGLU = bedr.sort.region(kozGLU, verbose = F)
kozGABA = bedr.sort.region(kozGABA, verbose = F)

## KOZLENKOV DA-CREs
kozDA = import_list('KOZLENKOV_2020/pnas.2011884117.sd11.xlsx')

# Excitatory
kozGLU_HUP = kozDA[[1]][,1:3] %>% add_column(Evolution = 'HumanUP')
kozGLU_HDOWN = kozDA[[7]][,1:3] %>% add_column(Evolution = 'HumanDOWN')
kozGLU_HUMAN = rbind(kozGLU_HUP, kozGLU_HDOWN)
kozGLU_HUMAN = kozGLU_HUMAN[!is.na(kozGLU_HUMAN$chr),]
kozGLU_HUMAN = bedr.sort.region(kozGLU_HUMAN, verbose = F)

kozGLU_CUP = kozDA[[2]][,1:3] %>% add_column(Evolution = 'ChimpUP')
kozGLU_CDOWN = kozDA[[8]][,1:3] %>% add_column(Evolution = 'ChimpDOWN')
kozGLU_CHIMP = rbind(kozGLU_CUP, kozGLU_CDOWN)
kozGLU_CHIMP = kozGLU_CHIMP[!is.na(kozGLU_CHIMP$chr),]
kozGLU_CHIMP = bedr.sort.region(kozGLU_CHIMP, verbose = F)

kozGLU_MUP = kozDA[[3]][,1:3] %>% add_column(Evolution = 'MacaqueUP')
kozGLU_MDOWN = kozDA[[9]][,1:3] %>% add_column(Evolution = 'MacaqueDOWN')
kozGLU_MACAQUE = rbind(kozGLU_MUP, kozGLU_MDOWN)
kozGLU_MACAQUE = kozGLU_MACAQUE[!is.na(kozGLU_MACAQUE$chr),]
kozGLU_MACAQUE = bedr.sort.region(kozGLU_MACAQUE, verbose = F)

# Inhibitory
kozGABA_HUP = kozDA[[4]][,1:3] %>% add_column(Evolution = 'HumanUP')
kozGABA_HDOWN = kozDA[[10]][,1:3] %>% add_column(Evolution = 'HumanDOWN')
kozGABA_HUMAN = rbind(kozGABA_HUP, kozGABA_HDOWN)
kozGABA_HUMAN = kozGABA_HUMAN[!is.na(kozGABA_HUMAN$chr),]
kozGABA_HUMAN = bedr.sort.region(kozGABA_HUMAN, verbose = F)

kozGABA_CUP = kozDA[[5]][,1:3] %>% add_column(Evolution = 'ChimpUP')
kozGABA_CDOWN = kozDA[[11]][,1:3] %>% add_column(Evolution = 'ChimpDOWN')
kozGABA_CHIMP = rbind(kozGABA_CUP, kozGABA_CDOWN)
kozGABA_CHIMP = kozGABA_CHIMP[!is.na(kozGABA_CHIMP$chr),]
kozGABA_CHIMP = bedr.sort.region(kozGABA_CHIMP, verbose = F)

kozGABA_MUP = kozDA[[6]][,1:3] %>% add_column(Evolution = 'MacaqueUP')
kozGABA_MDOWN = kozDA[[12]][,1:3] %>% add_column(Evolution = 'MacaqueDOWN')
kozGABA_MACAQUE = rbind(kozGABA_MUP, kozGABA_MDOWN)
kozGABA_MACAQUE = kozGABA_MACAQUE[!is.na(kozGABA_MACAQUE$chr),]
kozGABA_MACAQUE = bedr.sort.region(kozGABA_MACAQUE[,1:3], verbose = F)

####
## OVERLAP WITH ALL KOZLENKOV CREs ##
####
library('RColorBrewer')

# All CREs
colnames(testEXC) = colnames(kozGLU)
merg = bedr.merge.region(rbind(kozGLU, testEXC))
bedPie(bedL = list(merg, kozGLU, testEXC), fn = 'EXC_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG = 364000)

colnames(testINH) =colnames(kozGABA)
merg = bedr.merge.region(rbind(kozGABA, testINH))
bedPie(list(merg, kozGABA, testINH), 'INH_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG = 364000)

# HS-HS CREs GLU
fisherBCG = mean(nrow(testEXC), nrow(kozGLU))

kozGLU_HUMAN = kozGLU_HUMAN[,1:3]
colnames(darsEXC_HS) = colnames(kozGLU_HUMAN)
merg = bedr.merge.region(rbind(kozGLU_HUMAN, darsEXC_HS))
bedPie(bedL = list(merg, kozGLU_HUMAN, darsEXC_HS), 'HS_EXC_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG = fisherBCG)

# HS-EXC, CS-KOZ GLU
kozGLU_CHIMP = kozGLU_CHIMP[,1:3]
colnames(darsEXC_HS) = colnames(kozGLU_CHIMP)
merg = bedr.merge.region(rbind(kozGLU_CHIMP, darsEXC_HS))
bedPie(list(merg, kozGLU_CHIMP, darsEXC_HS), 'HSCS_EXC_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG = fisherBCG)

# HS-HS CREs GABA
fisherBCG = round(mean(c(nrow(testINH), nrow(kozGABA))))

kozGABA_HUMAN = kozGABA_HUMAN[,1:3]
colnames(darsINH_HS) = colnames(kozGABA_HUMAN)
merg = bedr.merge.region(rbind(kozGABA_HUMAN, darsINH_HS))
bedPie(list(merg, kozGABA_HUMAN, darsINH_HS), 'HS_INH_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG)

# HS-EXC, CS-KOZ GABA
kozGABA_CHIMP = kozGABA_CHIMP[,1:3]
colnames(darsINH_HS) = colnames(kozGABA_CHIMP)
merg = bedr.merge.region(rbind(kozGABA_CHIMP, darsINH_HS))
bedPie(list(merg, kozGABA_CHIMP, darsINH_HS), 'HSCS_INH_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG)


# CS-CS CREs
kozGLU_CHIMP = kozGLU_CHIMP[,1:3]
colnames(darsEXC_CS) = colnames(kozGLU_CHIMP)
merg = bedr.merge.region(rbind(kozGLU_CHIMP, darsEXC_CS))
bedPie(list(merg, kozGLU_CHIMP, darsEXC_CS), 'CS_EXC_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG)

kozGABA_CHIMP = kozGABA_CHIMP[,1:3]
colnames(darsINH_CS) = colnames(kozGABA_CHIMP)
merg = bedr.merge.region(rbind(kozGABA_CHIMP, darsINH_CS))
bedPie(list(merg, kozGABA_CHIMP, darsINH_CS), 'CS_INH_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG)


# MvsHC CREs
kozGLU_MACAQUE = kozGLU_MACAQUE[,1:3]
colnames(darsEXC_MHC) = colnames(kozGLU_MACAQUE)
merg = bedr.merge.region(rbind(kozGLU_MACAQUE, darsEXC_MHC))
bedPie(list(merg, kozGLU_MACAQUE, darsEXC_MHC), 'MvsHC_EXC_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG)

kozGABA_MACAQUE = kozGABA_MACAQUE[,1:3]
colnames(darsINH_MHC) = colnames(kozGABA_MACAQUE)
merg = bedr.merge.region(rbind(kozGABA_MACAQUE, darsINH_MHC))
bedPie(list(merg, kozGABA_MACAQUE, darsINH_MHC), 'MvsHC_INH_BA23_KOZLENKOV_ALL_CREs', pienames = c('ONLY_BULK', 'ONLY_BA23', 'BOTH'), fisherBCG)

####
## EPIGENOME VS TRANSCRIPTOME
####

# Transcriptome
kozDE = import_list('KOZLENKOV_2020/pnas.2011884117.sd12.xlsx')

deUP = kozDE[[1]]
deDOWN = kozDE[[2]]

deUPCount = sapply(deUP, function(x){length(na.omit(x))})
deDOWNCount = sapply(deDOWN, function(x){length(na.omit(x))})

huU = deUPCount[grepl('RNA.Hu', names(deUPCount))] %>% sum
huD = deDOWNCount[grepl('RNA.Hu', names(deDOWNCount))] %>% sum

totHDE = huU + huD

chU = deUPCount[grepl('RNA.Ch', names(deUPCount))] %>% sum
chD = deDOWNCount[grepl('RNA.Ch', names(deDOWNCount))] %>% sum

totCDE = chU + chD

# Epigenome
kozDA = import_list('KOZLENKOV_2020/pnas.2011884117.sd11.xlsx')

daAllCount = sapply(kozDA, function(x){nrow(x)})

totHDA = daAllCount[grepl('Hu', names(daAllCount))] %>% sum
totCDA = daAllCount[grepl('Ch', names(daAllCount))] %>% sum

prop.test(c(totHDA, totCDA), c(sum(totHDA, totHDE), sum(totCDA, totCDE)))

# PLOT
toplot = data.frame(vals = c('HS/CS', 'HS/CS'), assay = c('Transcriptome', 'Epigenome'), ratios = c(totHDE/totCDE, totHDA/totCDA))

rio::export(toplot, 'HS_CS_Epigenome_Transcriptome_Total.xlsx')

toplot_koz = rio::import('HS_CS_Epigenome_Transcriptome_Total.xlsx')
toplot_koz$type = 'Kozlenkov'

toplot_our = rio::import('HS_CS_Epigenome_Transcriptome_Total.xlsx')
toplot_our$type = 'This_Study'

toplot = rbind(toplot_koz, toplot_our)
toplot$type = factor(toplot$type, levels = c('This_Study', 'Kozlenkov'))

pdf('HS_CS_Ours_vs_Kozlenkov.pdf', width = 18, height = 3)
ggbarplot(toplot, y = 'ratios', x = 'assay', size = 5, fill = 'assay', color = 'white', palette = c('grey', 'black')) +
facet_wrap(~type) +
theme(text=element_text(size=20)) + ylab('Comparisons') +
ylab('(# of HS) / (# of CS)') +
xlab('') +
NoLegend() +
coord_flip()
dev.off()














