rm(list = ls())
library(dplyr)
library(tidyverse)
library(tidyr)
library(Seurat)
library(ggpubr)
library(paletteer)
library(ggrastr)
library(GeneOverlap)
library(reshape2)
library('bedr')
library('rtracklayer')
source("utility_functions.R")

####
## Find TSS regions
####

gtfHuman = import('HomSap_GRCh38/genes/genes.gtf')
gtfHuman = gtfHuman[gtfHuman$gene_biotype == 'protein_coding',]
gtfHumanT = gtfHuman[gtfHuman$type == 'transcript',]
gtfHumanTPos = gtfHumanT[strand(gtfHumanT) == '+',]
gtfHumanTNeg = gtfHumanT[strand(gtfHumanT) == '-',]

# TSS sites
gtfHumanTPos = resize(gtfHumanTPos, width=2, fix='start')
gtfHumanTNeg = resize(gtfHumanTNeg, width=2, fix='end')

####
## EXPAND TSS
####

# Expansion around the TSS up to 1Mbp on both sites
expandWidth = c(seq(1,1001,100), seq(2000, 10000, 1000), seq(20000,100000,10000), seq(200000,2000000,100000))
gtfTSSL = list()
for(i in 1:length(expandWidth)){

	# Expand the TSS, combine positive and negative TSS and sort for overlap
	gtfPosTSS = resize(gtfHumanTPos, width=expandWidth[i], fix='center')
	gtfNegTSS = resize(gtfHumanTNeg, width=expandWidth[i], fix='center')

	gtfTSS = as.data.frame(c(gtfPosTSS, gtfNegTSS))
	gtfTSS = gtfTSS[,c('seqnames', 'start', 'end', 'gene_name')]
	gtfTSS$seqnames = as.character(gtfTSS$seqnames)
	
	# If expansion is outside of the chromosome, set to 0
	gtfTSS$start = ifelse(gtfTSS$start < 0, 0, gtfTSS$start)
	
	gtfTSS = bedr.sort.region(gtfTSS, verbose = F)
	gtfTSSL[[i]] = gtfTSS
	
	print(i)
}


####
## OVERLAP TSS WITH DARs
####

# Read and convert dars to bed
dars = readRDS('PSEUDOBULK_DARs_ALL.RDS')
darsUniq = unique(dars$Gene)
tmp = gsub(':|-', '_', darsUniq) %>% strsplit(., '_') %>% do.call(rbind,.) %>% as.data.frame
tmp$V2 = as.numeric(tmp$V2)
tmp$V3 = as.numeric(tmp$V3)
darsUniq = cbind(tmp, darsUniq) # retain the compact form
darsUniq = bedr.sort.region(darsUniq, verbose = F)

hsdars = dars[dars$Evolution == 'Human_Specific',]

# Read DEGs
degs = readRDS('PSEUDOBULK_DEGs_ALL.RDS')
hsdegs = degs[degs$Evolution == 'Human_Specific',]
hsdegs$CellType = gsub('OL', 'Oligodendrocyte', hsdegs$CellType)
degs$CellType = gsub('OL', 'Oligodendrocyte', degs$CellType)

# Find CRE - TSS overlaps
st = length(peakOvL) + 1

#peakOvL = list()
for(i in 1:length(gtfTSSL)){

	peakOvL[[i]] = bedr(input = list(a = darsUniq, b = gtfTSSL[[i]]), method = "intersect", params = "-loj", verbose = F)
	print(i)
	
}

####
## PLOT HS / NS
####

ctypes = unique(degs$CellType)
peakOvL3 = readRDS('peakOvL_FINAL.RDS')

expandWidth = c(seq(1,1001,100), seq(2000, 10000, 1000), seq(20000,100000,10000), seq(200000,2000000,100000))

# HUMAN UP - HUMAN OPEN
orUpOpen = list()
hsRatUpOpen = list()
for(i in 1:length(peakOvL3)){

	peakOv = peakOvL3[[i]]
	orL = list()
	hsRatL = list()
	for(j in 1:length(ctypes)){
	
		ctype = ctypes[j]
		
		# HS overlap
		darsCT = hsdars[hsdars$CellType == ctype & hsdars$Regulation == 'Human_UP', 'Gene']
		degsCT = hsdegs[hsdegs$CellType == ctype & hsdegs$Regulation == 'Human_UP', 'Gene']
		degsNSCT = degs[degs$CellType == ctype & degs$Evolution == 'NS', 'Gene']
		
		darGenes = peakOv[peakOv$darsUniq %in% darsCT, 'gene_name'] %>% unique
		hsRat = sum(darGenes %in% degsCT) / length(degsCT)
		nsRat = sum(darGenes %in% degsNSCT) / length(degsNSCT)
		orL[[j]] = (hsRat / nsRat)
		hsRatL[[j]] = hsRat
	}
	
	tmp = unlist(orL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	orUpOpen[[i]] = mean(tmp)
	
	tmp = unlist(hsRatL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	hsRatUpOpen[[i]] = mean(tmp)
	print(i)
}



# HUMAN DOWN - HUMAN CLOSE
orDownClose = list()
hsRatDownClose = list()
for(i in 1:length(peakOvL3)){

	peakOv = peakOvL3[[i]]
	orL = list()
	hsRatL = list()
	for(j in 1:length(ctypes)){
	
		ctype = ctypes[j]
		
		# HS overlap
		darsCT = hsdars[hsdars$CellType == ctype & hsdars$Regulation == 'Human_DOWN', 'Gene']
		degsCT = hsdegs[hsdegs$CellType == ctype & hsdegs$Regulation == 'Human_DOWN', 'Gene']
		degsNSCT = degs[degs$CellType == ctype & degs$Evolution == 'NS', 'Gene']
		
		darGenes = peakOv[peakOv$darsUniq %in% darsCT, 'gene_name'] %>% unique
		hsRat = sum(darGenes %in% degsCT) / length(degsCT)
		nsRat = sum(darGenes %in% degsNSCT) / length(degsNSCT)
		orL[[j]] = (hsRat / nsRat)
		hsRatL[[j]] = hsRat
	}
	
	tmp = unlist(orL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	orDownClose[[i]] = mean(tmp)
	
	tmp = unlist(hsRatL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	hsRatDownClose[[i]] = mean(tmp)
	print(i)
}



# HUMAN DOWN - HUMAN OPEN
orDownOpen = list()
hsRatDownOpen = list()
for(i in 1:length(peakOvL3)){

	peakOv = peakOvL3[[i]]
	orL = list()
	hsRatL = list()
	for(j in 1:length(ctypes)){
	
		ctype = ctypes[j]
		
		# HS overlap
		darsCT = hsdars[hsdars$CellType == ctype & hsdars$Regulation == 'Human_UP', 'Gene']
		degsCT = hsdegs[hsdegs$CellType == ctype & hsdegs$Regulation == 'Human_DOWN', 'Gene']
		degsNSCT = degs[degs$CellType == ctype & degs$Evolution == 'NS', 'Gene']
		
		darGenes = peakOv[peakOv$darsUniq %in% darsCT, 'gene_name'] %>% unique
		hsRat = sum(darGenes %in% degsCT) / length(degsCT)
		nsRat = sum(darGenes %in% degsNSCT) / length(degsNSCT)
		orL[[j]] = (hsRat / nsRat)
		hsRatL[[j]] = hsRat
	}
	
	tmp = unlist(orL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	orDownOpen[[i]] = mean(tmp)
	
	tmp = unlist(hsRatL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	hsRatDownOpen[[i]] = mean(tmp)
	print(i)
}



# HUMAN UP - HUMAN CLOSE
orUpClose = list()
hsRatUpClose = list()
for(i in 1:length(peakOvL3)){

	peakOv = peakOvL3[[i]]
	orL = list()
	hsRatL = list()
	for(j in 1:length(ctypes)){
	
		ctype = ctypes[j]
		
		# HS overlap
		darsCT = hsdars[hsdars$CellType == ctype & hsdars$Regulation == 'Human_DOWN', 'Gene']
		degsCT = hsdegs[hsdegs$CellType == ctype & hsdegs$Regulation == 'Human_UP', 'Gene']
		degsNSCT = degs[degs$CellType == ctype & degs$Evolution == 'NS', 'Gene']
		
		darGenes = peakOv[peakOv$darsUniq %in% darsCT, 'gene_name'] %>% unique
		hsRat = sum(darGenes %in% degsCT) / length(degsCT)
		nsRat = sum(darGenes %in% degsNSCT) / length(degsNSCT)
		orL[[j]] = (hsRat / nsRat)
		hsRatL[[j]] = hsRat
	}
	
	tmp = unlist(orL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	orUpClose[[i]] = mean(tmp)
	
	tmp = unlist(hsRatL)
	tmp = tmp[!(is.nan(tmp))]
	tmp = tmp[!(is.infinite(tmp))]
	hsRatUpClose[[i]] = mean(tmp)
	print(i)
}


tmp = orUpClose
tmp2 = orDownOpen
orUpClose = tmp2
orDownOpen = tmp

# PLOT ALL
dfUpOpen = data.frame(or = unlist(orUpOpen), hsrat = unlist(hsRatUpOpen), width = round(expandWidth / 2), type = 'UP_OPEN')
dfDownClose = data.frame(or = unlist(orDownClose), hsrat = unlist(hsRatDownClose), width = round(expandWidth / 2), type = 'DOWN_CLOSE')
dfUpClose = data.frame(or = unlist(orUpClose), hsrat = unlist(hsRatUpClose), width = round(expandWidth / 2), type = 'UP_CLOSE')
dfDownOpen = data.frame(or = unlist(orDownOpen), hsrat = unlist(hsRatDownOpen), width = round(expandWidth / 2), type = 'DOWN_OPEN')

toplot = do.call(rbind, list(dfUpOpen, dfDownClose, dfUpClose, dfDownOpen))
toplot$type = factor(toplot$type, levels = c('UP_OPEN', 'DOWN_CLOSE', 'UP_CLOSE', 'DOWN_OPEN'))
toplot$width2 = ifelse(toplot$width < 1000, paste0(toplot$width, 'bp'), paste0(toplot$width/1000, 'kb'))
toplot$width3 = ifelse(toplot$width2 %in% c('1kb', '10kb', '100kb', '500kb', '1000kb'), toplot$width2, '')

pdf('DAR_DEG_OVERLAP_TSS_DISTANCE.pdf', height = 5, width = 20)
ggbarplot(toplot, x = 'width2', y = 'or', color = 'type', fill = 'type') +
theme(text = element_text(size=20, face = 'bold'), axis.text.x = element_text(size=7, face = 'bold'), legend.pos = 'right') +
xlab('') + ylab('Odds Ratio (HS/NS)') +
geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
facet_wrap(~type, nrow = 1) +
NoLegend() +
rotate_x_text(90)
dev.off()

pdf('HS_DEG_RATIO_TSS_DISTANCE.pdf', height = 5, width = 20)
ggbarplot(toplot, x = 'width2', y = 'hsrat', color = 'type', fill = 'type') +
theme(text = element_text(size=20, face = 'bold'), axis.text.x = element_text(size=7, face = 'bold'), legend.pos = 'right') +
xlab('') + ylab('Ratio of HS-DEGs with HS-DAR') +
facet_wrap(~type, nrow = 1) +
NoLegend() +
rotate_x_text(90)
dev.off()



####
## FIND AND SCORE DEG-DAR IN THE SAME DIRECTION
####

expandWidth = c(seq(1,1001,100), seq(2000, 10000, 1000), seq(20000,100000,10000), seq(200000,2000000,100000))

# HUMAN UP - HUMAN OPEN
peakOv = peakOvL3[[38]]

dfL3 = list()
for(j in 21:length(ctypes)){

	ctype = ctypes[j]
	
	darsCT_ALL = dars[dars$CellType == ctype, 'Gene']
	
	# HS overlap
	darsCT = hsdars[hsdars$CellType == ctype & hsdars$Regulation == 'Human_UP', 'Gene']
	degsCT = hsdegs[hsdegs$CellType == ctype & hsdegs$Regulation == 'Human_UP', 'Gene']
	degsNSCT = degs[degs$CellType == ctype & degs$Evolution == 'NS', 'Gene']
	
	dfL2 = list()
	for(i in 1:length(degsCT)){
		
		hsgene = degsCT[i]
		qw = peakOv[peakOv$gene_name %in% hsgene, 'darsUniq'] %>% unique
		qw = intersect(qw, darsCT)
		
		if(length(qw) == 0){next}
		
		# Loop over the HS-CREs CONNECTED TO THE CONCOMITANT GENE EXPRESSION CHANGE. COMPUTE SCORE FOR EACH CRE
		dfL = list()
		for(k in 1:length(qw)){
			q1 = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, c('HC_logFC', 'HM_logFC')] %>% as.numeric %>% min
			q1 = 2^q1
			q2 = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, c('HC_logCPM', 'HM_logCPM')] %>% as.numeric %>% min
			
			dargns = peakOv[peakOv$darsUniq == qw[k], 'gene_name'] %>% unique
			degCount = sum(dargns %in% degsCT)
			
			# Compute the score and create the data frame
			score = q1*q2 / degCount
			dfL[[k]] = data.frame(CRE = qw[k], Gene = hsgene, Regulation = 'Human_UP',
						cre_hcLogFC = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HC_logFC'],
						cre_hmLogFC = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HM_logFC'],
						cre_hcLogCPM = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HC_logCPM'],
						cre_hmLogCPM = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HM_logCPM'],
						creDegHit = degCount, CellType = ctype, score = score)
		}
		
		dfL2[[i]] = do.call(rbind, dfL)
		print(paste0('%', round(i/length(degsCT) * 100, digits = 2)))
	}

	print(paste0('Finished CellType: ', ctype))
	dfL3[[j]] = do.call(rbind, dfL2)
}

allfinal = do.call(rbind, dfL3)

# HUMAN DOWN - HUMAN CLOSE
peakOv = peakOvL3[[38]]

dfL3 = list()
for(j in 1:length(ctypes)){

	ctype = ctypes[j]
	
	darsCT_ALL = dars[dars$CellType == ctype, 'Gene']
	
	# HS overlap
	darsCT = hsdars[hsdars$CellType == ctype & hsdars$Regulation == 'Human_DOWN', 'Gene']
	degsCT = hsdegs[hsdegs$CellType == ctype & hsdegs$Regulation == 'Human_DOWN', 'Gene']
	degsNSCT = degs[degs$CellType == ctype & degs$Evolution == 'NS', 'Gene']
	
	dfL2 = list()
	for(i in 1:length(degsCT)){
		
		hsgene = degsCT[i]
		qw = peakOv[peakOv$gene_name %in% hsgene, 'darsUniq'] %>% unique
		qw = intersect(qw, darsCT)
		
		if(length(qw) == 0){next}
		
		# Loop over the HS-CREs CONNECTED TO THE CONCOMITANT GENE EXPRESSION CHANGE. COMPUTE SCORE FOR EACH CRE
		dfL = list()
		for(k in 1:length(qw)){
			q1 = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, c('HC_logFC', 'HM_logFC')] %>% as.numeric %>% min
			q1 = 2^q1
			q2 = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, c('HC_logCPM', 'HM_logCPM')] %>% as.numeric %>% min
			
			dargns = peakOv[peakOv$darsUniq == qw[k], 'gene_name'] %>% unique
			degCount = sum(dargns %in% degsCT)
			
			# Compute the score and create the data frame
			score = abs(q1)*abs(q2) / degCount
			dfL[[k]] = data.frame(CRE = qw[k], Gene = hsgene, Regulation = 'Human_UP',
						cre_hcLogFC = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HC_logFC'],
						cre_hmLogFC = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HM_logFC'],
						cre_hcLogCPM = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HC_logCPM'],
						cre_hmLogCPM = dars[grepl(qw[k], dars$Gene) & dars$CellType == ctype, 'HM_logCPM'],
						creDegHit = degCount, CellType = ctype, score = score)
		}
		
		dfL2[[i]] = do.call(rbind, dfL)
		print(paste0('%', round(i/length(degsCT) * 100, digits = 2)))
	}

	print(paste0('Finished CellType: ', ctype))
	dfL3[[j]] = do.call(rbind, dfL2)
}

allfinal = do.call(rbind, dfL3)


####
## SIGNIFICANCE OF OVERLAP IN 10KB WINDOW
####

# Read and convert dars to bed
dars = readRDS('PSEUDOBULK_DARs_ALL.RDS')

# Read DEGs
degs = readRDS('PSEUDOBULK_DEGs_ALL.RDS')

# Peak to gene overlap
peakOvL3 = readRDS('peakOvL_FINAL.RDS')

# Promoter and near enhancers (10kb around TSS)
peakOv = peakOvL3[[21]]

## HS ##
ctypes = unique(degs$CellType)
resdfL = list()
for(i in 1:length(ctypes)){

	hsUP_CRE = dars[dars$Regulation == 'Human_UP' & dars$CellType == ctypes[i], 'Gene']
	hsUP_CRE_G = peakOv[peakOv$darsUniq %in% hsUP_CRE, 'gene_name'] %>% unique
	hsUP_DEG = degs[degs$Regulation == 'Human_UP' & degs$CellType == ctypes[i], 'Gene']

	hsDOWN_CRE = dars[dars$Regulation == 'Human_DOWN' & dars$CellType == ctypes[i], 'Gene']
	hsDOWN_CRE_G = peakOv[peakOv$darsUniq %in% hsDOWN_CRE, 'gene_name'] %>% unique
	hsDOWN_DEG = degs[degs$Regulation == 'Human_DOWN' & degs$CellType == ctypes[i], 'Gene']

	# Calculate background
	gnsTested = degs[degs$CellType == ctypes[i], 'Gene']
	cresTested = dars[dars$CellType == ctypes[i], 'Gene']

	bcg = length(union(gnsTested,
		peakOv[peakOv$darsUniq %in% cresTested, 'gene_name']))

	hsDEGL = list(hsUP_DEG, hsDOWN_DEG)
	hsCREL = list(hsUP_CRE_G, hsDOWN_CRE_G)
	names(hsDEGL) = c('DEG_UP', 'DEG_DOWN')
	names(hsCREL) = c('CRE_UP_GENE', 'CRE_DOWN_GENE')

	# Enrichment
	resgom = newGOM(hsDEGL, hsCREL, genome.size=bcg)
	pvalmat = getMatrix(resgom, name="pval")
	pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
	pvalmelt = melt(pvalmat, value.name = 'pval')
	pvalmelt$log10_FDR = -log10(pvalmelt$pval)
	oddsmat = getMatrix(resgom, name="odds.ratio")
	oddsmelt = melt(oddsmat, value.name = 'odds.ratio')
	resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio)
	resdf$CellType = ctypes[i]

	resdfL[[i]] = resdf
	
	print(i)
}

resdfAll = do.call(rbind, resdfL)
resdfAll$log10_round = round(resdfAll$log10_FDR, digits = 2)
resdfAll$OR_round = round(resdfAll$OR, digits = 2)
resdfAll$log10_round = ifelse(resdfAll$log10_round > 3, 3, resdfAll$log10_round)
resdfAll$OR_round = ifelse(resdfAll$OR_round > 5, 5, resdfAll$OR_round)

pdf('HS_DAR_DEG_ENRICH_FISHER_10KB_AROUND_TSS.pdf', height = 5, width = 25)
ggscatter(resdfAll, x = 'Var1', y = 'Var2', color = 'log10_round', size = 'OR_round') +
  facet_wrap(~CellType, nrow = 2) +
  labs(x="", y="") + scale_size_continuous(range = c(2,6)) +
  scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=16, face = 'bold')) +
  rotate_x_text(45)
dev.off()




## CS ##
ctypes = unique(degs$CellType)
resdfL = list()
for(i in 1:length(ctypes)){

	hsUP_CRE = dars[dars$Regulation == 'Chimp_UP' & dars$CellType == ctypes[i], 'Gene']
	hsUP_CRE_G = peakOv[peakOv$darsUniq %in% hsUP_CRE, 'gene_name'] %>% unique
	hsUP_DEG = degs[degs$Regulation == 'Chimp_UP' & degs$CellType == ctypes[i], 'Gene']

	hsDOWN_CRE = dars[dars$Regulation == 'Chimp_DOWN' & dars$CellType == ctypes[i], 'Gene']
	hsDOWN_CRE_G = peakOv[peakOv$darsUniq %in% hsDOWN_CRE, 'gene_name'] %>% unique
	hsDOWN_DEG = degs[degs$Regulation == 'Chimp_DOWN' & degs$CellType == ctypes[i], 'Gene']

	# Calculate background
	gnsTested = degs[degs$CellType == ctypes[i], 'Gene']
	cresTested = dars[dars$CellType == ctypes[i], 'Gene']

	bcg = length(union(gnsTested,
		peakOv[peakOv$darsUniq %in% cresTested, 'gene_name']))

	hsDEGL = list(hsUP_DEG, hsDOWN_DEG)
	hsCREL = list(hsUP_CRE_G, hsDOWN_CRE_G)
	names(hsDEGL) = c('DEG_UP', 'DEG_DOWN')
	names(hsCREL) = c('CRE_UP_GENE', 'CRE_DOWN_GENE')

	# Enrichment
	resgom = newGOM(hsDEGL, hsCREL, genome.size=bcg)
	pvalmat = getMatrix(resgom, name="pval")
	pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
	pvalmelt = melt(pvalmat, value.name = 'pval')
	pvalmelt$log10_FDR = -log10(pvalmelt$pval)
	oddsmat = getMatrix(resgom, name="odds.ratio")
	oddsmelt = melt(oddsmat, value.name = 'odds.ratio')
	resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio)
	resdf$CellType = ctypes[i]

	resdfL[[i]] = resdf
	
	print(i)
}

resdfAll = do.call(rbind, resdfL)


resdfAll$log10_round = round(resdfAll$log10_FDR, digits = 2)
resdfAll$OR_round = round(resdfAll$OR, digits = 2)
resdfAll$log10_round = ifelse(resdfAll$log10_round > 3, 3, resdfAll$log10_round)
resdfAll$OR_round = ifelse(resdfAll$OR_round > 5, 5, resdfAll$OR_round)

pdf('CS_DAR_DEG_ENRICH_FISHER_10KB_AROUND_TSS.pdf', height = 5, width = 25)
ggscatter(resdfAll, x = 'Var1', y = 'Var2', color = 'log10_round', size = 'OR_round') +
  facet_wrap(~CellType, nrow = 2) +
  labs(x="", y="") + scale_size_continuous(range = c(2,6)) +
  scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=16, face = 'bold')) +
  rotate_x_text(45)
dev.off()


## MS ##
ctypes = unique(degs$CellType)
resdfL = list()
for(i in 1:length(ctypes)){

	hsUP_CRE = dars[dars$Regulation == 'Macaque_UP' & dars$CellType == ctypes[i], 'Gene']
	hsUP_CRE_G = peakOv[peakOv$darsUniq %in% hsUP_CRE, 'gene_name'] %>% unique
	hsUP_DEG = degs[degs$Regulation == 'Macaque_UP' & degs$CellType == ctypes[i], 'Gene']

	hsDOWN_CRE = dars[dars$Regulation == 'Macaque_DOWN' & dars$CellType == ctypes[i], 'Gene']
	hsDOWN_CRE_G = peakOv[peakOv$darsUniq %in% hsDOWN_CRE, 'gene_name'] %>% unique
	hsDOWN_DEG = degs[degs$Regulation == 'Macaque_DOWN' & degs$CellType == ctypes[i], 'Gene']

	# Calculate background
	gnsTested = degs[degs$CellType == ctypes[i], 'Gene']
	cresTested = dars[dars$CellType == ctypes[i], 'Gene']

	bcg = length(union(gnsTested,
		peakOv[peakOv$darsUniq %in% cresTested, 'gene_name']))

	hsDEGL = list(hsUP_DEG, hsDOWN_DEG)
	hsCREL = list(hsUP_CRE_G, hsDOWN_CRE_G)
	names(hsDEGL) = c('DEG_UP', 'DEG_DOWN')
	names(hsCREL) = c('CRE_UP_GENE', 'CRE_DOWN_GENE')

	# Enrichment
	resgom = newGOM(hsDEGL, hsCREL, genome.size=bcg)
	pvalmat = getMatrix(resgom, name="pval")
	pvalmat = apply(pvalmat, 2, function(x){p.adjust(x, method = 'fdr')})
	pvalmelt = melt(pvalmat, value.name = 'pval')
	pvalmelt$log10_FDR = -log10(pvalmelt$pval)
	oddsmat = getMatrix(resgom, name="odds.ratio")
	oddsmelt = melt(oddsmat, value.name = 'odds.ratio')
	resdf = cbind(pvalmelt, OR = oddsmelt$odds.ratio)
	resdf$CellType = ctypes[i]

	resdfL[[i]] = resdf
	
	print(i)
}

resdfAll = do.call(rbind, resdfL)


resdfAll$log10_round = round(resdfAll$log10_FDR, digits = 2)
resdfAll$OR_round = round(resdfAll$OR, digits = 2)
resdfAll$log10_round = ifelse(resdfAll$log10_round > 3, 3, resdfAll$log10_round)
resdfAll$OR_round = ifelse(resdfAll$OR_round > 5, 5, resdfAll$OR_round)

pdf('MS_DAR_DEG_ENRICH_FISHER_10KB_AROUND_TSS.pdf', height = 5, width = 25)
ggscatter(resdfAll, x = 'Var1', y = 'Var2', color = 'log10_round', size = 'OR_round') +
  facet_wrap(~CellType, nrow = 2) +
  labs(x="", y="") + scale_size_continuous(range = c(2,6)) +
  scale_color_gradient2(midpoint = 1, low = 'blue', high = 'red')+
  theme_classic() +
  theme(text = element_text(size=16, face = 'bold')) +
  rotate_x_text(45)
dev.off()

####
## RATIO OF HS-GENES in HS-CREs
####

degs = readRDS('PSEUDOBULK_DEGs_ALL.RDS')
ovs = import('CRE_Gene_Association_500kb.xlsx')

ctypes = unique(degs$CellType)
ovL = list()
for(i in 1:length(ctypes)){
	tmp1 = degs[degs$CellType == ctypes[i] & degs$Regulation == 'Human_UP', 'Gene'] %>% unique
	tmp2 = ovs[ovs$CellType == ctypes[i] & ovs$Regulation == 'Human_UP', 'Gene'] %>% unique

	ovL[[i]] = sum(tmp1 %in% tmp2)/length(tmp1)
}

mean(unlist(ovL))

# Ratio of DARs
dars = readRDS('PSEUDOBULK_DARs_ALL.RDS')
ctypes = unique(dars$CellType)
ovL = list()
for(i in 1:length(ctypes)){
	tmp1 = dars[dars$CellType == ctypes[i] & dars$Regulation == 'Human_UP', 'Gene'] %>% unique
	tmp2 = ovs[ovs$CellType == ctypes[i] & ovs$Regulation == 'Human_UP', 'CRE'] %>% unique

	ovL[[i]] = sum(tmp1 %in% tmp2)/length(tmp1)
}

mean(unlist(ovL))







