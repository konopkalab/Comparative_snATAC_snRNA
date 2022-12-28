rm(list = ls())
require("rphast")
library(ape)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(parallel)
library(Biostrings)
library(ggpubr)
library(Seurat)
library(parallel)

####
## PREPARE THE DATA
####

# Take chromosome as an input
args = commandArgs(trailingOnly = TRUE)
chr = args[1]
startPoint = as.numeric(args[2])
endPoint = as.numeric(args[3])

# Load MSA
align = readRDS(paste0('HAR_30WAY/MSA/', chr, '_anthropoids_msa.RDS'))

# Read MSA offset
ofsDF = read.table(paste0("UCSC_MULTIZ30_2017/head_", chr, ".maf"), fill = T)
ofs = as.numeric(ofsDF[2,3])

# Load tree topology
treetop = read.table("UCSC_MULTIZ30_2017/treeforR.txt")
treetop = treetop$V1

# Keep only apes and monkeys that are syntenic.
tree = read.tree(text = treetop)
dropspecies1 = c('tarSyr2', 'micMur3', 'proCoq1', 'eulMac1', 'eulFla1', 'otoGar3', 'mm10', 'canFam3', 'dasNov3')
dropspecies2 = c('papAnu3', 'rhiBie1', 'rhiRox1', 'ponAbe2', 'nasLar1', 'panPan2')
dropspecies = union(dropspecies1, dropspecies2)

tree = drop.tip(tree, dropspecies)
treetop = write.tree(tree)
align2 = align[tree$tip,]

# Load peaks
pks = read.feat('/home2/s422159/project/FROM_WORKDIR/pr3/atacseq/human/consensus_peaks/04_consensus_peaks/merged_lifted_human_sorted.bed')
pks = pks[pks[,1] == chr,]
pks$seqname = "hg38"

#####
## KEEP SUFFICIENTLY CONSERVED PEAKS
#####

# Get regions within CREs without gaps in 3 species or in at least 8 species (50% of all species)
hasHCM = informative.regions.msa(align2, 3, c("hg38", "panTro5", "rheMac8"))
hasHalf = informative.regions.msa(align2, 8)

# Adjust peaks by offset
hasHCM[,4] = hasHCM[,4] + ofs
hasHCM[,5] = hasHCM[,5] + ofs
hasHalf[,4] = hasHalf[,4] + ofs
hasHalf[,5] = hasHalf[,5] + ofs

# Check what the coverage is within all peaks
infElm = coverage.feat(pks, hasHalf, hasHCM, get.feats=TRUE)
coverage.feat(infElm)/coverage.feat(pks)

# Split peaks into same lengths
splitLength = 150
splitElm = split.feat(infElm, f=splitLength, drop=TRUE)

# Remove if no room to pad. This usually doesn't remove anything.
splitElm = splitElm[(splitElm[,5] + 25000) < (ncol(align2) + ofs), ]
splitElm = splitElm[(splitElm[,4] - 25000) > (ofs), ]

#####
## RUN ACCELERATION
#####

# Test each peak for human acceleration
tmp = seq(startPoint, endPoint, by = 100)
inters1 = tmp
inters2 = c(tmp[2:length(tmp)], endPoint)
inters2[1:(length(inters2) - 1)] = inters2[1:(length(inters2) - 1)] - 1


print(paste0('TESTING ACCELERATION OF ', nrow(splitElm), ' ELEMENTS'))

for(i in 1:length(inters1)){

	print(paste0(inters1[i], '_', inters2[i]))

	obsPs = mclapply(inters1[i]:inters2[i], mc.cores = 10, function(x){

						# Pad 25kb on each side
						st = splitElm[x,4] - 25000
						end = splitElm[x,5] + 25000

						if(x%%50 == 0){print(x)}

						# Adjust for ofset in alignment
						st = st - ofset
						end = end - ofset

						# Extract, fit neutral model and test acceleration of human branch on the neutral model
						q = sub.msa(align2, start.col = st, end.col = end, pointer.only = T)
						fit = phyloFit(q, tree=treetop, subst.mod="SSREV", EM = T, nrates = 4)
						phyloP(fit, method = 'LRT', msa=align2, mode="ACC", features=splitElm[x,], subtree="hg38")})

	obsP = do.call(rbind, obsPs)
	saveRDS(obsP, paste0(chr, '_obsP_', inters1[i], '_', inters2[i], '.RDS'))

}


####
## MERGE ALL RESULTS
####

# Keep offsets in a data frame
ofsFl = list.files(path = "UCSC_MULTIZ30_2017", pattern = "head", full = T)
offsets = sapply(ofsFl, function(x){ofsDF = read.table(x, fill = T); ofs = as.numeric(ofsDF[2,3]); ofs})
offsetsDF = data.frame(offset = as.numeric(offsets), chromosome = gsub('.*_chr', 'chr', names(offsets)) %>% gsub('.maf', '', .))


fls = list.files(path = '~/workdir', pattern = 'obsPL_chr', full = T)
accDFL = list()
for(j in 1:length(fls)){

	tmp = readRDS(fls[j]) %>% do.call(rbind, .)
	chrom = gsub('.*chr', 'chr', fls[j]) %>% gsub('_[0-9].*', '', .) %>% gsub('.RDS', '', .)
	tmp$chr = chrom
	tmp[,2] = tmp[,2] + offsetsDF[offsetsDF$chromosome  == chrom, 'offset']
	tmp[,3] = tmp[,3] + offsetsDF[offsetsDF$chromosome  == chrom, 'offset']

	accDFL[[j]] = tmp
}

accDF = do.call(rbind, accDFL)
saveRDS(accDF, 'allAccDF.RDS')

####
## PREPARE CREs FOR OVERLAP
####

# Load background peaks
darsAll = readRDS('PSEUDOBULK_DARs_ALL.RDS')
darsAll$peak = darsAll$Gene

# HS-DARs
hsDARs = darsAll[darsAll$Evolution == 'Human_Specific', 'peak'] %>% unique %>%
	gsub(':', '-', .) %>% strsplit(., '-') %>% do.call(rbind, .) %>%
	as.data.frame %>% mutate(V1 = as.character(V1), V2 = as.numeric(V2), V3 = as.numeric(V3))

# NS-DARs
nsCREs = darsAll[darsAll$Evolution == 'NS', 'peak'] %>%
	unique %>% gsub(':', '-', .) %>% strsplit(., '-') %>% do.call(rbind, .) %>%
	as.data.frame %>% mutate(V1 = as.character(V1), V2 = as.numeric(V2), V3 = as.numeric(V3))

# Sort for overlap
hsDARs = bedr.sort.region(hsDARs)
nsCREs = bedr.sort.region(nsCREs)


####
## HAR associations with HS-CREs and background (NS-CREs)
####

# Read Adult HARs
adultHARAll = readRDS('allAccDF.RDS')

# Keep the ones with positive acceleration
loglik = seq(1,10,0.25)
adultHAR_Acc = adultHARAll[adultHARAll$hRat > 1,]

# HAR-HS overlap odds ratio with increasing cutoff
dfL = list()
for(i in 1:length(loglik)){

	# HARs with the given loglik cutoff
	har = adultHAR_Acc[adultHAR_Acc$hRat > loglik[i],]
	har = adultHAR_Acc[sample(1:nrow(adultHAR_Acc), nrow(har)),]
	har = har[, 1:3]
	har = bedr.sort.region(har, verbose = F)
	
	# HAR-HS overlap
	ov = bedr(input = list(a = hsDARs, b = har), method = "intersect", params = "-loj", verbose = F, check.valid = F)
	ov$creCoord = paste0(ov[,1], '_', ov[,2], '_', ov[,3])
	ovCreHar = unique(ov[ov$start.b != '-1', 'creCoord'])
	harHS_Sum = length(ovCreHar)
	harHS_Rat = harHS_Sum / nrow(hsDARs)

	# HAR-NS overlap
	ov = bedr(input = list(a = nsCREs, b = har), method = "intersect", params = "-loj", verbose = F, check.valid = F)
	ov$creCoord = paste0(ov[,1], '_', ov[,2], '_', ov[,3])
	ovCreHar = unique(ov[ov$start.b != '-1', 'creCoord'])
	harNS_Sum = length(ovCreHar)
	harNS_Rat = harNS_Sum / nrow(nsCREs)

	dfL[[i]] = data.frame(HS_Ov = harHS_Sum, NS_Ov = harNS_Sum, HS_Ratio = harHS_Rat,
			NS_Ratio = harNS_Rat, OR = harHS_Rat / harNS_Rat, LogLik = loglik[i])

	print(i)
}

# Combine
toplot = do.call(rbind, dfL)

# Read Published HARs
allHars = rio::import('Doan_Combined_HARs_hg38.xlsx')
allHars = allHars[, c(1:3)]

# HAR-HS overlap in published HARs
ov = bedr(input = list(a = hsDARs, b = allHars), method = "intersect", params = "-loj", verbose = F, check.valid = F)
ov$creCoord = paste0(ov[,1], '_', ov[,2], '_', ov[,3])
ovCreHar = unique(ov[ov$start.b != '-1', 'creCoord'])
harHS_Sum = length(ovCreHar)
harHS_Rat = harHS_Sum / nrow(hsDARs)

# HAR-NS overlap in published HARs
ov = bedr(input = list(a = nsCREs, b = allHars), method = "intersect", params = "-loj", verbose = F, check.valid = F)
ov$creCoord = paste0(ov[,1], '_', ov[,2], '_', ov[,3])
ovCreHar = unique(ov[ov$start.b != '-1', 'creCoord'])
harNS_Sum = length(ovCreHar)
harNS_Rat = harNS_Sum / nrow(nsCREs)

pubHAR = data.frame(HS_Ov = harHS_Sum, NS_Ov = harNS_Sum, HS_Ratio = harHS_Rat,
		NS_Ratio = harNS_Rat, OR = harHS_Rat / harNS_Rat, LogLik = 12)

# Combine cortical HAR and published HAR stats to plot together
toplot2 = rbind(toplot, pubHAR)
toplot2$type = ifelse(toplot2$LogLik == 12, 'Published', 'This_Study')

pdf("HS_HAR_VS_NS_CRE_OddsRatio_WITH_PUBLISHED_HAR.pdf")
ggscatter(toplot2, x = 'LogLik', y = 'OR', color = 'type', size = 4) +
geom_vline(xintercept = 4.774837, linetype = 'dashed', color = 'red') +
theme(text=element_text(size=20, face = 'bold')) +
rotate_x_text(90) +
ylab('(HS & HAR) / (NS & HAR)') +
xlab('Human Sequence Divergence\n(Log likelihood cutoffs)')
dev.off()







