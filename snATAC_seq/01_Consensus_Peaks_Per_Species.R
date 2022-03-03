rm(list = ls())
library("bedr")
library(ggpubr)
options(scipen = 999)

####
## This script finds consensus pekas across samples of the same species
####


# Read peaks identified per sample
hum_p1 = read.table("hsa1_simplepeaks.bed", stringsAsFactors = F)
hum_p2 = read.table("hsa2_simplepeaks.bed", stringsAsFactors = F)
hum_p3 = read.table("hsa3_simplepeaks.bed", stringsAsFactors = F)
hum_p4 = read.table("hsa4_simplepeaks.bed", stringsAsFactors = F)

# Read peaks identified after pooling samples
hum_pool = read.table("hsa_all_simplepeaks.bed", stringsAsFactors = F)


# This function intersects two regions and returns the overlap even if there is zero overlap.
bedr_ov = function(merged, initial){

	# Get intersecting peaks
	overlap <- bedr(input = list(a = merged, b = initial), 
		       method = "intersect", 
		       params = "-sorted -loj")

	# Convert positions to numeric
	overlap[, c(2,3,5,6)] = lapply(c(2,3,5,6), function(x){as.numeric(overlap[,x])})

	# Calculate overlap
	overlap$ov =  as.numeric(pmin(overlap[,3], overlap[,6])) - as.numeric(pmax(overlap[,2], overlap[,5]))
	overlap[overlap$ov < 0, "ov"] = 0

	# Calculate overlap percentage
	overlap$ratio = ifelse((overlap[,6] - overlap[,5]) > 0, overlap$ov/(overlap[,6] - overlap[,5]), 0)


	# Keep best overlapping match
	overlap$peak = paste(overlap$V1, overlap$V2, overlap$V3, sep = "-")
	overlap = overlap[order(overlap$peak, -overlap$ratio),]
	overlap = overlap[!duplicated(overlap$peak),]

	return(overlap)
}


# Calculate pooled peaks that overlap peaks from individual samples.
hum_inter1 = bedr_ov(hum_pool, hum_p1)
hum_inter2 = bedr_ov(hum_pool, hum_p2)
hum_inter3 = bedr_ov(hum_pool, hum_p3)
hum_inter4 = bedr_ov(hum_pool, hum_p4)

hum_pool$peak = paste(hum_pool$V1, hum_pool$V2, hum_pool$V3, sep = "-")
hum_pool = hum_pool[match(hum_inter1$peak, hum_pool$peak),] # match order of the peaks


hum_pool_cons = hum_pool[(hum_inter1$ratio > 0.5 & hum_inter2$ratio > 0.5 &
			  hum_inter3$ratio > 0.5) |
			 (hum_inter1$ratio > 0.5 & hum_inter2$ratio > 0.5 &
			  hum_inter4$ratio > 0.5) |
			 (hum_inter2$ratio > 0.5 & hum_inter3$ratio > 0.5 &
			  hum_inter4$ratio > 0.5) |
			(hum_inter1$ratio > 0.5 & hum_inter3$ratio > 0.5 &
			  hum_inter4$ratio > 0.5) ,c(1:3)]

# Write consensus peaks
write.table(hum_pool_cons,
		"human_consensus.bed",
		col.names = F,
		row.names = F,
		quote = F,
		sep = "\t")


# Barplot of peak numbers
vars = c("Human_1", "Human_2", "Human_3", "Human_4", "Pooled", "Human1_in_Pooled", "Human2_in_Pooled", "Human3_in_Pooled", "Human4_in_Pooled", "Consensus")
vals = c(nrow(hum_p1), nrow(hum_p2), nrow(hum_p3), nrow(hum_p4), nrow(hum_pool), sum(hum_inter1$ratio > 0.5), sum(hum_inter2$ratio > 0.5), sum(hum_inter3$ratio > 0.5), sum(hum_inter4$ratio > 0.5), nrow(hum_pool_cons))

toplot = data.frame(var = vars, val = vals)

toplot$var = factor(toplot$var, levels = c("Human_1", "Human_2", "Human_3", "Human_4", "Pooled", "Human1_in_Pooled", "Human2_in_Pooled", "Human3_in_Pooled", "Human4_in_Pooled", "Consensus"))

pdf("Human_consensus_peaks.pdf")
ggbarplot(toplot, x = 'var', y = 'val', fill = 'var') + rotate_x_text(45)
dev.off()

