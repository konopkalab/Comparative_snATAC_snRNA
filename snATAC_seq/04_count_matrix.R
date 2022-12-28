rm(list = ls())
library(reshape2)
library(data.table)
library(plyr)
library(Matrix)
library(parallel)

####
## Custom approach (split by cell and append to file. column number stays the same)
####

args = commandArgs(trailingOnly = TRUE)
sp = 'human'
samp = 'h1'

ov_peak = fread(paste0("consensus_peaks/", sp, samp, "_peak_ov.bed"), stringsAsFactors = F)
trbc = read.table("atac_cellranger_trusted_bcs.txt")
ov_peak = as.data.frame(ov_peak)

# Metadata Peak #

# Total reads in peaks (rip).
rip = table(ov_peak$V1)

# Get total reads per barcode.
barc_count = fread(paste0("consensus_peaks/", sp,  samp, "_reads.bed"), stringsAsFactors = F)
barc_count = as.data.frame(barc_count)
barc_table = table(barc_count$V5)
barc_table_filt = barc_table[names(barc_table) %in% names(rip)]

# Update metrics so all barcodes are found in total barcodes
rip = rip[match(names(barc_table_filt), names(rip))]

# frip: fraction of reads in peaks
frip = rip / barc_table_filt

# Metadata
metadata = data.frame(barcode = names(barc_table_filt),
		      total = as.numeric(barc_table_filt),
		      rip = as.numeric(rip),
		      frip = as.numeric(frip))

# Keep barcodes > 1000 total and trusted
trusted = gsub("_.*", "", metadata$barcode) %in% trbc$V1
metadata = metadata[metadata$total > 1000 & trusted, ]
ov_peak = ov_peak[ov_peak$V1 %in% metadata$barcode, ]

writedir = paste0(sp, samp, "/")
dir.create(writedir)
write.table(metadata, file = paste0(writedir, sp, samp, "_metadata.tsv"), sep = "\t", quote = F, row.names = F)

# Release RAM
rm(barc_count, barc_table, frip, metadata)
gc()

# Reshape and aggregate fragment overlap per cell per peak
ov_peak$V5 = paste0(ov_peak$V2, ":", ov_peak$V3, "-", ov_peak$V4)
ov_peak = ov_peak[, -c(2,3,4)]
ov_peak$V6 = 1
ov_peak = data.table(ov_peak)
ov_peak = ov_peak[,list(V6 = sum(V6)), by = 'V1,V5']
ov_peak = as.data.frame(ov_peak)

# Split into cells per peak
ov_split = split(ov_peak, ov_peak$V1)

# Reshape into vector with peak in row, cells in column
ov_matlist = mclapply(ov_split, function(x){tmp = t(x[,3]);
					  colnames(tmp) = x[,2];		
					  rownames(tmp) = x[1,1];tmp
					 })

# Create a pseudo cell that has all the peaks. This will make sure each matrix of 10k cells has all the peaks.
peakun = unique(ov_peak$V5)
tmpdf = data.frame(rep(0, length(peakun)))
tmpdf = t(tmpdf)
colnames(tmpdf) = peakun
rownames(tmpdf) = "mock"

# Divide each peak vector into 100 small vectors for easier handling.
chunk = 100
ov_submatlist = split(ov_matlist, rep(1:ceiling(length(ov_matlist)/chunk), each = chunk)[1:length(ov_matlist)])

# Write each peak as a separate matrix market file
tmp = mclapply(names(ov_submatlist), mc.cores = 15, function(x){
						  ov_df1 = rbind.fill.matrix(c(list(tmpdf), ov_submatlist[[x]]));
					  	  ov_df1 = ov_df1[-1,];
					  	  ov_df1[is.na(ov_df1)] = 0;
					  	  ov_df1 <- as(ov_df1 , "dgCMatrix");
					 	  rownames(ov_df1) = names(ov_submatlist[[x]]);
					 	  writefil = paste0(writedir, "tmp_",x,"_.mtx");
					 	  writeMM(ov_df1, file = writefil)
						  gc()
						}
	      )


### Merge mtx files ###
# Adjust row numbers
fls = list.files(writedir, pattern = "tmp_")

tmp_long = list()
for(i in fls){

	cf = as.numeric(strsplit(i,"_")[[1]][2])
	tmp = read.table(paste0(writedir, i), skip = 2)
	tmp$V1 = tmp$V1 + ((cf - 1) * chunk)
	tmp_long[[cf]] = tmp
}

# Total non-zero elements
tot = lapply(tmp_long, function(x){nrow(x)})
tot = sum(unlist(tot))

# First line
firstline = c(length(ov_split), length(peakun), tot)

# Merge all
allmtx = do.call(rbind, c(list(firstline), tmp_long))

# Write MatrixMarket header
fc = file(paste0(writedir, "pc_mat.mtx"))
writeLines("%%MatrixMarket matrix coordinate integer general", fc)
close(fc)

# Append the matrix
options(scipen = 999)
write.table(allmtx, file = paste0(writedir, "pc_mat.mtx"), append = T, row.names = F, col.names = F, quote = F)

# Write barcodes
write.table(names(ov_split), file = paste0(writedir, "barcodes.tsv"), row.names = F, col.names = F, quote = F)

# Write peaks
tmp_bed = do.call(rbind, strsplit(peakun, ":"))
tmp_bed2 = do.call(rbind, strsplit(tmp_bed[,2], "-"))
peak_bed = cbind(tmp_bed[,1], tmp_bed2)
write.table(peakun, file = paste0(writedir, "peaks_notab.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(peak_bed, file = paste0(writedir, "peaks.bed"), row.names = F, col.names = F, quote = F, sep = "\t")


