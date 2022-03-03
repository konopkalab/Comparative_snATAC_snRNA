bedPie = function(bedL, fn, pienames = c('Only_First', 'Only_Second', 'Both'), fisherBCG){

	require(scales)
	require('bedr')

	refBED = bedL[[1]]
	refOv1 = bedr(input = list(a = refBED, b = bedL[[2]]), method = "intersect", params = "-loj", verbose = F)
	refOv1$refbed = paste0(refOv1[,1], ':', refOv1[,2], '-', refOv1[,3])

	refOv2 = bedr(input = list(a = refBED, b = bedL[[3]]), method = "intersect", params = "-loj", verbose = F)
	refOv2$refbed = paste0(refOv2[,1], ':', refOv2[,2], '-', refOv2[,3])


	ovALL = length( intersect( unique(refOv1[refOv1[,5] != '-1', 'refbed']), unique(refOv2[refOv2[,5] != '-1', 'refbed']) ) )
	only1 = length( setdiff( unique(refOv1[refOv1[,5] != '-1', 'refbed']), unique(refOv2[refOv2[,5] != '-1', 'refbed']) ) )
	only2 = length( setdiff( unique(refOv2[refOv2[,5] != '-1', 'refbed']), unique(refOv1[refOv1[,5] != '-1', 'refbed']) ) )

	# Compute percentage and position of labels
	toplot = data.frame(vals = c(only1,only2,ovALL), vars = pienames)
	toplot$vals = toplot$vals/sum(toplot$vals) * 100
	toplot <- toplot %>% 
	  arrange(desc(vars)) %>%
	  mutate(prop = vals / sum(toplot$vals) *100) %>%
	  mutate(ypos = cumsum(prop)- 0.5*prop )

	pval = fisher.test(matrix(c(ovALL, only1, only2, fisherBCG),2,2), alternative = 'greater')$p.value
	pval = ifelse(pval < 2.2e-16, 'p-value < 2.2e-16', paste0('p-value = ', round(pval, digits = 2)))

	blank_theme <- theme_minimal()+
	  theme(
	  axis.title.x = element_blank(),
	  axis.title.y = element_blank(),
	  panel.border = element_blank(),
	  panel.grid=element_blank(),
	  axis.ticks = element_blank(),
	  plot.title=element_text(size=14, face="bold")
	  )

	pdf(paste0(fn, '.pdf'), height = 10, width = 10)
	print( ggplot(toplot, aes(x="", y=vals, fill=vars))+
	geom_bar(width = 1, stat = "identity") +
	coord_polar("y", start=0) +
	scale_fill_manual(values = brewer.pal(n = 4, name = "OrRd")[c(4,3,2)]) + blank_theme +
	  theme(axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5))+
	  ggtitle(pval) +
	  geom_text(aes(y = ypos, label = percent(vals/100, accuracy = 1)), size=7, fontface = 'bold') )
	dev.off()

}

GOenrich = function(gns, uni){

	require(org.Hs.eg.db)
	require(clusterProfiler)
	gene_df <- bitr(gns, fromType = "SYMBOL",
		toType = c("ENTREZID"),
		OrgDb = org.Hs.eg.db)

	uni_df <- bitr(uni, fromType = "SYMBOL",
		toType = c("ENTREZID"),
		OrgDb = org.Hs.eg.db)
	
	goterms = c('MF', 'BP', 'CC')
	egos = list()
	for(i in 1:length(goterms)){
		ego <- enrichGO(gene          = gene_df$ENTREZID,
				universe      = uni_df$ENTREZID,
				OrgDb         = org.Hs.eg.db,
				ont           = goterms[i],
				pAdjustMethod = "BH",
				pvalueCutoff  = 0.05,
				qvalueCutoff  = 0.2)
		ego = setReadable(ego, OrgDb = org.Hs.eg.db)
		ego = as.data.frame(ego)
		if(nrow(ego) == 0){next}

		ego$GO = goterms[i]
		egos[[i]] = ego
	}
	egodf = do.call(rbind, egos)
	
	return(egodf)
}

geneOv = function(x,y,bcg){
	require(GeneOverlap)
	go.obj <- newGeneOverlap(x, y, genome.size=bcg)
	go.obj <- testGeneOverlap(go.obj)
	return(go.obj)
}

annotValid = function(seurobj, group1, group2, fn){

	# Calculate total number of cells for given predicted-cluster id match
	filtlabel = seurobj[[c(group1, group2)]]
	colnames(filtlabel) = c('group1', 'group2')
	filtlabel$total = 1
	tmp = aggregate(. ~ group1 + group2, filtlabel, sum)

	# Add total number of cells with the given predicted id
	predid_tab = table(filtlabel$group2)
	matchind = match(tmp$group2, names(predid_tab))
	tmp$group2_tot = predid_tab[matchind]

	# Add total number of cells with the given cluster id
	group1_tab = table(filtlabel$group1)
	matchind = match(tmp$group1, names(group1_tab))
	tmp$group1_tot = group1_tab[matchind]

	# Calculate percentage of label for the given cluster. Keep only the ones that are > 20%
	tmp$match_query_ratio = tmp$total / tmp$group1_tot
	tmp$match_ref_ratio = tmp$total / tmp$group2_tot


	# Heatmap of annotation
	toplot = tmp[, c('group2', 'group1', 'match_query_ratio', 'match_ref_ratio')]
	toplot$group2 = factor(toplot$group2)
	toplot$group1 = factor(toplot$group1)
	toplot$match_query_ratio = as.numeric(as.character(toplot$match_query_ratio))
	toplot$match_ref_ratio = as.numeric(as.character(toplot$match_ref_ratio))

	plt = ggplot(toplot, aes(x=group1, y=group2, color = match_query_ratio, size = match_ref_ratio)) +
		geom_point() +
		scale_fill_gradient2(low = 'white', high = 'blue', midpoint = 0.1, aesthetics = "color") +
		theme_classic() +
		xlab(group1) +
		ylab(group2) +
		theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
			rotate_x_text(90)

	pdf(fn)
	print(plt)
	dev.off()

}



ltpred = function(seurobj, vars = c("predicted.id", "seurat_clusters"), fn, heatmap = T, textsize = 20){

# Calculate total number of cells for given predicted-cluster id match
filtlabel = seurobj[[vars]]
colnames(filtlabel) = c("predicted.id", "seurat_clusters")
filtlabel$total = 1
tmp = aggregate(. ~ predicted.id + seurat_clusters, filtlabel, sum)

# Add total number of cells with the given predicted id
predid_tab = table(filtlabel$predicted.id)
matchind = match(tmp$predicted.id, names(predid_tab))
tmp$predicted_id_tot = predid_tab[matchind]

# Add total number of cells with the given cluster id
seurid_tab = table(filtlabel$seurat_clusters)
matchind = match(tmp$seurat_clusters, names(seurid_tab))
tmp$seurat_id_tot = seurid_tab[matchind]

# Calculate percentage of label for the given cluster. Keep only the ones that are > 20%
tmp$perc_inclus = tmp$total / tmp$seurat_id_tot
tmp2 = tmp[tmp$perc_inclus > 0.05,]


# Heatmap of annotation
toplot = tmp2[, c('predicted.id', 'seurat_clusters', 'perc_inclus')]
toplot$predicted.id = factor(toplot$predicted.id)
toplot$seurat_clusters = factor(toplot$seurat_clusters)
toplot$perc_inclus = as.numeric(as.character(toplot$perc_inclus))

plt1 = ggplot(toplot, aes(seurat_clusters, predicted.id, fill = perc_inclus)) +
	geom_tile(color = "white") +
	geom_text(label = round(toplot$perc_inclus, digits = 2)) +
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
	midpoint = 0.05, limit = c(0,1)) +
	theme(axis.text.x = element_text(size=textsize),
		axis.text.y = element_text(size=textsize),
		axis.title = element_text(size=textsize)) +
		rotate_x_text(90)

plt2 = ggscatter(toplot, x = 'predicted.id', y = 'seurat_clusters', color = 'perc_inclus', size = 10) +
	scale_color_gradient2(low = "white", high = "blue", mid = "white", midpoint = 0, limit = c(0,1)) +
	xlab('Prediction') + ylab('Annotation') +
	theme(axis.text.x = element_text(size=textsize),
		axis.text.y = element_text(size=textsize),
		axis.title = element_text(size=textsize)) +
		rotate_x_text(90) + grids(linetype = "dashed", color = 'grey')

if(heatmap == T){
	pdf(paste0(fn, '.pdf'), width = 20, height =10)
	print(plt1)
	dev.off()
} else{
	pdf(paste0(fn, '.pdf'), width = 10, height =10)
	print(plt2)
	dev.off()
}

}



stackedbarplot = function(meta, groupx, groupfill, fn, horizontal = F, bold = F){

clkeys = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_keys %>% as.data.frame
clkeys$size = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_size 

# Add all cell numbers
spkeys = meta %>% group_by(meta[,groupfill]) %>% group_keys %>% as.data.frame
colnames(spkeys) = groupfill
spkeys$size = meta %>% group_by(meta[, groupfill]) %>% group_size
spkeys = spkeys %>% add_column(tmp = 'AllCells')
colnames(spkeys)[3] = groupx
clkeys = rbind(clkeys, spkeys)

colnames(clkeys) = c('cluster', 'variable', 'value')

plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
	    geom_bar(position = "fill", stat = "identity") +
	    xlab("") +
	    ylab("Percentage") +
	    theme_classic() +
	    #scale_fill_manual(values = colors) +
	    scale_y_continuous(labels = scales::percent_format()) +
	    theme(text=element_text(size=30), axis.text.y = element_text(face = 'bold')) +
	    rotate_x_text(45) + coord_flip()

if(horizontal == T){
plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
	    geom_bar(position = "fill", stat = "identity") +
	    xlab("") +
	    ylab("Percentage") +
	    theme_classic() +
	    #scale_fill_manual(values = colors) +
	    scale_y_continuous(labels = scales::percent_format()) +
	    theme(text=element_text(size=30)) +
	    rotate_x_text(45)
}

pdf(paste0(fn, '.pdf'), width = 10, height = 10)
print(plt)
dev.off()

}


transferPlot = function(qry, ref, fn){

# Plot scATACseq and scRNAseq side by side
p1 <- DimPlot(qry, label = TRUE, repel = TRUE, label.size = 7, pt.size = 0.4) +
		ggtitle("Query") + 
		NoLegend() +
		scale_colour_hue(drop = F)

p2 <- DimPlot(ref, label = TRUE, repel = TRUE, label.size = 7, pt.size = 0.4) +
		ggtitle("Reference") +
		NoLegend() +
		scale_colour_hue(drop = F)
plt = p1 + p2

pdf(paste0(fn, '.pdf'), width = 20, height = 10)
print(plt)
dev.off()

}

convertHumanGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}


find_degs = function(seur, tocompare, ctypes, dir = '~', equalize = F, suffix = '', qNorm = F){

for(i in 1:length(ctypes)){

	sub = subset(seur, subset = newannot == ctypes[i])
	#Idents(sub) = sub$Species

	# Prepare data for MAST #
	mat <-  sub@assays$RNA@data
	metause = sub@meta.data

	# Keep genes with 10% expression in at least one species
	cells1 = WhichCells(sub, idents = sps[1])
	cells2 = WhichCells(sub, idents = sps[2])
	cells3 = WhichCells(sub, idents = sps[3])

	# Subset cells to make them equal
	if(equalize == T){
		minSize = min(length(cells1), length(cells2), length(cells3))
		cells1 = cells1[sample(1:length(cells1), minSize)]
		cells2 = cells2[sample(1:length(cells2), minSize)]
		cells3 = cells3[sample(1:length(cells3), minSize)]
	}

	if(qNorm == T){
		rat1 = apply(mat[,cells1], 1, function(x){sum(x > 0)/length(x)})
		rat2 = apply(mat[,cells2], 1, function(x){sum(x > 0)/length(x)})
		rat3 = apply(mat[,cells3], 1, function(x){sum(x > 0)/length(x)})

		rat1Q = ecdf(rat1)(0.1)
		rat2Cut = quantile(rat2, rat1Q)
		rat3Cut = quantile(rat3, rat1Q)

		pass1 = rat1 > 0.1
		pass2 = rat2 > rat2Cut
		pass3 = rat3 > rat3Cut
	} else{

		pass1 = apply(mat[,cells1], 1, function(x){sum(x > 0)/length(x)}) > 0.1
		pass2 = apply(mat[,cells2], 1, function(x){sum(x > 0)/length(x)}) > 0.1
		pass3 = apply(mat[,cells3], 1, function(x){sum(x > 0)/length(x)}) > 0.1
	}

	# Pairwise species comparison matrix.
	# Note that genes surviving only in the other species will also be tested.
	testcells1 = WhichCells(sub, idents = tocompare[1])
	testcells2 = WhichCells(sub, idents = tocompare[2])

	mat = mat[pass1 | pass2 | pass3, c(testcells1, testcells2)]

	print(paste0('Gene#: ', nrow(mat), ' ', sum(pass1), ' ', sum(pass2), ' ', sum(pass3)))
	print(paste0('Human#: ', length(cells1)))
	print(paste0('Chimp#: ', length(cells2)))
	print(paste0('Macaque#: ', length(cells3)))

	metause = metause[colnames(mat),]
	metause$Species = droplevels(metause$Species)

	# Create MAST object
	sca <- MAST::FromMatrix(exprsArray = as.matrix(x = mat),
				cData = metause,
				fData = data.frame(rownames(mat)))

	# Scale number of detected genes (cngeneson. Same as nFeature_RNA only after filtering)
	cdr2 <- colSums(assay(sca)>0)
	colData(sca)$cngeneson <- scale(cdr2)

	# Test with MAST
	mastfix = MAST::zlm(~Species + cngeneson + human_age + lib_batch + sex, sca)
	summaryCond <- summary(object = mastfix, doLRT = paste0('Species',tocompare[2]))
	summaryDt <- summaryCond$datatable
	p_val <- summaryDt[summaryDt$component == "H", 4]
	genes.return <- summaryDt[summaryDt$component == "H", 1]
	to.return <- data.frame(p_val, row.names = genes.return$primerid)
	fix_res = to.return

	# Calculate log fold change
	avg_logfc <- log(rowMeans(expm1(mat[,testcells1])) + 1) - log(rowMeans(expm1(mat[,testcells2])) + 1)
	fix_res$avg_logfc = avg_logfc[rownames(fix_res)]
	colnames(fix_res)[1] = 'p_value'

	# Calculate FDR
	fix_res$adj_p_value = p.adjust(fix_res$p_value, method = 'BH')

	saveRDS(fix_res, paste0(dir, '/', paste0(tocompare[1], '_', tocompare[2], '_', ctypes[i], '_DEGs', suffix, '.RDS')))

	print(paste0('Completed cluster ', ctypes[i]))
}
}




