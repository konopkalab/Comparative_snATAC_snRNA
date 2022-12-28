
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




