

da_peaks = function(clustermat, clustermeta, sf, covars){

	# This function fits two models: with and without 'Species' covariate.
	# Then uses likelihood ratio test to compute the p-value.

	# Get cells
	cells1 = clustermeta[clustermeta$Species == names(table(clustermeta$Species))[1], "barcode"]
	cells2 = clustermeta[clustermeta$Species == names(table(clustermeta$Species))[2], "barcode"]

	print(paste0("Testing: ", nrow(clustermat), " peaks"))

	# Check sample size before running
	print(table(clustermeta$Species))
	print(table(clustermeta$Sample))

	# Fit two models
	glmfit_h0 = list()
	glmfit_h1 = list()
	lr_res = list()
	acc1s = list()
	acc2s = list()
	deltas = list()

	fcs = list()
	for(i in 1:nrow(clustermat)){

		reg = rownames(clustermat)[i]

		acc1s[[i]] = mean(clustermat[i, cells1])
		acc2s[[i]] = mean(clustermat[i, cells2])

		# Calculate delta and fold change. Always normalize second group to first group.
		deltas[[i]] = mean(clustermat[i, cells1]) - (mean(clustermat[i, cells2]) * sf)
		fcs[[i]] = mean(clustermat[i, cells1]) / (mean(clustermat[i, cells2]) * sf)

		# Only total accessibility as covariate
		glmfit_h0[[i]] = glm(as.formula(paste("clustermat[reg, ] ~ ", paste(covars[2:length(covars)], collapse = '+'))), data = clustermeta, family="binomial")

		# Only species as binomial predictor in alternative model
		glmfit_h1[[i]] = glm(as.formula(paste("clustermat[reg, ] ~ ", paste(covars, collapse = '+'))), data = clustermeta, family="binomial")

		# Does adding the species result in significantly better fit?
		# Note that h0 will always give equal or worse likelihood than h1, therefore we can use lr test.
		lr_res[[i]] = lmtest::lrtest(glmfit_h0[[i]], glmfit_h1[[i]])

		if((i %% 1000) == 0){print(i)}
	}

	# Create results matrix
	pval = lapply(lr_res, function(x){x[2,5]}) # Likelihood test p-value

	# Create results. Convert fold change to log2 scale for easier interpretation. Adjust p-value with BH.
	res = data.frame(pval = unlist(pval), fc = unlist(fcs),
			 acc1 = unlist(acc1s), acc2 = unlist(acc2s), delta = unlist(deltas))

	#res = data.frame(pval = unlist(pval), fc = tmp)
	res$adj_pval = p.adjust(res$pval, method = 'BH', n = length(res$pval))
	rownames(res) = rownames(clustermat)[1:nrow(res)]

	res$condition = ifelse(res$delta > 0, "UP", "DOWN") # >0 is up in human

	res_sign = res[res$adj_pval < 0.05 &
			((res$delta > sd(res$delta)*1.5) |
			(res$delta < sd(res$delta)*-1.5)),]



	## Plot the results ##

	# Which regions are open in one or other?
	r1 = rowSums(clustermat[,cells1]) / length(cells1)
	r2 = rowSums(clustermat[,cells2]) / length(cells2)

	r1 = r1[names(r1) %in% rownames(clustermat)]
	r2 = r2[names(r2) %in% rownames(clustermat)]

	tmpdf = data.frame(peaks = names(r1), Species_1 = r1, Species_2 = r2, nums = 1:length(r1))

	tmpdf$is_sign = ifelse(tmpdf$peaks %in% rownames(res_sign), "1", "0")
	tmpdf$is_sign = factor(tmpdf$is_sign)
	
	plt = ggscatter(tmpdf %>% arrange(is_sign), x = 'Species_1', y = 'Species_2', color = 'is_sign',
		alpha = 0.5) + scale_colour_manual(values = c("grey", "red")) +
		xlab(names(table(clustermeta$Species))[1]) +
		ylab(names(table(clustermeta$Species))[2])

	pdf(paste(paste(names(table(clustermeta$Species)), collapse = "_"),
	    clustermeta$newannot[1], "scatter_differential.pdf", sep = '_'))
	print(plt)
	dev.off()

	rm(glmfit_h1, glmfit_h0); gc()

	res

}


# Same as above, but can be parallelized easier
da_peaksPar = function(allcomb){

	# This function fits two models: with and without 'Species' covariate.
	# Then uses likelihood ratio test to compute the p-value.

	# Get variables from the list
	require(Matrix)
	clustermat = as(allcomb[[1]], 'dgCMatrix')

	clustermeta = allcomb[[2]]
	sf = allcomb[[3]]
	covars = allcomb[[4]]

	rm(allcomb)
	gc()

	# Get cells
	cells1 = clustermeta[clustermeta$Species == names(table(clustermeta$Species))[1], "barcode"]
	cells2 = clustermeta[clustermeta$Species == names(table(clustermeta$Species))[2], "barcode"]

	print(paste0("Testing: ", nrow(clustermat), " peaks"))

	# Check sample size before running
	print(table(clustermeta$Species))
	print(table(clustermeta$Sample))

	# Fit two models
	glmfit_h0 = list()
	glmfit_h1 = list()
	lr_res = list()
	acc1s = list()
	acc2s = list()
	deltas = list()

	fcs = list()
	for(i in 1:nrow(clustermat)){

		reg = rownames(clustermat)[i]

		acc1s[[i]] = mean(clustermat[i, cells1])
		acc2s[[i]] = mean(clustermat[i, cells2])

		# Calculate delta and fold change. Always normalize second group to first group.
		deltas[[i]] = mean(clustermat[i, cells1]) - (mean(clustermat[i, cells2]) * sf)
		fcs[[i]] = mean(clustermat[i, cells1]) / (mean(clustermat[i, cells2]) * sf)

		# Only total accessibility as covariate
		glmfit_h0 = glm(as.formula(paste("clustermat[reg, ] ~ ", paste(covars[2:length(covars)], collapse = '+'))), data = clustermeta, family="binomial")

		# Only species as binomial predictor in alternative model
		glmfit_h1 = glm(as.formula(paste("clustermat[reg, ] ~ ", paste(covars, collapse = '+'))), data = clustermeta, family="binomial")

		# Does adding the species result in significantly better fit?
		# Note that h0 will always give equal or worse likelihood than h1, therefore we can use lr test.
		lr_res[[i]] = lmtest::lrtest(glmfit_h0, glmfit_h1)

		if((i %% 10) == 0){print(i)}
	}

	# Create results matrix
	pval = lapply(lr_res, function(x){x[2,5]}) # Likelihood test p-value

	# Create results. Convert fold change to log2 scale for easier interpretation. Adjust p-value with BH.
	res = data.frame(pval = unlist(pval), fc = unlist(fcs),
			 acc1 = unlist(acc1s), acc2 = unlist(acc2s), delta = unlist(deltas))

	#res = data.frame(pval = unlist(pval), fc = tmp)
	res$adj_pval = p.adjust(res$pval, method = 'BH', n = length(res$pval))
	rownames(res) = rownames(clustermat)[1:nrow(res)]

	res$condition = ifelse(res$delta > 0, "UP", "DOWN") # >0 is up in human

	rm(glmfit_h1, glmfit_h0); gc()

	res

}

