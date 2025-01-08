test_scHSQ <- function(){
	data(testis_3sp)
	checkTrue(is(testis_3sp, "cell_data_set"))
	
	dros_testis <- new_scHSQ(testis_3sp, "mel", "spec", "known_type")
	dros_testis <- genMarkGenes(dros_testis)
	dros_testis <- applyMarkGenes(dros_testis)
	checkTrue(is(dros_testis, "scHSQ"))
	
	new_CDS <- getReducedCDS(dros_testis)
	checkTrue(is(new_CDS, "cell_data_set"))
	
	gene_list <- getMarkGenes(dros_testis)
	checkTrue(is.character(gene_list))
	
	# png("scHSQ_ANOVA.png", width=1500, height=1500, res=300, pointsize=5)
	# plotANOVA(dros_testis)
	# dev.off()
		
	ctplots <- plotCTAssign(dros_testis)
	checkTrue(is(ctplots, "list"))
	
	# png("scHSQ_UMAP_ref.png", width=1500, height=1500, res=300, pointsize=5)
	# ctplots[[1]]
	# dev.off()
	
	# png("scHSQ_UMAP_all.png", width=1500, height=1500, res=300, pointsize=5)
	# ctplots[[2]]
	# dev.off()
	
	old_clu <- as.character(1:10)
	new_clu <- c("Late spermatocytes", "Early spermatocytes", "Early spermatids", "Early spermatocytes", "Early spermatids", "Late spermatogonia", "Somatic", "Late spermatids", "ananassae spermatocytes", "GSC Early spermatogonia")
	
	cell_types <- data.frame(HSQ_clu=old_clu, HSQ_clu1=new_clu)
	
	dros_testis <- remapCTAssign(dros_testis, cell_types)
	
	ctplots_reassign <- plotCTAssign(dros_testis, "HSQ_clu1")
	checkTrue(is(ctplots_reassign, "list"))

	# png("scHSQ_UMAP_reassign.png", width=1500, height=1500, res=300, pointsize=5)
	# ctplots_reassign[[2]]
	# dev.off()

	gene_exprs <- plotGeneExprs(dros_testis, "Rbp4", alpha=0.5, norm_method="size_only", min_expr=40)
	checkTrue(is(gene_exprs, "gg"))
	checkTrue(is(gene_exprs, "ggplot"))
	
	# png("scHSQ_UMAP_Rbp4.png", width=1500, height=1500, res=300, pointsize=5)
	# gene_exprs
	# dev.off()
	
}