require(Biobase)
require(monocle3)
require(ggplot2)
	
#define scHSQ class containing all information for HSQ application
setClass("scHSQ", 
	representation(fullCDS="cell_data_set", n_species="numeric",
		n_celltype="numeric", ref_spp="character", spp_labs="character", idx_ref="numeric",
		celltype_labs="character", col_scheme="character", p_species="numeric",
		p_celltype="numeric", anova_stats="matrix", marker_list="character",
		reducedCDS="cell_data_set", hsq_status="logical"))

#method definitions
setGeneric("new_scHSQ", signature=c("in_CDS", "ref_species", "mdat_lab_spec", "mdat_lab_ct"), 
		function(in_CDS, ref_species, mdat_lab_spec, mdat_lab_ct, ...)
		standardGeneric("new_scHSQ"))

setGeneric("genMarkGenes", signature = c("in_scHSQ"), function(in_scHSQ, ...) standardGeneric("genMarkGenes"))

setGeneric("getMarkGenes", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getMarkGenes"))

setGeneric("getFullCDS", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getFullCDS"))

setGeneric("getReducedCDS", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getReducedCDS"))

setGeneric("setReducedCDS", signature = c("in_scHSQ", "in_CDS"), function(in_scHSQ, in_CDS) standardGeneric("setReducedCDS"))

setGeneric("getFullPData", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getFullPData"))

setGeneric("getFullExprs", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getFullExprs"))

setGeneric("getRefIdx", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getRefIdx"))

setGeneric("setRefIdx", signature = c("in_scHSQ", "idx"), function(in_scHSQ, idx) standardGeneric("setRefIdx"))

setGeneric("getRefSpecExprs", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getRefSpecExprs"))

setGeneric("getSppLabs", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getSppLabs"))

setGeneric("getCTLabs", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getCTLabs"))

setGeneric("setCellType", signature = c("in_scHSQ", "ct_labels", "ct_field_new"), function(in_scHSQ, ct_labels, ct_field_new) standardGeneric("setCellType"))

setGeneric("getRefCellType", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getRefCellType"))

setGeneric("getSpp", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getSpp"))

setGeneric("setPThresh", signature = c("in_scHSQ"),
	function(in_scHSQ, ...) standardGeneric("setPThresh"))

setGeneric("findPercentile", signature = c("dat", "p"), function(dat, p) standardGeneric("findPercentile"))

setGeneric("setANOVAStats", signature = c("in_scHSQ", "in_anova"), function(in_scHSQ, in_anova) standardGeneric("setANOVAStats"))

setGeneric("getANOVAStats", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("getANOVAStats"))

setGeneric("calcANOVA", signature = c("in_scHSQ"), function(in_scHSQ, ...) standardGeneric("calcANOVA"))

setGeneric("remapCTAssign", signature = c("in_scHSQ", "ct_map"), function(in_scHSQ, ct_map) standardGeneric("remapCTAssign"))

setGeneric("calcANOVACor", signature = c("in_scHSQ"), function(in_scHSQ) standardGeneric("calcANOVACor"))

setGeneric("plotANOVA", signature = c("in_scHSQ"), function(in_scHSQ, ...) standardGeneric("plotANOVA"))

setGeneric("applyMarkGenes", signature = c("in_scHSQ"), function(in_scHSQ, ...) standardGeneric("applyMarkGenes"))

setGeneric("reclusterCDS", signature = c("in_scHSQ"), function(in_scHSQ, ...) standardGeneric("reclusterCDS"))

setGeneric("plotCTAssign", signature = c("in_scHSQ"), function(in_scHSQ, ...) standardGeneric("plotCTAssign"))

setGeneric("plotGeneExprs", signature = c("in_scHSQ", "gene"), function(in_scHSQ, gene, ...) standardGeneric("plotGeneExprs"))


#initializing new scHSQ object
setMethod("new_scHSQ", c(in_CDS="cell_data_set", ref_species="character", mdat_lab_spec="character",
			mdat_lab_ct="character"), 
	function(in_CDS, ref_species, mdat_lab_spec, mdat_lab_ct, color_scheme="default", p_spec=0.5, p_ct=0.5){
		this_scHSQ_obj <- new(Class="scHSQ")
		this_scHSQ_obj@fullCDS <- in_CDS
		this_scHSQ_obj@anova_stats <- matrix(0, nrow=2, ncol=2)
		this_scHSQ_obj@marker_list <- as.character(c(""))
		this_scHSQ_obj@reducedCDS <- new(Class="cell_data_set")
		this_scHSQ_obj@hsq_status <- FALSE

		ref_spec_idx <- NA
		
		out_err <- "\nNon-valid scHSQ object. Aborting."

		#checks if mdat_lab_spec is valid (proper label in pData)
		tryCatch(
			expr = {
				if(! mdat_lab_spec %in% colnames(monocle3::pData(in_CDS))) stop("mdat_lab_spec not found in pData(in_CDS)")
				this_scHSQ_obj@spp_labs <- mdat_lab_spec
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)

		#checks if mdat_lab_ct is valid (proper label in pData)
		tryCatch(
			expr = {
				if(! mdat_lab_ct %in% colnames(monocle3::pData(in_CDS))) stop("mdat_lab_ct not found in pData(in_CDS)")
				this_scHSQ_obj@celltype_labs <- mdat_lab_ct
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)		
		
		#checks if p_spec is valid (is on [0, 1])
		tryCatch(
			expr = {
				if(! ((p_spec >= 0) && (p_spec <= 1))) stop("p_spec not valid")
				this_scHSQ_obj@p_species <- p_spec
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)
		
		#checks if p_ct is valid (is on [0,1])
		tryCatch(
			expr = {
				if(! ((p_ct >= 0) && (p_ct <= 1))) stop("p_ct not valid")
				this_scHSQ_obj@p_celltype <- p_ct
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)

		#checks if ref_species is in mdat_lab_spec
		tryCatch(
			expr = {
				if(!ref_species %in% as.character(monocle3::pData(in_CDS)[,this_scHSQ_obj@spp_labs])) stop("ref_species not in pData(in_CDS)[, spp_labs]")
				
				ref_spec_idx <- which(as.character(monocle3::pData(in_CDS)[,this_scHSQ_obj@spp_labs]) == ref_species)
				# message(paste("number of reference species cells: ", length(ref_spec_idx), sep=""))
				
				this_scHSQ_obj@ref_spp <- ref_species
				this_scHSQ_obj@idx_ref <- ref_spec_idx

			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)

		#checks if reference species has cell type labels
		tryCatch(
			expr = {
				if(length(levels(as.factor(as.character(monocle3::pData(in_CDS)[this_scHSQ_obj@idx_ref, mdat_lab_ct])))) <= 1) stop("Only one annotated cell type in reference species")
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)
		
		#calculate number of species
		tryCatch(
			expr = {
				this_scHSQ_obj@n_species <- length(levels(as.factor(as.character(monocle3::pData(in_CDS)[,mdat_lab_spec]))))
				if(this_scHSQ_obj@n_species <= 1) stop("Only one annotated species in data set")
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)
		
		#calculate number of cell types
		tryCatch(
			expr = {
				this_ctlabs <- data.frame(idx=this_scHSQ_obj@idx_ref, ct=monocle3::pData(in_CDS)[this_scHSQ_obj@idx_ref, as.character(this_scHSQ_obj@celltype_labs)])
				
				# remove_idx <- which(this_ctlabs[,2] %in% c(NA, "unknown", "unassigned"))
				# this_ctlabs <- this_ctlabs[-remove_idx,]
				
				this_scHSQ_obj@idx_ref <- this_ctlabs[,1]
				
				this_scHSQ_obj@n_celltype <- length(levels(as.factor(as.character(this_ctlabs[,2]))))
				if(this_scHSQ_obj@n_celltype<= 1) stop("Only one known cell type in reference species")
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)			
		
		#check if color_scheme is valid and generate color_scheme if "default"
		tryCatch(
			expr = {
				if(color_scheme == "default"){
					gg_color_hue <- function(n) {
						hues = seq(15, 375, length = n + 1)
						hcl(h = hues, l = 65, c = 100)[1:n]
					}

					color_scheme <- gg_color_hue(this_scHSQ_obj@n_celltype)
				}
				else{
					areColors <- function(x) {
						sapply(x, function(X) {
						tryCatch(is.matrix(col2rgb(X)), 
							error = function(e) FALSE)
						})
					}
				
					if(!sum(areColors(color_scheme)) == this_scHSQ_obj@n_celltype) stop("color_scheme is not valid")
				}
				
				this_scHSQ_obj@col_scheme <- color_scheme
			},
			error = function(e){
				message(e)
				message(out_err)
				return(NA)
			}
		)
		
		return(this_scHSQ_obj)
	})



#genMarkGenes
setMethod("genMarkGenes", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, run_ANOVA=TRUE){
		if(!in_scHSQ@hsq_status){
			message("\nANOVA scores not available.")
			if(!run_ANOVA){
				message("\nNothing to do. Returning input.")
				return(in_scHSQ)
			}
			else{
				message("\nCalculating ANOVA scores.")
				in_scHSQ <- calcANOVA(in_scHSQ)
			}
		}
		anova_stats <- getANOVAStats(in_scHSQ)
		
		genes_spec <- rownames(in_scHSQ@anova_stats)[which(in_scHSQ@anova_stats[,getSppLabs(in_scHSQ)] <= findPercentile(in_scHSQ@anova_stats[,getSppLabs(in_scHSQ)], in_scHSQ@p_species))]
		genes_ct <- rownames(in_scHSQ@anova_stats)[which(in_scHSQ@anova_stats[,getCTLabs(in_scHSQ)] >= findPercentile(in_scHSQ@anova_stats[,getCTLabs(in_scHSQ)], in_scHSQ@p_celltype))]
		
		in_scHSQ@marker_list <- intersect(genes_spec, genes_ct)
		return(in_scHSQ)
})

##setter and getter functions

#getMarkGenes	
setMethod("getMarkGenes", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		if(!in_scHSQ@hsq_status){
			warning("ANOVA scores not available.")
			warning("Nothing to do. Returning empty list.")
			return(c(""))
		}
		
		return(in_scHSQ@marker_list)
})

#getXCDS functions
setMethod("getFullCDS", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(in_scHSQ@fullCDS)
})
	
setMethod("getReducedCDS", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
	if(!in_scHSQ@hsq_status){
		warning("Reduced CDS not available.")
		warning("Nothing to do. Returning empty CDS.")
		return(new(Class="cell_data_set"))
	}
	
	return(in_scHSQ@reducedCDS)
})

setMethod("setReducedCDS", c(in_scHSQ="scHSQ", in_CDS="cell_data_set"),
	function(in_scHSQ, in_CDS){
		in_scHSQ@reducedCDS <- in_CDS
		return(in_scHSQ)
})

#getFullPData	
setMethod("getFullPData", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(monocle3::pData(in_scHSQ@fullCDS))
})

#getFullExprs	
setMethod("getFullExprs", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(monocle3::exprs(in_scHSQ@fullCDS))
})

#getRefIdx	
setMethod("getRefIdx", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(in_scHSQ@idx_ref)
})

#setRefIdx
setMethod("setRefIdx", c(in_scHSQ="scHSQ", idx="numeric"),
	function(in_scHSQ, idx){
		in_scHSQ@idx_ref <- idx
		return(in_scHSQ)
})

#getRefSpecExprs
setMethod("getRefSpecExprs", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(monocle3::exprs(getFullCDS(in_scHSQ))[, in_scHSQ@idx_ref])
})

#getSppLabs	
setMethod("getSppLabs", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(in_scHSQ@spp_labs)
})

#getCTLabs	
setMethod("getCTLabs", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(in_scHSQ@celltype_labs)
})

#setCellType
setMethod("setCellType", c(in_scHSQ="scHSQ", ct_labels="character", ct_field_new="character"),
	#assuming names of ct_labels correspond to UMIs
	function(in_scHSQ, ct_labels, ct_field_new){
		new_CDS <- getReducedCDS(in_scHSQ)
		
		new_pData <- cbind(monocle3::pData(new_CDS), ct_labels)
		names(new_pData)[ncol(new_pData)] <- ct_field_new
		
		monocle3::pData(new_CDS) <- new_pData
		
		in_scHSQ <- setReducedCDS(in_scHSQ, new_CDS) 

		return(in_scHSQ)
})

#getRefCellType	
setMethod("getRefCellType", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(monocle3::pData(getFullCDS(in_scHSQ))[getRefIdx(in_scHSQ), getCTLabs(in_scHSQ)])
})

#getSpp	
setMethod("getSpp", c(in_scHSQ="scHSQ"),
	function(in_scHSQ){
		return(monocle3::pData(getFullCDS(in_scHSQ))[, getSppLabs(in_scHSQ)])
})

#setPThresh	
setMethod("setPThresh", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, p_spec=0.5, p_ct=0.5, gen_list=TRUE){
		in_scHSQ@p_species <- p_species
		in_scHSQ@p_celltype <- p_ct
		
		if(in_scHSQhsq_status && gen_list){
			message("\nRe-generating marker list with new thresholds")
			in_scHSQ <- genMarkGenes(in_scHSQ)			
		}
		
		return(in_scHSQ)
})

#findPercentile
setMethod("findPercentile", c(dat="numeric", p="numeric"),
		function(dat, p){
			dat <- sort(dat)
			e_cdf <- 1:length(dat) / length(dat)
			
			return(dat[which(e_cdf >= p)[1]])
})

#setANOVAStats	
setMethod("setANOVAStats", c(in_scHSQ="scHSQ", in_anova="matrix"),
	function(in_scHSQ, in_anova){
		in_scHSQ@anova_stats <- in_anova
		in_scHSQ@hsq_status <- TRUE
		
		return(in_scHSQ)
})

#getANOVAStats
setMethod("getANOVAStats", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ){
		if(!in_scHSQ@hsq_status){
			warning("ANOVA scores not available. Try running calcANOVA() first.")
			warning("Nothing to do. Returning blank matrix.")
		}
	
		return(in_scHSQ@anova_stats)
})

#calcANOVA	
setMethod("calcANOVA", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, verbose=TRUE){
		cds_in <- getFullCDS(in_scHSQ)
		
		message("Setting up reference for cell type labels")
		expr_matrix_ingroup <- getRefSpecExprs(in_scHSQ)
		expr_matrix_ingroup_labs <- getRefCellType(in_scHSQ)
		
		message("Setting up reference for species...")
		expr_matrix <- getFullExprs(in_scHSQ)
		expr_matrix_labs <- getSpp(in_scHSQ)
		
		anova_statistics <- matrix(0, nrow=nrow(getFullExprs(in_scHSQ)), ncol=2)
		rownames(anova_statistics) <- rownames(getFullExprs(in_scHSQ))
		colnames(anova_statistics) <- c(getSppLabs(in_scHSQ), getCTLabs(in_scHSQ))
		
		for(i in 1:nrow(anova_statistics)){
			if(verbose) message(paste(c("(", i, ")/(", nrow(anova_statistics),  ") Calculating ANOVA stats for ", rownames(anova_statistics)[i], sep="")))
			
			anova_statistics[i, 1] <- oneway.test(as.numeric(expr_matrix[i,]) ~ as.character(expr_matrix_labs), var.equal=F)$statistic
			anova_statistics[i, 2] <- oneway.test(as.numeric(expr_matrix_ingroup[i,]) ~ as.character(expr_matrix_ingroup_labs), var.equal=F)$statistic
		}
		
		message("\n\nCalculations complete")
		
		in_scHSQ <- setANOVAStats(in_scHSQ, anova_statistics)
		
		return(in_scHSQ)
})

#remapCTAssign
setMethod("remapCTAssign", c(in_scHSQ="scHSQ", ct_map="data.frame"), 
#ct_map is assumed to have first column as list of original cell type labels, second column as new labels
#colnames of ct_map must correlate to labels in pData
	function(in_scHSQ, ct_map){
		ct_assign_old <- monocle3::pData(getReducedCDS(in_scHSQ))[,colnames(ct_map)[1]]
		ct_assign_new <- ct_assign_old

		if(colnames(ct_map)[2] %in% colnames(monocle3::pData(getReducedCDS(in_scHSQ)))){
			message("Name for reassigned labels already present in input scHSQ object")
			message("Returning input")
			return(in_scHSQ)
		}

		for(i in 1:nrow(ct_map)){
			ct_assign_new[which(ct_assign_old == ct_map[i, 1])] <- ct_map[i, 2] 
		}
		
		names(ct_assign_new) <- rownames(ct_map)
		
		return(setCellType(in_scHSQ, ct_assign_new, colnames(ct_map)[2]))
})

#calcANOVACor
setMethod("calcANOVACor", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ){
		anova_stat <- getANOVAStats(in_scHSQ)
		
		idx <- intersect(which(!is.na(anova_stat[,1])), which(!is.na(anova_stat[,2])))
		anova_cor <- cor(log(anova_stat[idx,1]), log(anova_stat[idx,2]))
		anova_cor <- round(anova_cor, digits=3)

		return(anova_cor)
})

#plotANOVA
setMethod("plotANOVA", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, verbose=F, cex=2){
		if(verbose) message("Recommend minimum settings for raster graphics: width=1500, height=1500, res=300, pointsize=5")
	
		anova_stat <- getANOVAStats(in_scHSQ)
		idx <- intersect(which(!is.na(anova_stat[,1])), which(!is.na(anova_stat[,2])))
		anova_cor <- calcANOVACor(in_scHSQ)
		
		
		anova_main <- bquote("ANOVA log(F-Statistic) by Gene, " ~ rho == .(anova_cor))
		plot(log(anova_stat[idx,2]), log(anova_stat[idx,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75*cex/2, cex.main=cex)
		mtext("Predictive of Cell Type", side=1, line=3.2*cex/2, cex=cex)
		mtext("Predictive of Species", side=2, line=2.6*cex/2, cex=cex)
		abline(v=log(findPercentile(getANOVAStats(in_scHSQ)[,getCTLabs(in_scHSQ)], in_scHSQ@p_celltype)))
		abline(h=log(findPercentile(getANOVAStats(in_scHSQ)[,getSppLabs(in_scHSQ)], in_scHSQ@p_species)))
		
		return()
})

#applyMarkGenes
setMethod("applyMarkGenes", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, leiden_res=3e-4, new_clu="HSQ_clu"){
		if(!in_scHSQ@hsq_status){
			warning("ANOVA scores not available. Try running calcANOVA() first.")
			warning("Nothing to do. Returning unprocessed input.")
			return(in_scHSQ)
		}
		
		rd_share <- DataFrame(row.names=rownames(getFullExprs(in_scHSQ)))
		rd_share$id <- rownames(rd_share)
		rd_share$gene_short_name <- rownames(rd_share)	
		
		cds_comb_diffgene <- new_cell_data_set(getFullExprs(in_scHSQ)[getMarkGenes(in_scHSQ),], cell_metadata=colData(getFullCDS(in_scHSQ)), gene_metadata=rd_share[getMarkGenes(in_scHSQ),])
		
		if(new_clu %in% colnames(monocle3::pData(cds_comb_diffgene))){
			warning("Invalid new_clu label: already exists in pData")
			warning("Nothing to do. Returning unprocessed input.")
			return(in_scHSQ)
		}
		
		cds_comb_diffgene <- preprocess_cds(cds_comb_diffgene)
		cds_comb_diffgene <- align_cds(cds_comb_diffgene, alignment_group=getSppLabs(in_scHSQ))
		cds_comb_diffgene <- reduce_dimension(cds_comb_diffgene)
		cds_comb_diffgene <- cluster_cells(cds_comb_diffgene, resolution=leiden_res)
		
		clu_umap <- cds_comb_diffgene@clusters$UMAP$clusters
		monocle3::pData(cds_comb_diffgene)[,new_clu] <- as.character(clu_umap)
		colData(cds_comb_diffgene)[,new_clu] <- as.character(clu_umap)
		
		message(paste("New cluster ID: ", new_clu, sep=""))

		in_scHSQ <- setReducedCDS(in_scHSQ, cds_comb_diffgene)
		
		return(in_scHSQ)
})

#reclusterCDS
setMethod("reclusterCDS", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, leiden_res=3e-4, new_clu="HSQ_clu1"){
		cds_comb_diffgene <- getReducedCDS(in_scHSQ)
		cds_comb_diffgene <- cluster_cells(cds_comb_diffgene, resolution=leiden_res)
		
		clu_umap <- cds_comb_diffgene@clusters$UMAP$clusters
		
		if(new_clu %in% colnames(monocle3::pData(cds_comb_diffgene))){
			warning("Invalid new_clu label: already exists in pData")
			warning("Nothing to do. Returning unprocessed input.")
			return(in_scHSQ)
		}
		
		monocle3::pData(cds_comb_diffgene)[,new_clu] <- as.character(clu_umap)
		colData(cds_comb_diffgene)[,new_clu] <- as.character(clu_umap)
		
		message(paste("New cluster ID: ", new_clu, sep=""))
		
		in_scHSQ <- setReducedCDS(in_scHSQ, cds_comb_diffgene)
		
		return(in_scHSQ)
})

#plotCTAssign
setMethod("plotCTAssign", c(in_scHSQ="scHSQ"), 
	function(in_scHSQ, new_labels="HSQ_clu", group_labs=TRUE, verbose=FALSE, ...){
		if(verbose) message("Recommend minimum settings for raster graphics: width=1500, height=1500, res=300, pointsize=5")
		
		lab_size <- 3
		if(!group_labs) lab_size <- 0
		
		all_specs <- levels(as.factor(getFullPData(in_scHSQ)[, getSppLabs(in_scHSQ)]))
		this_spp <- getFullPData(in_scHSQ)[,getSppLabs(in_scHSQ)]
		
		plots_list <- list()
		
		cds_ref <- getReducedCDS(in_scHSQ)[,getRefIdx(in_scHSQ)]

		plots_list[[1]] <-  plot_cells(getReducedCDS(in_scHSQ)[, getRefIdx(in_scHSQ)], color_cells_by=getCTLabs(in_scHSQ), group_cells_by=getCTLabs(in_scHSQ), group_label_size = lab_size, ...) + ggtitle("Known Annotation in Reference Species") #+ scale_color_manual(in_scHSQ@col_scheme)
		
		plots_list[[2]] <-  plot_cells(getReducedCDS(in_scHSQ), color_cells_by=new_labels, group_cells_by=new_labels, group_label_size = lab_size, ...) + ggtitle("HSQ Clusters in All Species") #+ scale_color_manual(in_scHSQ@col_scheme)
		
		for(i in 1:length(all_specs)){
			plots_list[[i+2]] <- plot_cells(getReducedCDS(in_scHSQ)[,this_spp == all_specs[i]], color_cells_by=new_labels, group_cells_by=new_labels, group_label_size = lab_size, ...) + ggtitle(paste("HSQ Clusters in ", all_specs[i], sep="")) #+ scale_color_manual(in_scHSQ@col_scheme)
		}
	
		return(plots_list)
})

#plotGeneExprs
setMethod("plotGeneExprs", c(in_scHSQ="scHSQ", gene="character"), 
	function(in_scHSQ, gene, species=NULL, label_cell_groups=FALSE, alpha=0.35, ...){
		all_specs <- levels(as.factor(getFullPData(in_scHSQ)[, getSppLabs(in_scHSQ)]))
		this_spp <- getFullPData(in_scHSQ)[,getSppLabs(in_scHSQ)]
		
		this_idx <- 1:length(this_spp)
		if(!is.null(species)) this_idx <- which(this_spp == species)

		this_redCDS <- getReducedCDS(in_scHSQ)[, this_idx]
		this_redExprs <- monocle3::exprs(this_redCDS)
		this_fullExprs <- getFullExprs(in_scHSQ)[, this_idx]
		
		if(gene %in% getMarkGenes(in_scHSQ)){
			this_p <- plot_cells(this_redCDS, genes=gene)
			return(this_p)
		}
		
		this_redExprs[1,] <- this_fullExprs[gene,]
		rownames(this_redExprs)[1] <- gene
		
		this_metadata <- DataFrame(row.names=rownames(this_redExprs))
		this_metadata$id <- rownames(this_metadata)
		this_metadata$gene_short_name <- rownames(this_metadata)

		this_newCDS <- new_cell_data_set(this_redExprs, cell_metadata=colData(this_redCDS), gene_metadata=this_metadata)
		
		this_redCDS@assays <- this_newCDS@assays
		this_redCDS@elementMetadata <- this_newCDS@elementMetadata
		this_redCDS@int_elementMetadata <- this_newCDS@int_elementMetadata
		this_redCDS@rowRanges <- this_newCDS@rowRanges
		
		this_p <- plot_cells(this_redCDS, genes=gene, label_cell_groups=label_cell_groups, alpha=alpha, ...)
		return(this_p)	
})
