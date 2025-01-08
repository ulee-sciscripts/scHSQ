\name{plotGeneExprs}
\alias{plotGeneExprs}

\title{Plot gene expression in cross-species UMAP}

\description{
Plot gene expression in cross-species UMAP. Returned object is a UMAP projection showing expression for a single gene. Note that arguments for monocle3's plot_cells function can be passed to plotCTAssign function and will automatically be applied to all elements of output. Rasterized figures are recommended to use width=1500, height=1500, res=300, pointsize=5.
}

\usage{
plotGeneExprs(in_scHSQ, gene, species=NULL, label_cell_groups=FALSE, alpha=0.35, ...)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{gene}{
		a character specifying which gene to plot
		}
	\item{species}{
		a character specifying which species to plot. Defaults to NULL for all species.
		}
	\item{label_cell_groups}{
		a logical parameter passed to plot_cells. Defaults to FALSE.
		}
	\item{alpha}{
		a numerical value passed to plot_cells. Defaults to 0.35.
		}
}

\value{A ggplot2 object.}

\author{UnJin Lee}

\examples{\dontrun{
## Generate figures on cross-species UMAP

gene_exprs <- plotGeneExprs(dros_testis, "Rbp4", alpha=0.5, norm_method="size_only", min_expr=40)
	
png("scHSQ_UMAP_Rbp4.png", width=1500, height=1500, res=300, pointsize=5)
gene_exprs
dev.off()

}}