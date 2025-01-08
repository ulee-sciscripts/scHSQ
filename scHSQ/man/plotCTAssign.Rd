\name{plotCTAssign}
\alias{plotCTAssign}

\title{Plot cell type assignments in cross-species UMAP}

\description{
Plot cell type assignments in cross-species UMAP. Returned object is a list of UMAP projections containing: 1) known cell type labels for the reference species, 2) new clusters for all cells across all species, and 3) clusters for cells in each respective species. Note that arguments for monocle3's plot_cells function can be passed to plotCTAssign function and will automatically be applied to all elements of output. Rasterized figures are recommended to use width=1500, height=1500, res=300, pointsize=5.
}

\usage{
plotCTAssign(in_scHSQ, new_labels, group_labs=TRUE, verbose=TRUE, ...)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{new_labels}{
		a character corresponding to a metadata field in the processed cell_data_set's pData slot. Defaults to "HSQ_clu".
		}
	\item{group_labs}{
		a logical for whether or not group labels are included in resulting plots. Defaults to TRUE.
		}
	\item{verbose}{
		a logical for whether or not extra messages are outputted to the console. Defaults to TRUE.
		}
}

\value{A list of ggplot2 objects.}

\author{UnJin Lee}

\examples{\dontrun{
## Generate figures on cross-species UMAP

ctplots <- plotCTAssign(dros_testis)
	
png("scHSQ_UMAP_ref.png", width=1500, height=1500, res=300, pointsize=5)
ctplots[[1]]
dev.off()
	
png("scHSQ_UMAP_all.png", width=1500, height=1500, res=300, pointsize=5)
ctplots[[2]]
dev.off()

}}