\name{applyMarkGenes}
\alias{applyMarkGenes}

\title{Apply marker genes to generate cross-species UMAP and associated clusterings}

\description{
Apply marker genes to generate cross-species UMAP and associated clusterings. Marker genes are applied by subsetting read count data to marker gene set contained in in_scHSQ object, then follows standard monocle3 processing pipeline. Specifically, process_cds, align_cds (aligning over original species labels), reduce_dimension, and cluster_cells are executed. Resulting processed cell_data_set object is stored in in_scHSQ and returned.
}

\usage{
applyMarkGenes(in_scHSQ, leiden_res=3e-4, new_clu="HSQ_clu")
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{leiden_res}{
		a numerical value corresponding to the resolution for leiden clustering. Arbitrarily defaults to 3e-4 for historical reasons.
		}
	\item{new_clu}{
		a character that will be used to label the new cluster. Defaults to "HSQ_clu".
		}
}

\value{A scHSQ object with newly processed cell_data_set containing cross-species UMAP.}

\author{UnJin Lee}

\examples{\dontrun{
## Apply newly identified marker genes passing HSQ criteria

dros_testis <- genMarkGenes(dros_testis)
dros_testis <- applyMarkGenes(dros_testis)
}}