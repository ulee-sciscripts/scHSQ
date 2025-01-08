\name{reclusterCDS}
\alias{reclusterCDS}

\title{Re-cluster over processed cell_data_set}

\description{
Re-cluster over processed cell_data_set. Useful for adjusting resolution of clustering to find better cell type mappings.
}

\usage{
reclusterCDS(in_scHSQ, leiden_res=3e-4, new_clu="HSQ_clu1")
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{leiden_res}{
		a numerical value corresponding to the resolution for leiden clustering. Arbitrarily defaults to 3e-4 for historical reasons.
		}
	\item{new_clu}{
		a character that will be used to label the new cluster. Defaults to "HSQ_clu1".
		}
}

\value{A scHSQ object with newly processed cell_data_set containing cross-species UMAP.}

\author{UnJin Lee}

\examples{\dontrun{
## Apply newly identified marker genes passing HSQ criteria

dros_testis <- genMarkGenes(dros_testis)
dros_testis <- applyMarkGenes(dros_testis)
dros_testis <- reclusterCDS(dros_testis, leiden_res=1e-4)
}}