\name{getMarkGenes}
\alias{getMarkGenes}

\title{Query Marker Genes}

\description{
Queries marker gene set.
}

\usage{
getMarkGenes(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A character vector of marker genes.}

\author{UnJin Lee}

\examples{\dontrun{
## Query marker gene set
marker_genes <- getMarkGenes(dros_testis)
}}