\name{genMarkGenes}
\alias{genMarkGenes}

\title{Generate Marker Genes}

\description{
Applies HSQ criteria to generate a list of marker genes with high cell type specificity and low species specificity.
}

\usage{
genMarkGenes(in_scHSQ, run_ANOVA=TRUE)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{run_ANOVA}{
		a logical for calculating new ANOVA values if they are not present. Defaults to TRUE.
		}
}

\value{A scHSQ object that may be used for further processing.}

\author{UnJin Lee}

\examples{\dontrun{
## Generate new marker gene set
dros_testis <- genMarkGenes(dros_testis)
}}