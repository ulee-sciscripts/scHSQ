\name{calcANOVACor}
\alias{calcANOVACor}

\title{Calculates degree of cell type-specific evolutionary divergence}

\description{
Calculates degree of cell type-specific evolutionary divergence. This is done by calculating the Spearman correlation between log(F-statistic) value for species-of-origin labels and cell type assignment labels.
}

\usage{
calcANOVACor(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A numerical value computing the degree of cell type-specific evolutionary divergence.}

\author{UnJin Lee}

\examples{\dontrun{
## calculate degree of correlation between species-of-origin log(F-statistic) and cell type assignment log(F-statistic) across all genes.
anova_cor <- calcANOVACor(in_scHSQ)
}}