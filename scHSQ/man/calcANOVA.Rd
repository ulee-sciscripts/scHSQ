\name{calcANOVA}
\alias{calcANOVA}

\title{Calculates ANOVA statistics for original input cell_data_set}

\description{
Calculates ANOVA statistics for both species and cell type labels to explain variance in read count on a per gene basis.
}

\usage{
calcANOVA(in_scHSQ, verbose)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{verbose}{
		a logical for printing progress of computations. Defaults to TRUE.
		}
}

\value{A scHSQ objected with updated ANOVA F-statistic values for further processing.}

\author{UnJin Lee}

\examples{\dontrun{
## calculate ANOVA F-statistic values
dros_testis <- calcANOVA(dros_testis)
}}