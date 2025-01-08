\name{getFullCDS}
\alias{getFullCDS}

\title{Query original cell_data_set}

\description{
Queries original cell_data_set.
}

\usage{
getFullCDS(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A cell_data_set that was used for analysis.}

\author{UnJin Lee}

\examples{\dontrun{
## Query original cell_data_set
original_cds <- getFullCDS(dros_testis)
}}