\name{getRefCellType}
\alias{getRefCellType}

\title{Queries cell type information for reference species}

\description{
Queries cell type information for reference species.
}

\usage{
getRefCellType(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A character vector corresponding to the cell type information for the reference species.}

\author{UnJin Lee}

\examples{\dontrun{
## Query reference cell type information
mel_celltypes <- getRefCellType(dros_testis)
}}