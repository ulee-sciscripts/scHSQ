\name{getSpp}
\alias{getSpp}

\title{Queries metadata for species of origin}

\description{
Queries metadata species of origin for each cell. This is the species information contained in the pData of contained cell_data_set objects.
}

\usage{
getSpp(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A character vector corresponding to species of origin for each cell.}

\author{UnJin Lee}

\examples{\dontrun{
## Query species of origin metadata
spp_info <- getSpp(dros_testis)
}}