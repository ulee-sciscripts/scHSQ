\name{getSppLabs}
\alias{getSppLabs}

\title{Queries metadata label for species}

\description{
Queries metadata label for species. This is a column name of pData in the contained cell_data_set objects.
}

\usage{
getSppLabs(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A character corresponding to metadata label for species.}

\author{UnJin Lee}

\examples{\dontrun{
## Query metadata label for species
spp_labs <- getSppLabs(dros_testis)
}}