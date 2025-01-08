\name{getCTLabs}
\alias{getCTLabs}

\title{Queries metadata label for cell type}

\description{
Queries metadata label for cell type. This is a column name of pData in the contained cell_data_set objects.
}

\usage{
getCTLabs(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A character corresponding to metadata label for cell type.}

\author{UnJin Lee}

\examples{\dontrun{
## Query metadata label for cell type
ct_labs <- getCTLabs(dros_testis)
}}