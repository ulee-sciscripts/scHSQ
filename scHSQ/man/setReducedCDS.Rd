\name{setReducedCDS}
\alias{setReducedCDS}

\title{Set new processed cell_data_set}

\description{
Sets processed cell_data_set.
}

\usage{
setReducedCDS(in_scHSQ, in_CDS)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{in_CDS}{
		a cell_data_set object
		}
}

\value{A scHSQ object with updated cell_data_set.}

\author{UnJin Lee}

\examples{\dontrun{
## Set processed cell_data_set
dros_testis <- setReducedCDS(dros_testis, testis_3sp)
}}