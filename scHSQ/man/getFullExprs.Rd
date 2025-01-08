\name{getFullExprs}
\alias{getFullExprs}

\title{Query read counts for original cell_data_set}

\description{
Queries read counts for original cell_data_set saved in its exprs slot.
}

\usage{
getFullExprs(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A sparse matrix containing the read counts for the original cell_data_set.}

\author{UnJin Lee}

\examples{\dontrun{
## Query read counts for original cell_data_set
dros_exprs <- getFullExprs(dros_testis)
}}