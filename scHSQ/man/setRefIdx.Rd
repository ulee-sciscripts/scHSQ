\name{setRefIdx}
\alias{setRefIdx}

\title{Sets indices for cells in the reference species}

\description{
Sets indices for cells in the reference species. Useful for excluding particular cells or cell type classifications for applying HSQ criteria.
}

\usage{
setRefIdx(in_scHSQ, idx)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{idx}{
		numeric vector of indices
		}
}

\value{A scHSQ object with updated reference species indices.}

\author{UnJin Lee}

\examples{\dontrun{
## Set indices for reference species' cells
new_idx <- 1:2500
dros_testis <- setRefIdx(dros_testis, new_idx)
}}