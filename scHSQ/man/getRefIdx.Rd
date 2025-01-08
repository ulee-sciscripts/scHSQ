\name{getRefIdx}
\alias{getRefIdx}

\title{Query indices for cells in the reference species}

\description{
Queries indices for cells in the reference species.
}

\usage{
getRefIdx(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A numeric vector corresponding to the indices for cells derived from the original reference species.}

\author{UnJin Lee}

\examples{\dontrun{
## Query indices for reference species' cells
mel_idx <- getRefIdx(dros_testis)
}}