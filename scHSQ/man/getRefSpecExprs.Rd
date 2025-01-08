\name{getRefSpecExprs}
\alias{getRefSpecExprs}

\title{Queries read counts from cells in the reference species}

\description{
Queries read counts from cells in the reference species.
}

\usage{
getRefSpecExprs(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A sparse matrix with read counts from reference species cells.}

\author{UnJin Lee}

\examples{\dontrun{
## Query read counts from reference species cells
mel_exprs <- getRefSpecExprs(dros_testis)
}}