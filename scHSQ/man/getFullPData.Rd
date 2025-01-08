\name{getFullPData}
\alias{getFullPData}

\title{Query metadata for original cell_data_set}

\description{
Queries metadata for original cell_data_set saved in its pData slot.
}

\usage{
getFullPData(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A data.frame containing the pData for the original cell_data_set.}

\author{UnJin Lee}

\examples{\dontrun{
## Query metadata for original cell_data_set
dros_pdata <- getFullPData(dros_testis)
}}