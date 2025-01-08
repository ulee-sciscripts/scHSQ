\name{getReducedCDS}
\alias{getReducedCDS}

\title{Query processed cell_data_set}

\description{
Queries processed cell_data_set after application of HSQ criteria.
}

\usage{
getReducedCDS(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A cell_data_set containing data used for generating cross-species UMAP.}

\author{UnJin Lee}

\examples{\dontrun{
## Query newly processed cell_data_set
new_cds <- getReducedCDS(dros_testis)
}}