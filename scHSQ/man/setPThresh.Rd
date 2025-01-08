\name{setPThresh}
\alias{setPThresh}

\title{Sets probability thresholds for HSQ criteria}

\description{
Sets probability thresholds for HSQ criteria. Useful if ANOVA values are already calculated and thresholds are adjusted.
}

\usage{
setPThresh(in_scHSQ, p_spec=0.5, p_ct=0.5, gen_list=TRUE)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{p_spec}{
		a numeric value between 0 and 1 specifying the threshold probability for species. Defaults to 0.5.
		}
	\item{p_ct}{
		a numeric value between 0 and 1 specifying the threshold probability for cell type. Defaults to 0.5.
		}
	\item{gen_list}{
		a logical for regenerating marker gene list after adjusting thresholds. Defaults to TRUE.
		}
}

\value{A scHSQ object with updated probability thresholds.}

\author{UnJin Lee}

\examples{\dontrun{
## Update probability thresholds to change total number of marker genes
dros_testis <- setPThresh(dros_testis, p_spec=0.75, p_ct=0.75)
}}