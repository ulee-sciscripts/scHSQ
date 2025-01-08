\name{setANOVAStats}
\alias{setANOVAStats}

\title{Sets ANOVA statistics in processed cell_data_set}

\description{
Sets ANOVA statistics in processed cell_data_set slot. Useful for testing.
}

\usage{
setANOVAStats(in_scHSQ, in_anova)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{in_anova}{
		numeric matrix of ANOVA F-statistic values for each gene
		}
}

\value{A scHSQ object with updated ANOVA F-statistic values.}

\author{UnJin Lee}

\examples{\dontrun{
## changes ANOVA F-statistic values
anova_vals <- getANOVAStats(dros_testis)
anova_vals[,1] <- sample(anova_vals[,1], replace=FALSE)
anova_vals[,2] <- sample(anova_vals[,2], replace=FALSE)

dros_testis <- setANOVAStats(dros_testis, anova_vals)
}}