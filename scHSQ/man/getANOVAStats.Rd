\name{getANOVAStats}
\alias{getANOVAStats}

\title{Queries ANOVA statistics in processed cell_data_set}

\description{
Queries ANOVA statistics in processed cell_data_set slot. Useful for testing.
}

\usage{
getANOVAStats(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{A matrix with ANOVA F-statistic values for each gene.}

\author{UnJin Lee}

\examples{\dontrun{
## changes ANOVA F-statistic values
anova_vals <- getANOVAStats(dros_testis)
anova_vals[,1] <- sample(anova_vals[,1], replace=FALSE)
anova_vals[,2] <- sample(anova_vals[,2], replace=FALSE)

dros_testis <- setANOVAStats(dros_testis, anova_vals)
}}