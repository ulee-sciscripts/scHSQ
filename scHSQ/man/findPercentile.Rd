\name{findPercentile}
\alias{findPercentile}

\title{Identification of empirical cutoff for probability thresholds}

\description{
Identification of empirical cutoff for probability thresholds. Performed through simple sorting. Equivalent to the inverse of the empirical cumulative distribution function.
}

\usage{
findPercentile(dat, p)
}


\arguments{
	\item{dat}{
		numeric vector of empirical data
		}
	\item{p}{
		a numeric value between 0 and 1 specifying the threshold probability.
		}
}

\value{A numeric value representing the empirical distribution correlating to the threshold probability.}

\author{UnJin Lee}

\examples{\dontrun{
## find empirical value corresponding to probability threshold
dist <- rnorm(1000, 0, 1)
median_thresh <- findPercentile(dist, 0.5)
}}