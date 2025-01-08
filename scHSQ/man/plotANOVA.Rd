\name{plotANOVA}
\alias{plotANOVA}

\title{Generate plots summarizing ANOVA analysis}

\description{
Generate plots summarizing ANOVA analysis. Plots will show log(F-statistic) values for species-of-origin and cell type assignment on a per gene basis. Plots will also contain lines reflecting HSQ criteria probability thresholds in addition to correlation between log(F-statistic) values. Rasterized plots are recommended to have width=1500, height=1500, res=300, pointsize=5.
}

\usage{
plotANOVA(in_scHSQ)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
}

\value{NULL. Execution of function will output a plot to current R graphics device.}

\author{UnJin Lee}

\examples{\dontrun{
## Generate ANOVA plot

png("scHSQ_ANOVA.png", width=1500, height=1500, res=300, pointsize=5)
plotANOVA(dros_testis)
dev.off()
}}