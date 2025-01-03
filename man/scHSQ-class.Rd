\name{scHSQ-class}
\Rdversion{1.1}
\docType{class}
\alias{scHSQ-class}
\alias{new_scHSQ,cell_data_set,character,character,character-method}
\alias{genMarkGenes,scHSQ-method}
\alias{getMarkGenes,scHSQ-method}
\alias{getFullCDS,scHSQ-method}
\alias{getReducedCDS,scHSQ-method}
\alias{setReducedCDS,scHSQ,cell_data_set-method}
\alias{getFullPData,scHSQ-method}
\alias{getFullExprs,scHSQ-method}
\alias{getRefIdx,scHSQ-method}
\alias{setRefIdx,scHSQ,numeric-method}
\alias{getRefSpecExprs,scHSQ-method}
\alias{getSppLabs,scHSQ-method}
\alias{getCTLabs,scHSQ-method}
\alias{setCellType,scHSQ,character,character-method}
\alias{getRefCellType,scHSQ-method}
\alias{getSpp,scHSQ-method}
\alias{setPThresh,scHSQ-method}
\alias{findPercentile,numeric,numeric-method}
\alias{setANOVAStats,scHSQ,matrix-method}
\alias{getANOVAStats,scHSQ-method}
\alias{calcANOVA,scHSQ-method}
\alias{remapCTAssign,scHSQ,data.frame-method}
\alias{calcANOVACor,scHSQ-method}
\alias{plotANOVA,scHSQ-method}
\alias{applyMarkGenes,scHSQ-method}
\alias{reclusterCDS,scHSQ-method}
\alias{plotCTAssign,scHSQ-method}
\alias{plotGeneExprs,scHSQ,character-method}


\title{Class \code{"scHSQ"}}
\description{
The scHSQ object contains the necessary elements defining the data set on which the HSQ criteria will be applied.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new_scHSQ("cell_data_set", "character", "character", "character", ...)}.
Objects all contain 14 slots - fullCDS, n_species, n_celltype, ref_spp, spp_labs, idx_ref, celltype_labs, col_scheme, p_species, p_celltype, anova_stats, marker_list, reducedCDS, hsq_status.
}
\section{Slots}{
  \describe{
    \item{\code{fullCDS}:}{Object of class \code{"cell_data_set"} ~~ }
    \item{\code{n_species}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n_celltype}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ref_spp}:}{Object of class \code{"character"} ~~ }
	\item{\code{spp_labs}:}{Object of class \code{"character"} ~~ }
	\item{\code{idx_ref}:}{Object of class \code{"numeric"} ~~ }
	\item{\code{celltype_labs}:}{Object of class \code{"character"} ~~ }
	\item{\code{col_scheme}:}{Object of class \code{"character"} ~~ }
	\item{\code{p_species}:}{Object of class \code{"numeric"} ~~ }
	\item{\code{p_celltype}:}{Object of class \code{"numeric"} ~~ }
	\item{\code{anova_stats}:}{Object of class \code{"matrix"} ~~ }
	\item{\code{marker_list}:}{Object of class \code{"character"} ~~ }
	\item{\code{reducedCDS}:}{Object of class \code{"cell_data_set"} ~~ }
	\item{\code{hsq_status}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{new_scHSQ}{\code{signature(in_CDS = "cell_data_set", ref_species = "character", mdat_lab_spec = "character", mdat_lab_ct = "character")}: ... }
    \item{genMarkGenes}{\code{signature(in_scHSQ = "scHSQ")}: ... }
    \item{getMarkGenes}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getFullCDS}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getReducedCDS}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{setReducedCDS}{\code{signature(in_scHSQ = "scHSQ", in_CDS = "cell_data_set")}: ... }
	\item{getFullPData}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getFullExprs}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getRefIdx}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{setRefIdx}{\code{signature(in_scHSQ = "scHSQ", idx = "numeric")}: ... }
	\item{getRefSpecExprs}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getSppLabs}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getCTLabs}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{setCellType}{\code{signature(in_scHSQ = "scHSQ", ct_labels = "character", ct_field_new = "character")}: ... }
	\item{getRefCellType}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{getSpp}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{setPThresh}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{findPercentile}{\code{signature(dat = "numeric", p = "numeric")}: ... }
	\item{setANOVAStats}{\code{signature(in_scHSQ = "scHSQ", in_anova="matrix")}: ... }
	\item{getANOVAStats}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{calcANOVA}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{remapCTAssign}{\code{signature(in_scHSQ = "scHSQ", ct_map = "data.frame")}: ... }
	\item{calcANOVACor}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{plotANOVA}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{applyMarkGenes}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{reclusterCDS}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{plotCTAssign}{\code{signature(in_scHSQ = "scHSQ")}: ... }
	\item{plotGeneExprs}{\code{signature(in_scHSQ = "scHSQ", gene = "character")}: ... }
	 }
}
\references{
Lee U, Li C, Langer CB, Svetec N, Zhao L (2025) Comparative Single Cell Analysis of Transcriptional Bursting Reveals the Role of Genome Organization on de novo Transcript Origination. bioRxiv doi:10.1101/2024.04.29.591771
}
\author{
UnJin Lee
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("scHSQ")
}
\keyword{classes}
