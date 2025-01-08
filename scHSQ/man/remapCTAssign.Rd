\name{remapCTAssign}
\alias{remapCTAssign}

\title{Remaps cell type assignments in cell type metadata for processed cell_data_set}

\description{
Remaps cell type assignments in cell type metadata for processed cell_data_set. This is useful when transferring known cell type labels onto newly identified clusters across species. ct_map is assumed to have the original cell type labels in the first column, and new labels in the second column. The first colname of ct_map must correlated to a metadata label in pData of the processed cell_data_set in in_scHSQ. The second colname of ct_map will be used to assign a new metadata label in pData of the processed cell_data_set.
}

\usage{
remapCTAssign(in_scHSQ, ct_map)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{ct_map}{
		a data.frame that maps cell type labels in the processed cell_data_set's pData onto new cell type labels.
		}
}

\value{A scHSQ objected with updated cell type labels.}

\author{UnJin Lee}

\examples{\dontrun{
## reassign newly generated clusters onto known cell types
old_clu <- as.character(1:10)
new_clu <- c("Late spermatocytes", "Early spermatocytes", "Early spermatids", "Early spermatocytes", "Early spermatids", "Late spermatogonia", "Somatic", "Late spermatids", "ananassae spermatocytes", "GSC Early spermatogonia")
	
cell_types <- data.frame(HSQ_clu=old_clu, HSQ_clu1=new_clu)
	
dros_testis <- remapCTAssign(dros_testis, cell_types)
}}