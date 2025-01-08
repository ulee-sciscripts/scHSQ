\name{setCellType}
\alias{setCellType}

\title{Sets new metadata label for updated cell type information}

\description{
Sets new metadata label for updated cell type information. Updated information is a new metadata label in contained processed cell_data_set object using the ct_field_new label. Names of ct_labels should correspond to UMIs.
}

\usage{
setCellType(in_scHSQ, ct_labels, ct_field_new)
}


\arguments{
	\item{in_scHSQ}{
		a scHSQ object
		}
	\item{ct_labels}{
		a character vector of updated cell type labels
		}
	\item{ct_field_new}{
		a character representing the metadata field for updated cell type labels
		}
}

\value{A scHSQ object with updated cell type labels in the processed cell_data_set.}

\author{UnJin Lee}

\examples{\dontrun{
## Set new cell type information
new_ct_labs <- rep(as.character(1:5000), 3)
names(new_ct_labs) <- rownames(getFullPData(dros_testis))

dros_testis <- setCellType(dros_testis, new_ct_labs, "new_labs")
}}