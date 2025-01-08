\name{new_scHSQ}
\alias{new_scHSQ}

\title{Initialize new scHSQ object}

\description{
Allows for the generation of a new scHSQ object containing relevant information and slots for applying the HSQ criteria on monocle3 cell_data_set objects.
}

\usage{
new_scHSQ(in_CDS, ref_species, mdat_lab_spec, mdat_lab_ct, color_scheme="default", p_spec=0.5, p_ct=0.5)
}


\arguments{
	\item{in_CDS}{
		monocle3 cell_data_set object with pData labels for species and cell type
		}
	\item{ref_species}{
		a character reflecting the reference species. ref_species should be a level of the pData label for species-of-origin
		}
	\item{mdat_lab_spec}{
		a character reflecting the metadata label for species information. mdat_lab_spec should be a column name of the pData of in_CDS
		}
	\item{mdat_lab_ct}{
		a character reflecting the metadata label for cell type information. mdat_lab_ct should be a column name of the pData of in_CDS
		}
	\item{color_scheme}{
		a character vector specifying a valid ggplot2 color scheme for visualization of cell type classifications. Leave blank for default color schemes.
		}
	\item{p_spec}{
		a numeric value between 0 and 1 specifying the threshold probability for species. Defaults to 0.5.
		}
	\item{p_ct}{
		a numeric value between 0 and 1 specifying the threshold probability for cell type. Defaults to 0.5.
		}
}

\value{A new scHSQ object that may be used for further processing.}

\author{UnJin Lee}

\examples{\dontrun{
require(Biobase, monocle3)
## Generate new scHSQ object
data(testis_3sp)
dros_testis <- new_scHSQ(testis_3sp, "mel", "spec", "known_type")
}}