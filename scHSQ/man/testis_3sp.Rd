\name{testis_3sp}
\docType{data}
\alias{testis_3sp}
\title{Breast Cancer 443 Data Set}

\description{
testis_3sp is an monocle3 cell_data_set object that contains scRNA-Seq counts for testis tissue obtained from D. melanogaster, D. yakuba, and D. ananassae. Only genes with one-to-one and reciprocal best orthology were retained. Data was merged using monocle3's combine_cds function. Minimal processing occurred, as monocle3's preprocess_cds function was applied. Additionally, while monocle3's reduce_dimension and cluster_cells functions were applied to this data set, there is no impact on downstream analysis performed in scHSQ.
}

\usage{data(testis_3sp)}
\format{cell_data_set with expression data and phenotypic meta-data}

