# scHSQ 

scHSQ is an R package for applying the HSQ criteria on single cell RNA-Seq data [Lee, et al. 2025](https://www.biorxiv.org/content/10.1101/2024.04.29.591771v2.full). This methodology is useful for transferring known cell type assignment across the same tissue that has undergone a large degree of cell type-specific evolutionary divergence. In these cases, the effect of transcriptome divergence is so strong that it masks cell type effects on global gene expression. 

To tackle this issue, genes passing the HSQ criteria have high cell type-specific effects on expression with low species-specfic effects. This is calculated using ANOVA F-statistics. 

Using this package, one can identify genes passing the HSQ criteria, project cells onto a single cross-species manifold (UMAP), perform clustering, transfer cell type labels from a reference species onto other species, and visualize gene expression on the cross-species UMAP.

Instructions for running the software package on a sample data set is found in the [vignette](https://github.com/ulee-sciscripts/scHSQ/blob/main/scHSQ/inst/doc/scHSQ.pdf). You may also find the [manual](https://github.com/ulee-sciscripts/scHSQ/blob/main/scHSQ/vignettes/scHSQ-manual.pdf) to be helpful.

## Installation

Note: R, Bioconductor, monocle3, and ggplot2 must be installed before using scHSQ.

scHSQ may be installed either through source:

```bash
wget https://github.com/ulee-sciscripts/scHSQ/blob/main/scHSQ/scHSQ_1.0.0.tar.gz
tar -xzvf scHSQ_1.0.0.tar.gz
R CMD INSTALL scHSQ
```

or through remotes' install_github function:
```bash
Rscript -e 'remotes::install_github("ulee-sciscripts/scHSQ/scHSQ")'
```

or 

```R
remotes::install_github("ulee-sciscripts/scHSQ/scHSQ")
```

# Docker image

A dockerfile is provided in this package for installation. A docker image containing R and all dependencies is hosted on [dockerhub](https://hub.docker.com/r/avianalter/schsq) as well.

```bash
docker pull avianalter/schsq:1.0.0
```


Please cite:
U. Lee, C. Li, C.B. Langer, N. Svetec, & L. Zhao, Comparative single-cell analysis of transcriptional bursting reveals the role of genome organization in de novo transcript origination, Proc. Natl. Acad. Sci. U.S.A. 122 (18) e2425618122, https://doi.org/10.1073/pnas.2425618122 (2025).
