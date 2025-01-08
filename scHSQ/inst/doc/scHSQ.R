### R code from vignette source 'scHSQ.Rnw'

###################################################
### code chunk number 1: scHSQ.Rnw:51-54
###################################################
library(Biobase)
library(monocle3)
library(ggplot2)


###################################################
### code chunk number 2: scHSQ.Rnw:58-60
###################################################
library(scHSQ)
data(testis_3sp)


###################################################
### code chunk number 3: scHSQ.Rnw:67-68
###################################################
dros_testis <- new_scHSQ(testis_3sp, "mel", "spec", "known_type")


###################################################
### code chunk number 4: scHSQ.Rnw:78-80
###################################################
dros_testis <- genMarkGenes(dros_testis)
dros_testis <- applyMarkGenes(dros_testis)


###################################################
### code chunk number 5: scHSQ.Rnw:85-86
###################################################
plotANOVA(dros_testis)


###################################################
### code chunk number 6: scHSQ.Rnw:97-101
###################################################
new_CDS <- getReducedCDS(dros_testis)
gene_list <- getMarkGenes(dros_testis)

ctplots <- plotCTAssign(dros_testis)


###################################################
### code chunk number 7: scHSQ.Rnw:122-131
###################################################
old_clu <- as.character(1:10)
new_clu <- c("Late spermatocytes", "Early spermatocytes", "Early spermatids", 
"Early spermatocytes", "Early spermatids", "Late spermatogonia", "Somatic", 
"Late spermatids", "ananassae spermatocytes", "GSC Early spermatogonia")
cell_types <- data.frame(HSQ_clu=old_clu, HSQ_clu1=new_clu)

dros_testis <- remapCTAssign(dros_testis, cell_types)

ctplots_reassign <- plotCTAssign(dros_testis, "HSQ_clu1")	


###################################################
### code chunk number 8: scHSQ.Rnw:147-149
###################################################
gene_exprs <- plotGeneExprs(dros_testis, "Rbp4", alpha=0.5, 
norm_method="size_only", min_expr=40)


