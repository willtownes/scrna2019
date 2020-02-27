# Feature Selection and Dimension Reduction for Single Cell RNA-Seq based on a Multinomial Model

[![DOI](https://zenodo.org/badge/174751869.svg)](https://zenodo.org/badge/latestdoi/174751869)

This repository contains supporting code to facilitate reproducible analysis. For details see the [biorxiv preprint](https://www.biorxiv.org/content/10.1101/574574v1). If you find bugs please create a github issue. 

GLM-PCA (dimension reduction for generalized linear model likelihoods) is now available as a [standalone R package](https://github.com/willtownes/glmpca). This method is highlighted in the paper as being suitable for single cell RNA-Seq data.

The [scry R package](https://github.com/kstreet13/scry) contains functions for feature selection using deviance,
computation of null residuals, and interfaces to apply these methods and GLM-PCA to Bioconductor objects 
such as SingleCellExperiment and SummarizedExperiment.

### Authors

Will Townes, Stephanie Hicks, Martin Aryee, and Rafa Irizarry

## Description of Repository Contents

### algs

Implementations of dimension reduction algorithms 
* *existing.R* - wrapper functions for PCA, tSNE, ZINB-WAVE, etc
* *glmpca.R* - placeholder file that just loads the [glmpca package](https://github.com/willtownes/glmpca).

### real

Analysis of various real scRNA-Seq datasets. The Rmarkdown files can be used to produce figures in the manuscript

### real_benchmarking

Systematic assessment of clustering performance of a variety of normalization, feature selection, and dimension reduction algorithms using ground-truth datasets.

**[Downloadable table of results from assessments](https://raw.githubusercontent.com/willtownes/scrna2019/master/real_benchmarking/results/cluster_accuracy.txt)**

### util

Utility functions. Please consider using the updated versions of these functions via the [scry R package](https://github.com/kstreet13/scry).

* *clustering.R* - wrappers for seurat clustering, model based clustering, and k-means
* *functions.R* - Poisson and Binomial deviance and residuals functions, a function for loading 10x read counts from molecule information files.
* *functions_genefilter.R* - convenience functions for gene filtering (feature selection) based on highly variable genes, highly expressed genes, and deviance.
