# Feature Selection and Dimension Reduction for Single Cell RNA-Seq based on a Multinomial Model

This repository contains supporting code to facilitate reproducible analysis. For details see the biorxiv preprint. If you find bugs please create a github issue. We are working on turning GLM-PCA and its fast approximations into a standalone R package. When it is available we will post a link here.

### Coauthors

Will Townes, Stephanie Hicks, Martin Aryee, and Rafa Irizarry

## Description of Repository Contents

### algs

Implementations of dimension reduction algorithms 
* *existing.R* - wrapper functions for PCA, tSNE, ZINB-WAVE, etc
* *glmpca.R* - implementation of PCA for generalized linear model likelihoods. This method is highlighted in the paper as being suitable for single cell RNA-Seq data.
* *ortho.R* - supporting functions for GLM-PCA that post-process latent factors so that the loadings are orthonormal (just like regular PCA).

### real

Analysis of various real scRNA-Seq datasets. The Rmarkdown files can be used to produce figures in the manuscript

### real_benchmarking

Systematic assessment of clustering performance of a variety of normalization, feature selection, and dimension reduction algorithms using ground-truth datasets.

**[Downloadable table of results from assessments](https://raw.githubusercontent.com/willtownes/scrna2019/master/real_benchmarking/results/cluster_accuracy.txt)**

### util

Utility functions. 

* *clustering.R* - wrappers for seurat clustering, model based clustering, and k-means
* *functions.R* - Poisson and Binomial deviance and residuals functions, a function for loading 10x read counts from molecule information files.
* *functions_genefilter.R* - convenience functions for gene filtering (feature selection) based on highly variable genes, highly expressed genes, and deviance.

