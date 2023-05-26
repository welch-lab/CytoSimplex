<img src="https://github.com/mvfki/scPlotSimplex/raw/main/man/figures/logo.png" width="120">

[![R_CMD_check](https://github.com/mvfki/scPlotSimplex/actions/workflows/R_CMD_check.yml/badge.svg?branch=main)](https://github.com/mvfki/scPlotSimplex/actions/workflows/R_CMD_check.yml)[![codecov](https://codecov.io/gh/mvfki/scPlotSimplex/branch/main/graph/badge.svg?token=AYU2AOE25I)](https://codecov.io/gh/mvfki/scPlotSimplex)[![Seurat](https://img.shields.io/badge/Seurat-4.3.0-green)](https://cran.r-project.org/web/packages/Seurat/index.html)[![sce](https://img.shields.io/badge/SingleCellExperiment-1.22.0-green)](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)

# scPlotSimplex

"scPlotSimplex" is an R package the creates simplex plot showing similarity between single-cells and clusters.

## Installation

For latest developmental version, run the following command in R console:

```R
if (!requireNamespace("devtools", quietly=TRUE)
    install.packages("devtools")
devtools::install_github("mvfki/scPlotSimplex")
```

## Quick Start

In the R package, we provided a small dataset for demonstrating and testing the methods. Users can load them with:

```R
library(scPlotSimplex)
data(rnaRaw)
data(rnaCluster)
```

For high dimensional single-cell transcriptomic data, reducing the dimensionality by selecting top differentially expressed genes for each terminal cluster is recommended. The genes can be selected in the following way:

```R
rnaNorm <- colNormalize(rnaRaw)
gene <- selectTopFeatures(rnaNorm, clusterVar = rnaCluster, vertices = c("OS", "RE", "CH"))
```

Then we can create a demonstrative ternary plot with log-transformed data:

```R
rnaLog <- colNormalize(rnaRaw, scaleFactor = 1e4, log = TRUE)
plotTernary(rnaLog[gene,], clusterVar = rnaCluster, vertices = c("OS", "RE", "CH"))
```
<p align="center">
  <img src="https://github.com/mvfki/scPlotSimplex/raw/main/man/figures/ternary_example.png" alt="Ternary Example"/>
</p>
