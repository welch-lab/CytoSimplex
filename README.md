# CytoSimplex <img src="man/figures/logo.png" align="right" width="120" />

[![R_CMD_check](https://github.com/welch-lab/CytoSimplex/actions/workflows/R_CMD_check.yml/badge.svg?branch=main)](https://github.com/welch-lab/CytoSimplex/actions/workflows/R_CMD_check.yml)[![codecov](https://codecov.io/gh/mvfki/CytoSimplex/branch/main/graph/badge.svg?token=AYU2AOE25I)](https://codecov.io/gh/mvfki/CytoSimplex)[![Seurat](https://img.shields.io/badge/Seurat-5.0.1-blue)](https://cran.r-project.org/web/packages/Seurat/index.html)[![sce](https://img.shields.io/badge/SingleCellExperiment-1.22.0-blue)](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)


"CytoSimplex" is an R package that creates simplex plot showing similarity between single-cells and terminals represented by clusters of cells. RNA velocity can be added as another layer of information.

For Python users, we have a Python package [CytoSimplex](https://github.com/welch-lab/pyCytoSimplex) that provides the same functionalities.

## Installation

We **will be** available on cran:

```R
install.packages("CytoSimplex")
```

For latest developmental version, run the following command in R console:

```R
if (!requireNamespace("devtools", quietly = TRUE)
    install.packages("devtools")
devtools::install_github("welch-lab/CytoSimplex")
```

## Quick start

Please checkout our [vignette](https://github.com/welch-lab/CytoSimplex/blob/main/vignettes/CytoSimplex.Rmd) which provides a comprehensive demonstration of the functionalities.

## Citation

If you used CytoSimplex in your work, please cite us with the following work:
