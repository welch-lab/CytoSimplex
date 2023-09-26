#' @importFrom rlang .data
#' @import ggplot2
#' @importFrom plot3D lines3D scatter3D text3D getplist arrows3D
#' @importFrom methods show setMethod setOldClass
#' @importFrom Matrix t colSums rowMeans
#' @importFrom Rcpp evalCpp
#' @useDynLib CytoSimplex
NULL

sqrt3 <- sqrt(3)
sqrt6 <- sqrt(6)
triangle <- as.matrix(data.frame(x = c(0, 0.5, 1),
                                 y = c(0, sqrt3/2, 0)))
tetra <- as.matrix(data.frame(x = c(-1/sqrt3, 2*sqrt3/3, -1/sqrt3, 0),
                              y = c(-1, 0, 1, 0),
                              z = c(0, 0, 0, 2*sqrt6/3)))
