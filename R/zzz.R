#' @importFrom rlang .data
#' @import ggplot2
#' @importFrom plot3D lines3D scatter3D text3D getplist
#' @importFrom methods show setMethod setOldClass
#' @useDynLib scPlotSimplex
NULL

triangle <- data.frame(x = c(0, 0.5, 1),
                       y = c(0, 3^0.5/2, 0))
