#' Create ternary plots
#' @description
#' Create ternary plots that show similarity between single cells and
#' selected terminals in a triangle 2D space.
#' @param object An object
#' @param ... Arguments passed to other methods.
#' @rdname plotTernary
#' @export plotTernary
#' @details
#' \bold{Argument inheritance} - For matrix/dgCMatrix input object, we first
#' calculate the similarity matrix and obtain a \code{dist.matrix} object. Then
#' we call the \code{plotTernary.dist.matrix} method. For data container objects
#' (e.g. liger), we obtain the correct data matrix first and then call the
#' matrix/dgCMatrix method. Therefore, the arguments inherits as the flow
#' described above.
#'
#' \bold{The calculation of similarity matrix} - TODO
#' @return For distMatrix method, a ggplot object. For other methods, a ggplot
#' object when \code{splitCluster = FALSE}, or a list of ggplot objects when
#' \code{splitCluster = TRUE}.
plotTernary <- function(object, ...) {
    UseMethod('plotTernary', object)
}

#' @rdname plotTernary
#' @param axisPortion Whether to show axis value by proportion, so that the
#' coordinates of each dot sum up to 100. Default \code{TRUE}.
#' @param dotSize,dotColor Dot aesthetics passed to
#' \code{\link[ggplot2]{geom_point}}. Default \code{0.6} and \code{"grey60"}.
#' @param labelColors Colors of the axis lines and vertex labels.
#' Default \code{c("#EE0000FF", "#3B4992FF", "#008B45FF")} (red, blue and green)
#' @param distMethod Method name used for calculating the dist.matrix object.
#' @param title Title text of the plot. Default \code{NULL}.
#' @export
#' @method plotTernary distMatrix
plotTernary.distMatrix <- function(
        object,
        axisPortion = TRUE,
        dotSize = 0.6,
        dotColor = "grey60",
        labelColors = c("#EE0000FF", "#3B4992FF", "#008B45FF"),
        distMethod = attributes(object)$method,
        title = NULL,
        ...
) {
    if (is.null(distMethod) && isFALSE(axisPortion)) {
        stop("Cannot identify distance method for visualizing axis.")
    }
    if (ncol(object) != 4) {
        stop("`distMatrix` object must have four columns for ternary plot, ",
             "where the first three are for vertices and the last for cluster ",
             "assignment.")
    }
    XLAB <- colnames(object)[1]
    YLAB <- colnames(object)[2]
    ZLAB <- colnames(object)[3]
    colnames(object)[1:3] <- c("x", "y", "z")
    x <- y <- z <- NULL
    # NEVER REMOVE "ggplot2::", the imported namespace has problem with
    # rlang ".data" pronoun
    p <- ggtern(object, aes(x, y, z)) +
        geom_point(size = dotSize, stroke = 0.2, color = dotColor) +
        labs(x = XLAB, y = YLAB, z = ZLAB) +
        theme_custom(base_size = 11, base_family = "",
                     col.T = labelColors[2],
                     col.L = labelColors[1],
                     col.R = labelColors[3]) +
        theme(
            tern.panel.background = element_rect(fill = "white"),
            tern.panel.mask.show = FALSE,
            tern.panel.grid.minor = element_line(colour = "grey80"),
            panel.grid.minor = element_line(colour = "grey80")
        )

    # changing the axis legends
    if (isFALSE(axisPortion)) {
        if (distMethod %in% c("pearson", "spearman")) {
            breaks <- seq(0, 2, by = 0.2)
        } else {
            breaks <- seq(0, 1, by = 0.1)
        }
        p <- p +
            scale_L_continuous(breaks = breaks, labels = breaks) +
            scale_T_continuous(breaks = breaks, labels = breaks) +
            scale_R_continuous(breaks = breaks, labels = breaks)
    } else {
        p <- p + custom_percent("Percent")
    }

    return(p)
}

#' @param clusterVar A vector/factor assigning the cluster variable to each
#' column of the matrix object. Or a character scalar selecting a cell metadata
#' variable from container object.
#' @param vertices A vector of THREE values specifying the clusters as the
#' terminals.
#' @param method Distance calculation method. Default \code{"euclidean"}.
#' Choose from \code{"euclidean"}, \code{"cosine"}, \code{"pearson"},
#' \code{"spearman"}.
#' @param normCluster Whether to normalize the distance matrix by clusters.
#' See Details. Default \code{FALSE}.
#' @param scale Whether to min-max scale the distance matrix by clusters.
#' Default \code{TRUE}.
#' @param splitCluster Logical, whether to return a list of plots where each
#' contains cells belonging to only one cluster. Default \code{FALSE}.
#' @param clusterTitle If \code{splitCluster = TRUE}, whether to title each
#' subplot by the name of each cluster. Default \code{TRUE}.
#' @rdname plotTernary
#' @export
#' @method plotTernary default
plotTernary.default <- function(
        object,
        clusterVar,
        vertices,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        normCluster = FALSE,
        scale = TRUE,
        splitCluster = FALSE,
        clusterTitle = TRUE,
        ...
) {
    method <- match.arg(method)
    if (length(unique(vertices)) != 3) {
        stop("Must only specify 3 different vertices")
    }
    distMat <- calcDist(object, clusterVar = clusterVar,
                        vertices = vertices,
                        method = method, normCluster = normCluster,
                        scale = scale)
    if (isFALSE(splitCluster)) plotTernary(object = distMat, ...)
    else {
        if (isTRUE(clusterTitle)) {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotTernary(distMat[distMat$Label == clust,], title = clust,
                           ...)
            })
        } else {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotTernary(distMat[distMat$Label == clust,], ...)
            })
        }
        names(plotList) <- levels(clusterVar)
        return(plotList)
    }
}

#' @rdname plotTernary
#' @param useDatasets For liger method, select datasets where the distance
#' calculation should only be limited within this range. Default \code{NULL}
#' uses all datasets.
#' @param features For container object methods. Valid row subsetting index that
#' selects features. Default \code{NULL}.
#' @export
#' @method plotTernary liger
plotTernary.liger <- function(
        object,
        clusterVar,
        features = NULL,
        useDatasets = NULL,
        ...
) {
    values <- .ligerPrepare(object, clusterVar, features, useDatasets)
    plotTernary(values[[1]], clusterVar = values[[2]], ...)
}

