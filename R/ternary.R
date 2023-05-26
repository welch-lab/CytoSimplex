#' Create ternary plots
#' @description
#' Create ternary plots that show similarity between single cells and
#' selected three terminals in a triangle 2D space.
#' @param object An object
#' @param ... Arguments passed to other methods.
#' @rdname plotTernary
#' @export plotTernary
#' @details
#' \bold{Argument inheritance} - For matrix/dgCMatrix input object ("default"
#' method), we first calculate the similarity matrix and obtain a "simMat"
#' object. Then we call the "simMat" method. For data container objects
#' (e.g. Seurat), we obtain the correct data matrix first and then call the
#' "default" method. Therefore, the arguments inherits as the flow
#' described above.
#'
#' \bold{The calculation of similarity matrix} - TODO
#' @return For "simMat" method, a ggplot object. For other methods, a ggplot
#' object when \code{splitCluster = FALSE}, or a list of ggplot objects when
#' \code{splitCluster = TRUE}.
#' @examples
#' rnaNorm <- colNormalize(rnaRaw)
#' gene <- selectTopFeatures(rnaNorm, rnaCluster, c("OS", "RE", "CH"))
#' rnaLog <- colNormalize(rnaRaw, 1e4, TRUE)
#' plotTernary(rnaLog[gene, ], rnaCluster, c("OS", "RE", "CH"))
plotTernary <- function(object, ...) {
    UseMethod('plotTernary', object)
}

#' @rdname plotTernary
#' @param title Title text of the plot. Default \code{NULL}.
#' @param veloMat Aggregated velocity matrix. Output of \code{aggrVeloGraph}.
#' @param nGrid Number of grids along the bottom side of the equilateral
#' triangle. Default \code{10}.
#' @param radius Arrow length of unit velocity. Lower this when arrows point
#' outside of the coordinate. Default \code{0.1}.
#' @param dotSize,dotColor Dot aesthetics passed to
#' \code{\link[ggplot2]{geom_point}}. Default \code{0.6} and \code{"grey60"}.
#' @param labelColors Colors of the axis lines and vertex labels.
#' Default \code{c("#3B4992FF", "#EE0000FF", "#008B45FF")} (blue, red and green)
#' @param vertexLabelSize,vertexLabelDrift Adjustment on the three vertex text
#' labels. Drift means the distance that the labels should be moved against the
#' center of the plot. Default size \code{5}, drifted distance \code{0.03}.
#' @param axisBreak Number of breaks to be labeled along axis. Default
#' \code{5}.
#' @param axisTextSize,axisTextDrift Similar to the vertex adjustment applied
#' to the text label along the axis breaks. Default size \code{4}, drifted
#' distance \code{0.02}.
#' @param gridLineAlpha Transparency of background grid lines. Default
#' \code{0.6}.
#' @param arrowLinewidth Arrow aesthetics. Default \code{0.25}.
#' @param titleSize Size of title text. Default \code{14}.
#' @param equilateral Logical, whether to always display the triangle as
#' equilateral. Default \code{TRUE}.
#' @param margin Margin allowed around of the triangle plotting region when
#' \code{equilateral = TRUE}
#' @export
#' @method plotTernary simMat
plotTernary.simMat <- function(
        object,
        title = NULL,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.1,
        dotSize = 0.6,
        dotColor = "grey60",
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF"),
        vertexLabelSize = 5,
        vertexLabelDrift = 0.03,
        axisBreak = 5,
        axisTextSize = 4,
        axisTextDrift = 0.02,
        gridLineAlpha = 0.6,
        arrowLinewidth = 0.25,
        titleSize = 14,
        equilateral = TRUE,
        margin = 0.1,
        ...
) {
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    df <- as.data.frame(as.matrix(object[,1:3]) %*% triangle)
    triangle <- as.data.frame(triangle)
    ternCoord <- ggplot(df, aes(x = .data$x, y = .data$y)) +
        annotate("segment", x = triangle$x[1], xend = triangle$x[2],
                          y = triangle$y[1], yend = triangle$y[2],
                          colour = labelColors[1]) +
        annotate("segment", x = triangle$x[2], xend = triangle$x[3],
                          y = triangle$y[2], yend = triangle$y[3],
                          colour = labelColors[2]) +
        annotate("segment", x = triangle$x[3], xend = triangle$x[1],
                          y = triangle$y[3], yend = triangle$y[1],
                          colour = labelColors[3]) +
        annotate("text", x = -0.01 - vertexLabelDrift,
                 y = 0.01 - vertexLabelDrift, label = colnames(object)[1],
                 colour = labelColors[1], size = vertexLabelSize) +
        annotate("text", x = 0.5, y = 3^0.5/2 + vertexLabelDrift,
                 label = colnames(object)[2], colour = labelColors[2],
                 size = vertexLabelSize) +
        annotate("text", x = 1.01 + vertexLabelDrift,
                 y = 0.01 - vertexLabelDrift, label = colnames(object)[3],
                 colour = labelColors[3], size = vertexLabelSize) +
        labs(title = title) +
        theme_void() +
        theme(plot.title = element_text(face = "bold", size = titleSize,
                                        hjust = 0.5))

    vecTLBreak <- c(triangle$x[1] - triangle$x[2],
                    triangle$y[1] - triangle$y[2]) / axisBreak
    vecRTBreak <- c(triangle$x[2] - triangle$x[3],
                    triangle$y[2] - triangle$y[3]) / axisBreak
    vecLRBreak <- c(triangle$x[3] - triangle$x[1],
                    triangle$y[3] - triangle$y[1]) / axisBreak
    for (i in seq(axisBreak - 1)) {
        perc <- round(100 / axisBreak * i, digits = 1)
        ternCoord <- ternCoord +
            # Adding axis break number label
            annotate("text", label = perc,
                     x = triangle$x[2] + vecTLBreak[1]*i - axisTextDrift,
                     y = triangle$y[2] + vecTLBreak[2]*i + axisTextDrift,
                     color = labelColors[1], size = axisTextSize) +
            # Adding grid line
            annotate("segment", x = triangle$x[2] + vecTLBreak[1]*i,
                     xend = triangle$x[3] - vecLRBreak[1]*i,
                     y = triangle$y[2] + vecTLBreak[2]*i,
                     yend = triangle$y[3] - vecLRBreak[2]*i,
                     linetype = 4, color = labelColors[1],
                     alpha = gridLineAlpha) +

            annotate("text", label = perc,
                     x = triangle$x[3] + vecRTBreak[1]*i + axisTextDrift,
                     y = triangle$y[3] + vecRTBreak[2]*i + axisTextDrift,
                     color = labelColors[2], size = axisTextSize) +
            annotate("segment", x = triangle$x[3] + vecRTBreak[1]*i,
                     xend = triangle$x[1] - vecTLBreak[1]*i,
                     y = triangle$y[3] + vecRTBreak[2]*i,
                     yend = triangle$y[1] - vecTLBreak[2]*i,
                     linetype = 4, color = labelColors[2],
                     alpha = gridLineAlpha) +

            annotate("text", label = perc,
                     x = triangle$x[1] + vecLRBreak[1]*i,
                     y = triangle$y[1] + vecLRBreak[2]*i - axisTextDrift,
                     color = labelColors[3], size = axisTextSize) +
            annotate("segment", x = triangle$x[1] + vecLRBreak[1]*i,
                     xend = triangle$x[2] - vecRTBreak[1]*i,
                     y = triangle$y[1] + vecLRBreak[2]*i,
                     yend = triangle$y[2] - vecRTBreak[2]*i,
                     linetype = 4, color = labelColors[3],
                     alpha = gridLineAlpha)

    }
    p <- ternCoord + geom_point(color = dotColor, size = dotSize)

    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = object, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            subcoords <- arrowCoords[[i]]
            p <- p +
                annotate("segment", x = subcoords[,1], y = subcoords[,2],
                         xend = subcoords[,3], yend = subcoords[,4],
                         color = labelColors[i],
                         arrow = arrow(angle = 20, length = unit(.1, "cm"),
                                       type = "closed"))
        }
    }
    if (isTRUE(equilateral)) {
        suppressMessages({
            p <- p +
                coord_fixed(xlim = c(0 - margin, 1 + margin),
                            ylim = c(0 - margin, 3^0.5/2 + margin))
        })
    }
    return(p)
}

#' @param clusterVar A vector/factor assigning the cluster variable to each
#' column of the matrix object. For "Seurat" method, leave \code{NULL} for using
#' default "Idents", or can also be a variable in \code{meta.data} slot. For
#' "SingleCellExperiment" method, leave \code{NULL} for using "colLabels", or
#' can be a variable in \code{colData} slot.
#' @param vertices Vector of three unique cluster names that will be used for
#' plotting. Or a named list that groups clusters as three terminal vertices.
#' There must not be any overlap between groups.
#' @param veloGraph Cell x cell dgCMatrix object containing velocity
#' information. Shows velocity grid-arrow layer when specified. Default
#' \code{NULL} does not show velocity.
#' @param method Distance calculation method. Default \code{"euclidean"}.
#' Choose from \code{"euclidean"}, \code{"cosine"}, \code{"pearson"},
#' \code{"spearman"}.
#' @param force Whether to force calculate the distance when more then 500
#' features are detected, which is generally not recommended. Default
#' \code{FALSE}.
#' @param sigma Gaussian kernel parameter that controls the effect of variance.
#' Only effective when using a distance metric (i.e. \code{method} is
#' \code{"euclidian"} or \code{"cosine"}). Larger value tighten the dot
#' spreading on figure. Default \code{0.08}.
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
        veloGraph = NULL,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        force = FALSE,
        #distKernel = c("gaussian", "log"),
        sigma = 0.08,
        scale = TRUE,
        splitCluster = FALSE,
        clusterTitle = TRUE,
        dotColor = "grey60",
        ...
) {
    method <- match.arg(method)
    vcheck <- .checkVertex(object, clusterVar, vertices, n = 3)
    vertClust <- vcheck[[1]]
    vertices <- vcheck[[2]]
    if (length(dotColor) == 1) dotColor <- rep(dotColor, ncol(object))
    if (length(dotColor) != ncol(object)) {
        stop("`dotColor` need to be either 1 scalar or match the number of ",
             "samples in `object`.")
    }
    simMat <- calcSim(object, clusterVar = vertClust,
                       vertices = vertices, method = method,
                       scale = scale, force = force, sigma = sigma)
    # simMat <- calcDist(object, clusterVar = vertClust,
    #                     vertices = vertices, method = method,
    #                     scale = scale, force = force)

    veloMat <- NULL
    if (!is.null(veloGraph)) {
        if (ncol(veloGraph) != nrow(veloGraph) ||
            !all(rownames(simMat) %in% rownames(veloGraph))) {
            stop("`veloGraph must be of shape N x N and has dimnames covering ",
                 "all cells in `object`.")
        }
        veloGraph <- veloGraph[rownames(simMat), rownames(simMat)]
        veloMat <- aggrVeloGraph(veloGraph, clusterVar = vertClust,
                                 vertices = vertices)
    }

    if (isFALSE(splitCluster)) plotTernary(object = simMat,
                                           veloMat = veloMat,
                                           dotColor = dotColor,
                                           ...)
    else {
        if (isTRUE(clusterTitle)) {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotTernary(simMat[clusterVar == clust,],
                            veloMat = veloMat[clusterVar == clust,],
                            dotColor = dotColor[clusterVar == clust],
                            title = clust, ...)
            })
            plotList$allCells <- plotTernary(simMat,
                                             veloMat = veloMat,
                                             dotColor = dotColor,
                                             title = "All cells", ...)
        } else {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotTernary(simMat[clusterVar == clust,],
                            veloMat = veloMat[clusterVar == clust,],
                            dotColor = dotColor[clusterVar == clust],
                            ...)
            })
            plotList$allCells <- plotTernary(simMat,
                                             veloMat = veloMat,
                                             dotColor = dotColor, ...)
        }
        names(plotList) <- c(levels(clusterVar), "allCells")
        return(plotList)
    }
}

#' @param features For "Seurat" method. Valid row subsetting index that
#' selects features. Default \code{NULL} uses all available features.
#' @param slot For "Seurat" method, choose from \code{"data"},
#' \code{"scale.data"} or \code{"counts"}. Default \code{"data"}.
#' @param assay For "Seurat" method, the specific assay to get data from.
#' Default \code{NULL} to the default assay.
#' @rdname plotTernary
#' @export
#' @method plotTernary Seurat
#' @examples
#'
#' # Seurat example
#' if (FALSE) {
#'     srt <- CreateSeuratObject(rnaRaw)
#'     Idents(srt) <- rnaCluster
#'     srt <- colNormalize(srt)
#'     gene <- selectTopFeatures(srt, vertices = c("OS", "RE", "CH"))
#'     srt <- colNormalize(srt, scaleFactor = 1e4, log = TRUE)
#'     plotTernary(srt, features = geneSel, vertices = c("OS", "RE", "CH"))
#' }
plotTernary.Seurat <- function(
        object,
        features = NULL,
        slot = "data",
        assay = NULL,
        clusterVar = NULL,
        ...
) {
    values <- .getSeuratData(object, features = features,
                             slot = slot, assay = assay,
                             clusterVar = clusterVar)
    plotTernary(values[[1]], clusterVar = values[[2]], ...)
}

#' @param subset.row For "SingleCellExperiment" methods. Valid row subsetting
#' index that selects features. Default \code{NULL} uses all available features.
#' @param assay.type For "SingleCellExperiment" methods. Which assay to use for
#' calculating the similarity. Default \code{"logcounts"}.
#' @rdname plotTernary
#' @export
#' @method plotTernary SingleCellExperiment
#' @examples
#'
#' # SingleCellExperiment example
#' if (FALSE) {
#'     library(SingleCellExperiment)
#'     sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#'     colLabels(sce) <- rnaCluster
#'     sce <- colNormalize(sce)
#'     gene <- selectTopFeatures(sce, vertices = c("OS", "RE", "CH"))
#'     sce <- colNormalize(sce, scaleFactor = 1e4, log = TRUE)
#'     plotTernary(sce, subset.row = gene, vertices = c("OS", "RE", "CH"))
#' }
plotTernary.SingleCellExperiment <- function(
        object,
        assay.type = "logcounts",
        subset.row = NULL,
        clusterVar = NULL,
        ...
) {
    values <- .getSCEData(object, subset.row = subset.row,
                          assay.type = assay.type,
                          clusterVar = clusterVar)
    plotTernary(values[[1]], clusterVar = values[[2]], ...)
}

#' #' @rdname plotTernary
#' #' @param useDatasets For liger method, select datasets where the distance
#' #' calculation should only be limited within this range. Default \code{NULL}
#' #' uses all datasets.
#' #' @param features For container object methods. Valid row subsetting index that
#' #' selects features. Default \code{NULL}.
#' #' @export
#' #' @method plotTernary liger
#' plotTernary.liger <- function(
#'         object,
#'         clusterVar,
#'         features = NULL,
#'         useDatasets = NULL,
#'         ...
#' ) {
#'     values <- .ligerPrepare(object, clusterVar, features, useDatasets)
#'     plotTernary(values[[1]], clusterVar = values[[2]], ...)
#' }
