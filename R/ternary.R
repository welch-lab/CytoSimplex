#' Create ternary plots
#' @description
#' Create ternary plots that show similarity between single cells and
#' selected three terminals in a ternary baricentric coordinate.
#' @param x Input data. Can be a \code{matrix} or \code{dgCMatrix} object with
#' cells as columns, a \code{Seurat} or \code{SingleCellExperiment} object.
#' "simMat" method takes intermediate values.
#' @param ... Arguments passed to other methods.
#' @rdname plotTernary
#' @export plotTernary
#' @details
#' \bold{Argument inheritance} - For matrix/dgCMatrix ("default" method), we
#' first calculate the similarity matrix and obtain a "simMat" object. Then the
#' "simMat" method is internally called. For data container objects (e.g.
#' Seurat), we obtain the correct data matrix first and then call the "default"
#' method. The arguments inherits as the flow described above.
#'
#' \bold{The calculation of similarity matrix} - The similarity is calculated
#' either by converting a distance metric ("euclidean" or "cosine") with
#' Gaussian kernel, or directly computed with correlation metrics ("pearson" or
#' "spearman"). The centroid of each terminal is obtained first, and the
#' specified metric from each cell to each terminal is calculated. The
#' similarity matrix (n cells by v terminals) is lastly normalized to sum to 1
#' for each cell, so it becomes a baricentric coordinate.
#'
#' \bold{Arrow aesthetics parameters} - The shape of arrows is controlled by 3
#' arguments. Considering an arrow as the combination of a line segment and a
#' triangle, \code{arrowLinewidth} controls the width of the line as well as
#' the edge line of the triangle; \code{arrowAngle} equals to angle of the
#' arrow-tip vertex of the triangle devided by 2 (e.g. the triangle is
#' equilateral when \code{arrowAngle = 20}); \code{arrowLen} controls the
#' absolute length from the arrow-tip vertex to its opposite edge.
#' @return For "simMat" method, a ggplot object. For other methods, a ggplot
#' object when \code{splitCluster = FALSE}, or a list of ggplot objects when
#' \code{splitCluster = TRUE}.
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("OS", "RE", "CH"))
#' plotTernary(rnaRaw, rnaCluster, c("OS", "RE", "CH"), gene)
plotTernary <- function(x, ...) {
    UseMethod('plotTernary', x)
}

#' @param clusterVar A vector/factor assigning the cluster variable to each
#' column of the matrix object. For "Seurat" method, \code{NULL} (default) for
#' \code{Idents(x)}, or a variable name in \code{meta.data} slot. For
#' "SingleCellExperiment" method, \code{NULL} (default) for \code{colLabels(x)},
#' or a variable name in \code{colData} slot.
#' @param vertices Vector of three unique cluster names that will be used for
#' plotting. Or a named list that groups clusters as three terminal vertices.
#' There must not be any overlap between groups.
#' @param features Valid matrix row subsetting index to select features for
#' similarity calculation. Default \code{NULL} uses all available features.
#' @param byCluster Default \code{NULL} to generate one plot with all cells.
#' Set \code{"all"} to split cells in plot by cluster and returns a list of
#' subplots for each cluster as well as the plot including all cells. Otherwise,
#' a vector of cluster names to generate a list of subplots for the specified
#' clusters.
#' @param processed Logical. Whether the input matrix is already processed.
#' \code{TRUE} will bypass internal preprocessing and input matrix will be
#' directly used for similarity calculation. Default \code{FALSE} and raw count
#' input is recommended. If missing in call, using \code{slot = "counts"} in
#' "Seurat" method or using \code{assay.type = "counts"} in
#' "SingleCellExperiment" method will force this argument to be \code{FALSE} and
#' others for \code{TRUE}.
#' @param veloGraph Cell x cell \code{dgCMatrix} object containing velocity
#' information. Shows velocity grid-arrow layer when specified. Default
#' \code{NULL} does not show velocity.
#' @param method Similarity calculation method. Default \code{"euclidean"}.
#' Choose from \code{"euclidean"}, \code{"cosine"}, \code{"pearson"},
#' \code{"spearman"}.
#' @param force Whether to force calculate the similarity when more then 500
#' features are detected, which is generally not recommended. Default
#' \code{FALSE}.
#' @param sigma Gaussian kernel parameter that controls the effect of variance.
#' Only effective when using a distance metric (i.e. \code{method} is
#' \code{"euclidian"} or \code{"cosine"}). Larger values tighten the dot
#' spreading on figure. Default \code{0.08}.
#' @param scale Whether to min-max scale the distance matrix by clusters.
#' Default \code{TRUE}.
#' @param returnData Logical. Whether to return similarity and aggregated
#' velocity data if applicable instead of generating plot. Default \code{FALSE}.
#' @rdname plotTernary
#' @export
#' @method plotTernary default
plotTernary.default <- function(
        x,
        clusterVar,
        vertices,
        features = NULL,
        veloGraph = NULL,
        byCluster = NULL,
        processed = FALSE,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        force = FALSE,
        sigma = 0.08,
        scale = TRUE,
        dotColor = "grey60",
        returnData = FALSE,
        ...
) {
    method <- match.arg(method)
    vcheck <- .checkVertex(x, clusterVar, vertices, n = 3)
    vertClust <- vcheck[[1]]
    vertices <- vcheck[[2]]
    if (length(dotColor) == 1) dotColor <- rep(dotColor, ncol(x))
    if (length(dotColor) != ncol(x)) {
        stop("`dotColor` need to be either 1 scalar or match the number of ",
             "samples in `x`.")
    }
    if (isFALSE(processed) && !is.rawCounts(x)) {
        warning("Input matrix is not raw counts (integers). ",
                "Results may be affected.")
    }
    if (isFALSE(processed)) {
        x <- colNormalize(x, scaleFactor = 1e4, log = TRUE)
    }
    if (!is.null(features)) x <- x[features,]
    simMat <- calcSim(x, clusterVar = vertClust,
                       vertices = vertices, method = method,
                       scale = scale, force = force, sigma = sigma)
    # simMat <- calcDist(x, clusterVar = vertClust,
    #                     vertices = vertices, method = method,
    #                     scale = scale, force = force)

    veloMat <- NULL
    if (!is.null(veloGraph)) {
        if (ncol(veloGraph) != nrow(veloGraph) ||
            !all(rownames(simMat) %in% rownames(veloGraph))) {
            stop("`veloGraph must be of shape N x N and has dimnames covering ",
                 "all cells in `x`.")
        }
        veloGraph <- veloGraph[rownames(simMat), rownames(simMat)]
        veloMat <- aggrVeloGraph(veloGraph, clusterVar = vertClust,
                                 vertices = vertices)
    }

    if (isTRUE(returnData)) return(list(sim = simMat, velo = veloMat))

    if (is.null(byCluster)) {
        return(plotTernary(x = simMat, veloMat = veloMat,
                              dotColor = dotColor, ...))
    } else {
        if (identical(byCluster, "all")) {
            pl <- lapply(levels(clusterVar), function(clust) {
                plotTernary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               title = clust,
                               dotColor = dotColor[clusterVar == clust], ...)
            })
            pl$allCells <- plotTernary(simMat,
                                          veloMat = veloMat,
                                          dotColor = dotColor,
                                          title = "All cells", ...)
            names(pl) <- c(levels(clusterVar), "allCells")
        } else {
            if (any(!byCluster %in% levels(clusterVar))) {
                stop("`byCluster` must be either a vector of cluster name ",
                     "or \"all\".")
            }
            pl <- lapply(byCluster, function(clust) {
                plotTernary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               title = clust,
                               dotColor = dotColor[clusterVar == clust], ...)
            })
            names(pl) <- byCluster
        }
        return(pl)
    }
}

#' @param layer For "Seurat" method, which layer of the assay to be used.
#' Default \code{"counts"}.
#' @param assay For "Seurat" method, the specific assay to get data from.
#' Default \code{NULL} to the default assay.
#' @rdname plotTernary
#' @export
#' @method plotTernary Seurat
#' @examples
#' \donttest{
#' # Seurat example
#' library(Seurat)
#' srt <- CreateSeuratObject(rnaRaw)
#' Idents(srt) <- rnaCluster
#' gene <- selectTopFeatures(srt, vertices = c("OS", "RE", "CH"))
#' plotTernary(srt, features = gene, vertices = c("OS", "RE", "CH"))
#' }
plotTernary.Seurat <- function(
        x,
        layer = "counts",
        assay = NULL,
        clusterVar = NULL,
        processed = FALSE,
        ...
) {
    values <- .getSeuratData(x, layer = layer, assay = assay,
                             clusterVar = clusterVar)
    if (missing(processed)) {
        processed <- layer != "counts"
    }
    plotTernary(values[[1]], clusterVar = values[[2]], processed = processed,
                ...)
}

#' @param assay.type For "SingleCellExperiment" methods. Which assay to use for
#' calculating the similarity. Default \code{"counts"}.
#' @rdname plotTernary
#' @export
#' @method plotTernary SingleCellExperiment
#' @examples
#' \donttest{
#' # SingleCellExperiment example
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#' colLabels(sce) <- rnaCluster
#' gene <- selectTopFeatures(sce, vertices = c("OS", "RE", "CH"))
#' plotTernary(sce, features = gene, vertices = c("OS", "RE", "CH"))
#' }
plotTernary.SingleCellExperiment <- function(
        x,
        assay.type = "counts",
        clusterVar = NULL,
        processed = FALSE,
        ...
) {
    values <- .getSCEData(x, assay.type = assay.type, clusterVar = clusterVar)
    if (missing(processed)) {
        processed <- assay.type != "counts"
    }
    plotTernary(values[[1]], clusterVar = values[[2]], processed = processed,
                ...)
}

# #' @rdname plotTernary
# #' @param useDatasets For liger method, select datasets where the distance
# #' calculation should only be limited within this range. Default \code{NULL}
# #' uses all datasets.
# #' @param features For container object methods. Valid row subsetting index that
# #' selects features. Default \code{NULL}.
# #' @export
# #' @method plotTernary liger
# plotTernary.liger <- function(
#         x,
#         clusterVar,
#         features = NULL,
#         useDatasets = NULL,
#         ...
# ) {
#     values <- .ligerPrepare(x, clusterVar, features, useDatasets)
#     plotTernary(values[[1]], clusterVar = values[[2]], ...)
# }

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
#' @param axisTextShow Logical, whether to show axis text. Default \code{TRUE}.
#' @param axisTextSize,axisTextDrift Similar to the vertex adjustment applied
#' to the text label along the axis breaks. Default size \code{4}, drifted
#' distance \code{0.02}.
#' @param gridLineAlpha Transparency of background grid lines. Default
#' \code{0.6}.
#' @param arrowLinewidth,arrowAngle,arrowLen Arrow aesthetics, see Details.
#' @param titleSize Size of title text. Default \code{14}.
#' @param equilateral Logical, whether to always display the triangle as
#' equilateral. Default \code{TRUE}.
#' @param margin Margin allowed around of the triangle plotting region when
#' \code{equilateral = TRUE}
#' @export
#' @method plotTernary simMat
plotTernary.simMat <- function(
        x,
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
        axisTextShow = TRUE,
        axisTextSize = 4,
        axisTextDrift = 0.02,
        gridLineAlpha = 0.6,
        arrowLinewidth = 0.25,
        arrowAngle = 20,
        arrowLen = 0.2,
        titleSize = 14,
        equilateral = TRUE,
        margin = 0.1,
        ...
) {
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    df <- as.data.frame(as.matrix(x[,1:3]) %*% triangle)
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
        annotate("text", x = -0.01 - vertexLabelDrift, hjust = 0,
                 y = 0.01 - vertexLabelDrift, label = colnames(x)[1],
                 colour = labelColors[1], size = vertexLabelSize,
                 fontface = "bold") +
        annotate("text", x = 0.5, y = 3^0.5/2 + vertexLabelDrift, hjust = 0.5,
                 label = colnames(x)[2], colour = labelColors[2],
                 size = vertexLabelSize, fontface = "bold") +
        annotate("text", x = 1.01 + vertexLabelDrift, hjust = 1,
                 y = 0.01 - vertexLabelDrift, label = colnames(x)[3],
                 colour = labelColors[3], size = vertexLabelSize,
                 fontface = "bold") +
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

            # Adding grid line
            annotate("segment", x = triangle$x[2] + vecTLBreak[1]*i,
                     xend = triangle$x[3] - vecLRBreak[1]*i,
                     y = triangle$y[2] + vecTLBreak[2]*i,
                     yend = triangle$y[3] - vecLRBreak[2]*i,
                     linetype = 4, color = labelColors[1],
                     alpha = gridLineAlpha) +
            annotate("segment", x = triangle$x[3] + vecRTBreak[1]*i,
                     xend = triangle$x[1] - vecTLBreak[1]*i,
                     y = triangle$y[3] + vecRTBreak[2]*i,
                     yend = triangle$y[1] - vecTLBreak[2]*i,
                     linetype = 4, color = labelColors[2],
                     alpha = gridLineAlpha) +
            annotate("segment", x = triangle$x[1] + vecLRBreak[1]*i,
                     xend = triangle$x[2] - vecRTBreak[1]*i,
                     y = triangle$y[1] + vecLRBreak[2]*i,
                     yend = triangle$y[2] - vecRTBreak[2]*i,
                     linetype = 4, color = labelColors[3],
                     alpha = gridLineAlpha)
        if (axisTextShow) {
            ternCoord <- ternCoord +
                annotate("text", label = perc,
                         x = triangle$x[2] + vecTLBreak[1]*i - axisTextDrift,
                         y = triangle$y[2] + vecTLBreak[2]*i + axisTextDrift,
                         color = labelColors[1], size = axisTextSize) +

                annotate("text", label = perc,
                         x = triangle$x[3] + vecRTBreak[1]*i + axisTextDrift,
                         y = triangle$y[3] + vecRTBreak[2]*i + axisTextDrift,
                         color = labelColors[2], size = axisTextSize) +

                annotate("text", label = perc,
                         x = triangle$x[1] + vecLRBreak[1]*i,
                         y = triangle$y[1] + vecLRBreak[2]*i - axisTextDrift,
                         color = labelColors[3], size = axisTextSize)
        }
    }
    p <- ternCoord + geom_point(color = dotColor, size = dotSize)

    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = x, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            subcoords <- arrowCoords[[i]]
            p <- p +
                annotate("segment", x = subcoords[,1], y = subcoords[,2],
                         xend = subcoords[,3], yend = subcoords[,4],
                         color = labelColors[i], linewidth = arrowLinewidth,
                         arrow = arrow(angle = arrowAngle,
                                       length = unit(arrowLen, "cm"),
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
