#' Create ternary simplex plots
#' @description
#' Create ternary plots that show similarity between single cells and
#' selected three terminals in a ternary baricentric coordinate.
#' @param x Input data. Can be a \code{matrix} or \code{dgCMatrix} object with
#' cells as columns, a \code{Seurat} or \code{SingleCellExperiment} object.
#' "simMat" method takes intermediate values.
#' @param ... Arguments passed to other S3 methods.
#' @rdname plotTernary
#' @export plotTernary
#' @family plotTernary
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
#' @return By default, a "ggplot" object when \code{byCluster} is not specified,
#' a list of "ggplot" object when \code{byCluster} is specified. When
#' \code{interactive = TRUE}, a "plotly" object is returned. When
#' \code{returnData = TRUE}, a list of similarity matrix and aggregated velocity
#' matrix is returned.
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("OS", "RE", "CH"))
#' plotTernary(rnaRaw, rnaCluster, c("OS", "RE", "CH"), gene)
plotTernary <- function(x, ...) {
    UseMethod('plotTernary', x)
}

#' @rdname plotTernary
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
#' @param dotColorBy A vector/factor for coloring dots, can be either categorical
#' (must be character or factor) or continuous. Default \code{NULL}.
#' @param dotColor Character vector of color codes. When \code{dotColorBy} is
#' \code{NULL}, use one or as many colors as the number of cells. If
#' \code{dotColorBy} is categorical, specify as many colors
#' as the number of categories in \code{dotColorBy} or ggplot2 categorical
#' color palette is used by default. If \code{dotColorBy} is continuous, specify
#' together with \code{breaks} argument.
#' @param palette Color palette to use when \code{dotColorBy} is given. Default
#' \code{"D"} (viridis) for continuous value and ggplot2 default for categorical
#' value. See detail for alternatives.
#' @param direction Sets the order of colors in the scale. Default \code{1}
#' orders as palette default. If \code{-1}, the order of colors is reversed.
#' @param breaks Number of breaks for continuous color scale passed to
#' non-interactive "plot3D::scatter3D" call. Default \code{NULL}.
#' @param legendTitle Title on the legend/colorbar. Default \code{NULL} uses
#' \code{"cluster"} if \code{dotColorBy} is missing (default); user-end variable
#' expression if \code{dotColorBy} is directly specified from
#' plotQuaternary.default method; variable name if \code{dotColorBy} is
#' specified from Seurat or SingleCellExperiment method.
#' @param returnData Logical. Whether to return similarity and aggregated
#' velocity data if applicable instead of generating plot. Default \code{FALSE}.
#' @inheritDotParams plotTernary.simMat title nGrid radius dotSize dotShuffle labelColors vertexLabelSize vertexLabelDrift axisBreak axisTextShow axisTextSize axisTextDrift gridLineAlpha arrowLinewidth arrowAngle arrowLen titleSize equilateral margin interactive
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
        dotColorBy = NULL,
        dotColor = NULL,
        palette = "D",
        direction = 1,
        breaks = NULL,
        legendTitle = NULL,
        returnData = FALSE,
        ...
) {
    method <- match.arg(method)
    vcheck <- .checkVertex(x, clusterVar, vertices, n = 3)
    vertClust <- vcheck[[1]]
    vertices <- vcheck[[2]]

    colorArg <- .checkDotColor(
        n = ncol(x), cluster = clusterVar, colorBy = dotColorBy,
        colors = dotColor, palette = palette, direction = direction,
        breaks = breaks,
        legendTitle = legendTitle
    )
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

    if (isTRUE(returnData)) {
        return(list(
            sim = simMat,
            velo = veloMat,
            colorArg = colorArg
        ))
    }

    if (is.null(byCluster)) {
        return(plotTernary(
            x = simMat,
            veloMat = veloMat,
            colorArg = colorArg,
            legendTitle = legendTitle,
            ...
        ))
    } else {
        if (identical(byCluster, "all")) {
            pl <- lapply(levels(clusterVar), function(clust) {
                plotTernary(simMat[clusterVar == clust,],
                            veloMat = veloMat[clusterVar == clust,],
                            title = clust,
                            colorArg = colorArg[clusterVar == clust], ...)
            })
            pl$allCells <- plotTernary(simMat,
                                       veloMat = veloMat,
                                       colorArg = colorArg,
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
                            colorArg = colorArg[clusterVar == clust], ...)
            })
            names(pl) <- byCluster
        }
        return(pl)
    }
}

#' Create ternary simplex plot with Seurat objects
#' @inherit plotTernary description return
#' @param x A Seurat object
#' @param layer Layer in the specified assay to use. Default \code{"counts"}.
#' @param assay The assay to get data from. Default \code{NULL} uses the default
#' assay.
#' @param clusterVar A variable name in meta.data (\code{x[[]]}). Default
#' \code{NULL} uses \code{Idents(x)}.
#' @inheritParams plotTernary.default
#' @inheritDotParams plotTernary.default vertices features veloGraph byCluster method force sigma scale dotColor palette direction breaks legendTitle returnData
#' @inheritDotParams plotTernary.simMat title nGrid radius dotSize dotShuffle labelColors vertexLabelSize vertexLabelDrift axisBreak axisTextShow axisTextSize axisTextDrift gridLineAlpha arrowLinewidth arrowAngle arrowLen titleSize equilateral margin interactive
#'
#' @export
#' @family plotTernary
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
        dotColorBy = NULL,
        legendTitle = NULL,
        ...
) {
    values <- .getSeuratData(x, layer = layer, assay = assay,
                             clusterVar = clusterVar)
    if (missing(processed)) {
        processed <- layer != "counts"
    }
    if (!is.null(dotColorBy)) {
        if (is.character(dotColorBy) && length(dotColorBy) == 1) {
            if (!dotColorBy %in% names(x[[]])) {
                cli::cli_abort(
                    c(x = "Variable {.val {dotColorBy}} does not exist in meta.data {.code x[[]]}.",
                      i = "Available ones are: {.val {names(x[[]])}}.")
                )
            }
            legendTitle <- legendTitle %||% dotColorBy
            dotColorBy <- x[[dotColorBy, drop = TRUE]]
        }
        legendTitle <- legendTitle %||% deparse(substitute(dotColorBy))
    }
    plotTernary(values[[1]], clusterVar = values[[2]], processed = processed,
                dotColorBy = dotColorBy, legendTitle = legendTitle,
                ...)
}

#' Create ternary simplex plot with SingleCellExperiment objects
#' @inherit plotTernary description return
#' @param x A SingleCellExperiment object.
#' @param assay.type Assay to use for calculating the similarity. Default
#' \code{"counts"}.
#' @param clusterVar A variable name in \code{colData(x)}. Default \code{NULL}
#' uses \code{colLabels(x)}.
#' @inheritParams plotTernary.default
#' @inheritDotParams plotTernary.default vertices features veloGraph byCluster method force sigma scale dotColor palette direction breaks legendTitle returnData
#' @inheritDotParams plotTernary.simMat title nGrid radius dotSize dotShuffle labelColors vertexLabelSize vertexLabelDrift axisBreak axisTextShow axisTextSize axisTextDrift gridLineAlpha arrowLinewidth arrowAngle arrowLen titleSize equilateral margin interactive
#' @family plotTernary
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
        dotColorBy = NULL,
        legendTitle = NULL,
        ...
) {
    values <- .getSCEData(x, assay.type = assay.type, clusterVar = clusterVar)
    if (missing(processed)) {
        processed <- assay.type != "counts"
    }
    if (!is.null(dotColorBy)) {
        if (is.character(dotColorBy) && length(dotColorBy) == 1) {
            if (!dotColorBy %in% names(SummarizedExperiment::colData(x))) {
                cli::cli_abort(
                    c(x = "Variable {.val {dotColorBy}} does not exist in {.code colData(x)}.",
                      i = "Available ones are: {.val {names(x@colData)}}.")
                )
            }
            legendTitle <- legendTitle %||% dotColorBy
            dotColorBy <- x[[dotColorBy]]
        }
        legendTitle <- legendTitle %||% deparse(substitute(dotColorBy))
    }
    plotTernary(values[[1]], clusterVar = values[[2]], processed = processed,
                dotColorBy = dotColorBy, legendTitle = legendTitle,
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

#' Create quaternary simplex plot with pre-calculated similarity matrix
#' @param x simMat object, n cells by 4 vertices, each row summing to 1.
#' @param title Title text of the plot. Default \code{NULL}.
#' @param veloMat Aggregated velocity matrix. Output of \code{aggrVeloGraph}.
#' @param nGrid Number of grids along the bottom side of the equilateral
#' triangle. Default \code{10}.
#' @param radius Arrow length of unit velocity. Lower this when arrows point
#' outside of the coordinate. Default \code{0.1}.
#' @param dotSize Dot aesthetics passed to \code{geom_point}. Default \code{0.6}
#' when not interactive, \code{4} when interactive.
#' @param dotShuffle Whether to shuffle the order of dots being added to the
#' plot, useful when categorical colors are used and mixing of categories is
#' expected. Default \code{NULL} does shuffle when \code{dotColorBy} given is
#' categorical and does not otherwise.
#' @param colorArg A "colorArg" object, internally prepared by
#' \code{\link{plotQuaternary.default}}. Default \code{NULL}.
#' @param labelColors Colors of the axis lines and vertex labels.
#' Default \code{c("#3B4992FF", "#EE0000FF", "#008B45FF")} (blue, red and green)
#' @param vertexLabelSize Size of vertex labels. Default \code{6} when not
#' interactive, \code{16} when interactive.
#' @param vertexLabelDrift Position adjustment of the vertex labels, only
#' applied to non-interactive view. Default \code{0.03}.
#' @param axisBreak Number of breaks to be labeled along axis. Default
#' \code{5}.
#' @param axisTextShow Logical, whether to show axis text. Default \code{TRUE}.
#' @param axisTextSize Size of text along each axis break. Default \code{4} for
#' non-interactive view, \code{12} for interactive view.
#' @param axisTextDrift Position adjustment of the axis text, only applied to
#' non-interactive view. Default \code{0.01}.
#' @param gridLineAlpha Transparency of background grid lines. Default
#' \code{0.6}.
#' @param arrowLinewidth Line width of the velocity arrows. Default \code{0.25}
#' for non-interactive view, \code{2} for interactive view.
#' @param arrowAngle Controls the angle of the arrowhead, only applied to
#' non-interactive view. Default \code{20}.
#' @param arrowLen Control length in centimetre from arrow tip to arrow tail,
#' only applied to non-interactive view. Default \code{0.2}.
#' @param titleSize Size of title text. Default \code{14} for non-interactive
#' view, \code{20} for interactive view.
#' @param equilateral Logical, whether to always display the triangle as
#' equilateral. Default \code{TRUE}.
#' @param margin Margin allowed around of the triangle plotting region when
#' \code{equilateral = TRUE}
#' @param interactive Logical. Whether to display plotly interactive view.
#' Default \code{FALSE}.
#' @param ... Not used
#' @method plotTernary simMat
#' @export
plotTernary.simMat <- function(
        x,
        title = NULL,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.1,
        dotSize = NULL,
        dotShuffle = NULL,
        colorArg = NULL,
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF"),
        vertexLabelSize = NULL,
        vertexLabelDrift = 0.03,
        axisBreak = 5,
        axisTextShow = TRUE,
        axisTextSize = NULL,
        axisTextDrift = 0.01,
        gridLineAlpha = 0.6,
        arrowLinewidth = NULL,
        arrowAngle = 20,
        arrowLen = 0.2,
        titleSize = NULL,
        equilateral = TRUE,
        margin = 0.1,
        interactive = FALSE,
        ...
) {
    if (isFALSE(interactive)) {
        .plotTernGG(
            x, title = title, veloMat = veloMat, nGrid = nGrid,
            radius = radius, dotSize = dotSize, colorArg = colorArg,
            dotShuffle = dotShuffle,
            labelColors = labelColors, vertexLabelSize = vertexLabelSize,
            vertexLabelDrift = vertexLabelDrift, axisBreak = axisBreak,
            axisTextShow = axisTextShow, axisTextSize = axisTextSize,
            axisTextDrift = axisTextDrift, gridLineAlpha = gridLineAlpha,
            arrowLinewidth = arrowLinewidth, arrowAngle = arrowAngle,
            arrowLen = arrowLen, titleSize = titleSize,
            equilateral = equilateral, margin = margin
        )
    } else {
        .plotTernPlotly(
            x, title = title, veloMat = veloMat, nGrid = nGrid,
            radius = radius, dotSize = dotSize, colorArg = colorArg,
            dotShuffle = dotShuffle,
            labelColors = labelColors, vertexLabelSize = vertexLabelSize,
            vertexLabelDrift = vertexLabelDrift, axisBreak = axisBreak,
            axisTextShow = axisTextShow, axisTextSize = axisTextSize,
            axisTextDrift = axisTextDrift, gridLineAlpha = gridLineAlpha,
            arrowLinewidth = arrowLinewidth, arrowAngle = arrowAngle,
            arrowLen = arrowLen, titleSize = titleSize,
            equilateral = equilateral, margin = margin
        )
    }
}


.plotTernGG <- function(
        x,
        title = NULL,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.1,
        dotSize = 0.6,
        colorArg = NULL,
        dotShuffle = NULL,
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF"),
        vertexLabelSize = 5,
        vertexLabelDrift = 0.03,
        axisBreak = 5,
        axisTextShow = TRUE,
        axisTextSize = 4,
        axisTextDrift = 0.01,
        gridLineAlpha = 0.6,
        arrowLinewidth = 0.25,
        arrowAngle = 20,
        arrowLen = 0.2,
        titleSize = 14,
        equilateral = TRUE,
        margin = 0.1
) {
    axisTextSize <- axisTextSize %||% 4
    vertexLabelSize <- vertexLabelSize %||% 6
    titleSize <- titleSize %||% 14
    dotSize <- dotSize %||% 0.6
    arrowLinewidth <- arrowLinewidth %||% 0.25
    dotShuffle <- dotShuffle %||% colorArg$type == "categorical"
    df <- as.data.frame(as.matrix(x[,1:3]) %*% triangle)
    legendTitle <- colorArg$palette$legendTitle
    if (colorArg$type == "continuous" ||
        colorArg$type == "categorical") {
        df[[legendTitle]] <- colorArg$colorBy
    }
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    triangle <- as.data.frame(triangle)
    if (isTRUE(dotShuffle)) df <- df[sample(nrow(df)),]
    if (colorArg$type == "continuous") {
        ternCoord <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[[legendTitle]])) +
            geom_point(size = dotSize) +
            switch(colorArg$palette$dep,
                viridis = viridis::scale_color_viridis(
                    option = colorArg$palette$palette,
                    direction = colorArg$palette$direction
                    ),
                RcolorBrewer = scale_color_distiller(
                    palette = colorArg$palette$palette,
                    direction = colorArg$palette$direction,
                    type = "seq",
                    guide = "colourbar"
                    )
            )
        # if (colorArg$palette)
    } else if (colorArg$type == "categorical") {
        ternCoord <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[[legendTitle]])) +
            geom_point(size = dotSize) +
            scale_color_manual(values = colorArg$palette$colors, drop = FALSE) +
            # Dot size in the legend area should be larger
            guides(color = guide_legend(override.aes = list(size = 4)))
    } else if (colorArg$type == "single") {
        ternCoord <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
            geom_point(size = dotSize, color = colorArg$palette$colors)
    } else if (colorArg$type == "customize") {
        ternCoord <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
            geom_point(size = dotSize, color = colorArg$colors)
    }
    ternCoord <- ternCoord +
        annotate("segment", x = triangle$x[1], xend = triangle$x[2],
                 y = triangle$y[1], yend = triangle$y[2],
                 colour = labelColors[1]) +
        annotate("segment", x = triangle$x[2], xend = triangle$x[3],
                 y = triangle$y[2], yend = triangle$y[3],
                 colour = labelColors[2]) +
        annotate("segment", x = triangle$x[3], xend = triangle$x[1],
                 y = triangle$y[3], yend = triangle$y[1],
                 colour = labelColors[3]) +
        annotate("text", x = -0.06 - vertexLabelDrift, hjust = 0,
                 y = 0.0 - vertexLabelDrift, label = colnames(x)[1],
                 colour = labelColors[1], size = vertexLabelSize,
                 fontface = "bold") +
        annotate("text", x = 0.5, y = 3^0.5/2 + 0.01 + vertexLabelDrift,
                 hjust = 0.5,
                 label = colnames(x)[2], colour = labelColors[2],
                 size = vertexLabelSize, fontface = "bold") +
        annotate("text", x = 1.06 + vertexLabelDrift, hjust = 1,
                 y = 0.0 - vertexLabelDrift, label = colnames(x)[3],
                 colour = labelColors[3], size = vertexLabelSize,
                 fontface = "bold") +
        labs(title = title) +
        theme_void() +
        theme(plot.title = element_text(face = "bold", size = titleSize,
                                        hjust = 0.5),
              legend.title = element_text(face = "bold", size = titleSize - 2)
        )

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
                         hjust = 1, vjust = 0,
                         color = labelColors[1], size = axisTextSize) +

                annotate("text", label = perc,
                         x = triangle$x[3] + vecRTBreak[1]*i + axisTextDrift,
                         y = triangle$y[3] + vecRTBreak[2]*i + axisTextDrift,
                         hjust = 0, vjust = 0,
                         color = labelColors[2], size = axisTextSize) +

                annotate("text", label = perc,
                         x = triangle$x[1] + vecLRBreak[1]*i,
                         y = triangle$y[1] + vecLRBreak[2]*i - axisTextDrift,
                         hjust = 0.5, vjust = 1,
                         color = labelColors[3], size = axisTextSize)
        }
    }
    # p <- ternCoord + geom_point(color = dotColor, size = dotSize)

    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = x, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            subcoords <- arrowCoords[[i]]
            ternCoord <- ternCoord +
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
            ternCoord <- ternCoord +
                coord_fixed(xlim = c(0 - margin, 1 + margin),
                            ylim = c(0 - margin, 3^0.5/2 + margin))
        })
    }
    return(ternCoord)
}



# .plotTernPlotly(
#     x, title = title, veloMat = veloMat, nGrid = nGrid,
#     radius = radius, dotSize = dotSize, colorArg = colorArg,
#     labelColors = labelColors, vertexLabelSize = vertexLabelSize,
#     vertexLabelDrift = vertexLabelDrift, axisBreak = axisBreak,
#     axisTextShow = axisTextShow, axisTextSize = axisTextSize,
#     axisTextDrift = axisTextDrift, gridLineAlpha = gridLineAlpha,
#     arrowLinewidth = arrowLinewidth, arrowAngle = arrowAngle,
#     arrowLen = arrowLen, titleSize = titleSize,
#     equilateral = equilateral, margin = margin
# )

.plotTernPlotly <- function(
        x = x,
        title = NULL,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.2,
        dotSize = 4,
        colorArg = NULL,
        dotShuffle = NULL,
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF"),
        vertexLabelSize = 16,
        vertexLabelDrift = 0.03,
        axisBreak = 5,
        axisTextShow = TRUE,
        axisTextSize = 12,
        axisTextDrift = 0.02,
        gridLineAlpha = 0.6,
        arrowLinewidth = 2,
        arrowAngle = 20,
        arrowLen = 0.2,
        titleSize = 20,
        equilateral = TRUE,
        margin = 0.1
) {
    # Since plotly supports ternary coordinate, no need to convert to cartesian
    axisTextSize <- axisTextSize %||% 12
    vertexLabelSize <- vertexLabelSize %||% 16
    titleSize <- titleSize %||% 20
    dotSize <- dotSize %||% 4
    arrowLinewidth <- arrowLinewidth %||% 2
    df <- as.data.frame(as.matrix(x[,1:3]))
    hoverDF <- cbind(ID = rownames(x), df)
    legendTitle <- colorArg$palette$legendTitle
    if (colorArg$type == "categorical") {
        df[[legendTitle]] <- colorArg$colorBy
        fig <- plotly::plot_ly(
            data = df, type = "scatterternary", mode = "markers",
            a = I(df[[1]]), b = I(df[[2]]), c = I(df[[3]]),
            hovertemplate = .plotly_formatCustomData(hoverDF),
            color = stats::formula(paste0("~", legendTitle)),
            colors = colorArg$palette$colors,
            marker = list(
                size = dotSize
            )
        )
    } else if (colorArg$type == "continuous") {
        df[[legendTitle]] <- colorArg$colorBy
        colorScale <- .plotlyColorscale(colorArg$palette$colors)
        fig <- plotly::plot_ly(
            data = df, type = "scatterternary", mode = "markers",
            a = I(df[[1]]), b = I(df[[2]]), c = I(df[[3]]),
            hovertemplate = .plotly_formatCustomData(hoverDF),
            marker = list(
                color = stats::formula(paste0("~", legendTitle)),
                autocolorscale = FALSE,
                colorscale = colorScale,
                cmin = min(colorArg$colorBy),
                cmax = max(colorArg$colorBy),
                size = dotSize,
                showscale = TRUE,
                colorbar = list(
                    title = .htmlbold(legendTitle),
                    titleside = "top",
                    titlefont = list(size = 14),
                    tickfont = list(size = 12),
                    len = 0.6
                )
            )
        )
    } else {
        fig <- plotly::plot_ly(
            data = df, type = "scatterternary", mode = "markers",
            a = I(df[[1]]), b = I(df[[2]]), c = I(df[[3]]),
            hovertemplate = .plotly_formatCustomData(hoverDF),
            marker = list(
                color = colorArg$colors,
                size = dotSize
            )
        )
    }

    subplotList <- list(fig)
    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = x, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            legendGroup <- paste0("velo_", i)
            subcoords <- arrowCoords[[i]]
            potential <- sqrt(rowSums((subcoords[, c(3,4)] - subcoords[, c(1,2)])^2))/radius
            subcoordsStartBary <- cart2bary(triangle, subcoords[,c(1,2)])
            subcoordsEndBary <- cart2bary(triangle, subcoords[,c(3,4)])
            for (j in seq(nrow(subcoords))) {
                arrowBary <- rbind(subcoordsStartBary[j,], subcoordsEndBary[j,])
                arrowBary <- as.data.frame(arrowBary)
                subplotList[[paste0("arrow_", i, "_", j)]] <- plotly::plot_ly(
                    data = arrowBary,
                    type = "scatterternary",
                    mode = "lines+markers",
                    a = I(arrowBary[[1]]), b = I(arrowBary[[2]]),
                    c = I(arrowBary[[3]]),
                    # Hack to show only one legend item
                    showlegend = j == 1,
                    name = colnames(x)[i],
                    color = I(labelColors[i]),
                    line = list(width = arrowLinewidth),
                    hovertemplate = sprintf(
                        "Potential to %s: %.3f", colnames(x)[i], potential[j]
                    ),
                    # angleref = "previous" is the magic that makes the arrow
                    # pointing from start to end
                    # Must have later plotly javascript version than R plotly
                    # package default
                    marker = list(
                        size = 8,
                        symbol = "arrow",
                        angleref = "previous"
                    ),
                    legendgroup = legendGroup,
                    legendgrouptitle = if (i == 1) list(text = .htmlbold("Velocity to")) else NULL
                )
            }
        }
    }

    axisArg <- function(i) {
        list(
            title = list(
                text = .htmlbold(colnames(x)[i]),
                font = list(size = vertexLabelSize)
            ),
            color = labelColors[i],
            gridcolor = .colAlpha(labelColors[i], gridLineAlpha),
            linewidth = 1,
            ticklen = 5,
            tickcolor = labelColors[i],
            tickmode = 'linear',
            tick0 = 0,
            dtick = 1 / axisBreak,
            showline = TRUE,
            showgrid = TRUE,
            showticklabels = axisTextShow,
            tickfont = list(size = axisTextSize),
            gridcolor = labelColors[i],
            gridwidth = 1,
            gridalpha = gridLineAlpha
        )
    }
    fig <- plotly::subplot(subplotList)
    fig <- plotly::layout(
        fig,
        ternary = list(
            aaxis = axisArg(1),
            baxis = axisArg(2),
            caxis = axisArg(3),
            domain = list(x = c(0, 1), y = c(0, 0.9))
        ),
        title = list(
            text = .htmlbold(title),
            font = list(size = titleSize),
            xanchor = "center",
            xref = "paper",
            y = 1,
            pad = list(t = 10, b = 10)
        ),
        legend = list(
            title = list(text = .htmlbold(legendTitle)),
            itemsizing = "constant",
            y = 0.5
        ),
        legend2 = list(
            title = list(text = .htmlbold("Velocity")),
            itemsizing = "constant"
        )
    )
    # HACK for avoiding unreasonable warnings
    fig$x$layout <- fig$x$layout[grep('NA', names(fig$x$layout), invert = TRUE)]

    # LAST HACK FOR ALLOWING THE ARROW HEAD TO ROTATE with "angleref = 'previous'"
    # BLAME PLOTLY.R TEAM FOR NOT UPDATING THIS
    # THANKS TO Kat for https://stackoverflow.com/questions/77348179/how-can-i-set-arrow-markers-to-rotate-in-plotly-for-r
    dep <- list(
        name = "plotly-latest",
        version = "2.21.1",
        src = list(
            href = "https://cdn.plot.ly"
        ),
        meta = NULL,
        script = "plotly-2.21.0.min.js",
        stylesheet = NULL, head = NULL, attachment = NULL, package = NULL,
        all_files = TRUE
    )
    class(dep) <- "html_dependency"
    fig$dependencies <- append(fig$dependencies, list(dep))
    return(fig)
}
