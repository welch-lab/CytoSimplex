#' Create binary plots
#' @description
#' Create binary plots that show similarity between single cells and two
#' selected terminals in a barycentric coordinate. The two vertices are placed
#' at the left and right of a 2D plot where x-axis measures the similarity.
#' Y-axis is jittered for a clear view. A density (histogram) curve is added for
#' indicating the distribution.
#'
#' See \code{\link{plotTernary}} manual for more details.
#' @param x Input data. Can be a \code{matrix} or \code{dgCMatrix} object with
#' cells as columns, a \code{Seurat} or \code{SingleCellExperiment} object.
#' "simMat" method takes intermediate values.
#' @param ... Arguments passed to other methods.
#' @rdname plotBinary
#' @export plotBinary
#' @return For 'simMat' method, a ggplot object. For other methods, a ggplot
#' object when \code{splitCluster = FALSE}, or a list of ggplot objects when
#' \code{splitCluster = TRUE}.
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS"))
#' plotBinary(rnaRaw, rnaCluster, c("RE", "OS"), gene)
plotBinary <- function(x, ...) {
    UseMethod('plotBinary', x)
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
#' @param method Similarity calculation method. Default \code{"euclidean"}.
#' Choose from \code{"euclidean"}, \code{"cosine"}, \code{"pearson"},
#' \code{"spearman"}.
#' @param force Whether to force calculate the similarity when more then 500
#' features are detected, which is generally not recommended. Default
#' \code{FALSE}.
#' @param sigma Gaussian kernel parameter that controls the effect of variance.
#' Only effective when using a distance metric (i.e. \code{method} is
#' \code{"euclidian"} or \code{"cosine"}). Larger value tighten the dot
#' spreading on figure. Default \code{0.08}.
#' @param scale Whether to min-max scale the distance matrix by clusters.
#' Default \code{TRUE}.
#' @param returnData Logical. Whether to return similarity data instead of
#' generating plot. Default \code{FALSE}.
#' @rdname plotBinary
#' @export
#' @method plotBinary default
plotBinary.default <- function(
        x,
        clusterVar,
        vertices,
        features = NULL,
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
    vcheck <- .checkVertex(x, clusterVar, vertices, n = 2)
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

    if (isTRUE(returnData)) return(list(sim = simMat))

    if (is.null(byCluster)) {
        return(plotBinary(x = simMat, dotColor = dotColor, ...))
    } else {
        if (identical(byCluster, "all")) {
            pl <- lapply(levels(clusterVar), function(clust) {
                plotBinary(simMat[clusterVar == clust,], title = clust,
                           dotColor = dotColor[clusterVar == clust], ...)
            })
            pl$allCells <- plotBinary(simMat,
                                      dotColor = dotColor,
                                      title = "All cells", ...)
            names(pl) <- c(levels(clusterVar), "allCells")
        } else {
            if (any(!byCluster %in% levels(clusterVar))) {
                stop("`byCluster` must be either a vector of cluster name ",
                     "or \"all\".")
            }
            pl <- lapply(byCluster, function(clust) {
                plotBinary(simMat[clusterVar == clust,], title = clust,
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
#' @rdname plotBinary
#' @export
#' @method plotBinary Seurat
#' @examples
#' \donttest{
#' # Seurat example
#' library(Seurat)
#' srt <- CreateSeuratObject(rnaRaw)
#' Idents(srt) <- rnaCluster
#' gene <- selectTopFeatures(srt, vertices = c("OS", "RE"))
#' plotBinary(srt, features = gene, vertices = c("OS", "RE"))
#' }
plotBinary.Seurat <- function(
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
    plotBinary(x = values[[1]], clusterVar = values[[2]], processed = processed,
               ...)
}

#' @param assay.type For "SingleCellExperiment" methods. Which assay to use for
#' calculating the similarity. Default \code{"counts"}.
#' @rdname plotBinary
#' @export
#' @method plotBinary SingleCellExperiment
#' @examples
#' \donttest{
#' # SingleCellExperiment example
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#' colLabels(sce) <- rnaCluster
#' gene <- selectTopFeatures(sce, vertices = c("OS", "RE"))
#' plotBinary(sce, features = gene, vertices = c("OS", "RE"))
#' }
plotBinary.SingleCellExperiment <- function(
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
    plotBinary(values[[1]], clusterVar = values[[2]], processed = processed,
               ...)
}

#' #' @rdname plotBinary
#' #' @param useDatasets For liger method, select datasets where the distance
#' #' calculation should only be limited within this range. Default \code{NULL}
#' #' uses all datasets.
#' #' @param features For container object methods. Valid row subsetting index that
#' #' selects features. Default \code{NULL}.
#' #' @export
#' #' @method plotBinary liger
#' plotBinary.liger <- function(
        #'         x,
#'         clusterVar,
#'         features = NULL,
#'         useDatasets = NULL,
#'         ...
#' ) {
#'     values <- .ligerPrepare(x, clusterVar, features, useDatasets)
#'     plotBinary(values[[1]], clusterVar = values[[2]], ...)
#' }

#' @rdname plotBinary
#' @param dotSize,dotColor Dot aesthetics passed to
#' \code{\link[ggplot2]{geom_point}}. Default \code{0.6} and \code{"grey60"}.
#' @param densLinewidth Density plot line aesthetic. Default \code{0.8}.
#' @param labelColors Color of the axis lines and vertex labels. Default
#' \code{c("#3B4992FF", "#EE0000FF")} (blue and red).
#' @param title Title text of the plot. Default \code{NULL}.
#' @export
#' @method plotBinary simMat
plotBinary.simMat <- function(
        x,
        dotSize = 0.6,
        dotColor = "grey60",
        densLinewidth = 0.8,
        labelColors = c("#3B4992FF", "#EE0000FF"),
        title = NULL,
        ...
) {
    topLAB <- colnames(x)[1]
    bottomLAB <- colnames(x)[2]
    x$Y <- stats::runif(nrow(x))

    x[,1] <- 100 * x[,1]
    x[,2] <- 100 * x[,2]
    sumProp <- 100

    # top <- bottom <- NULL
    # NEVER REMOVE "ggplot2::", the imported namespace has problem with
    # rlang ".data" pronoun
    p <- ggplot(x, ggplot2::aes(x = .data[[bottomLAB]],
                                y = .data[["Y"]])) +
        geom_point(color = dotColor, size = dotSize, stroke = 0.2) +
        geom_density(
            ggplot2::aes(x = .data[[bottomLAB]],
                         y = .scaleMinMax(after_stat(.data[["density"]]))),
            inherit.aes = FALSE, linewidth = densLinewidth) +
        scale_x_continuous(
            limits = c(0, sumProp),
            breaks = seq(0, sumProp, 0.2 * sumProp),
            sec.axis = sec_axis(trans = function(x) sumProp - x,
                                name = topLAB,
                                breaks = seq(sumProp, 0, -0.2 * sumProp))
        ) +
        ggtitle(title) +
        theme(panel.background = element_rect(fill = NA),
              axis.line.x.bottom = element_line(colour = labelColors[2],
                                                arrow = arrow()),
              axis.line.x.top = element_line(colour = labelColors[1],
                                             arrow = arrow(ends = "first")),
              axis.title.x.bottom = element_text(face = "bold", hjust = 1,
                                                 colour = labelColors[2]),
              axis.title.x.top = element_text(face = "bold", hjust = 0,
                                              colour = labelColors[1]),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              panel.grid.major.x = element_line(colour = "grey80",
                                                linetype = 3))
    return(p)
}
