#' Create binary plots
#' @description
#' Create binary plots that show similarity between single cells and two
#' selected terminals in a 2D space. The two vertices are placed at the left and
#' right of a 2D plot where x-axis measures the similarity. Y-axis is jittered
#' for a clear view. A density (histogram) curve is added for indicating the
#' distribution.
#' @param object An object
#' @param ... Arguments passed to other methods.
#' @rdname plotBinary
#' @export plotBinary
#' @return For 'simMat' method, a ggplot object. For other methods, a ggplot
#' object when \code{splitCluster = FALSE}, or a list of ggplot objects when
#' \code{splitCluster = TRUE}.
#' @examples
#' rnaNorm <- colNormalize(rnaRaw)
#' gene <- selectTopFeatures(rnaNorm, rnaCluster, c("RE", "OS"))
#' rnaLog <- colNormalize(rnaRaw, 1e4, TRUE)
#' plotBinary(rnaLog[gene, ], rnaCluster, c("RE", "OS"))
plotBinary <- function(object, ...) {
    UseMethod('plotBinary', object)
}

#' @rdname plotBinary
#' @param dotSize,dotColor Dot aesthetics passed to
#' \code{\link[ggplot2]{geom_point}}. Default \code{0.6} and \code{"grey60"}.
#' @param densLinewidth Density plot line aesthetic. Default \code{0.8}.
#' @param labelColors Color of the axis lines and vertex labels. Default
#' \code{c("#3B4992FF", "#EE0000FF")} (blue and red).
#' @param distMethod Method name used for calculating the dist.matrix object.
#' @param title Title text of the plot. Default \code{NULL}.
#' @export
#' @method plotBinary simMat
plotBinary.simMat <- function(
        object,
        dotSize = 0.6,
        dotColor = "grey60",
        densLinewidth = 0.8,
        labelColors = c("#3B4992FF", "#EE0000FF"),
        distMethod = attributes(object)$method,
        title = NULL,
        ...
) {
    topLAB <- colnames(object)[1]
    bottomLAB <- colnames(object)[2]
    object$Y <- stats::runif(nrow(object))

    object[,1] <- 100 * object[,1]
    object[,2] <- 100 * object[,2]
    sumProp <- 100

    # top <- bottom <- NULL
    # NEVER REMOVE "ggplot2::", the imported namespace has problem with
    # rlang ".data" pronoun
    p <- ggplot(object, ggplot2::aes(x = .data[[bottomLAB]],
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

#' @param clusterVar A vector/factor assigning the cluster variable to each
#' column of the matrix object. For "Seurat" method, leave \code{NULL} for using
#' default "Idents", or can also be a variable in \code{meta.data} slot. For
#' "SingleCellExperiment" method, leave \code{NULL} for using "colLabels", or
#' can be a variable in \code{colData} slot.
#' @param vertices Vector of two unique cluster names that will be used for
#' plotting. Or a named list that groups clusters as two terminal vertices.
#' There must not be any overlap between groups.
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
#' @rdname plotBinary
#' @export
#' @method plotBinary default
plotBinary.default <- function(
        object,
        clusterVar,
        vertices,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        force = FALSE,
        sigma = 0.08,
        scale = TRUE,
        splitCluster = FALSE,
        clusterTitle = TRUE,
        dotColor = "grey60",
        ...
) {
    method <- match.arg(method)
    vcheck <- .checkVertex(object, clusterVar, vertices, n = 2)
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
    if (isFALSE(splitCluster)) plotBinary(object = simMat,
                                          dotColor = dotColor,
                                          ...)
    else {
        if (isTRUE(clusterTitle)) {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotBinary(simMat[clusterVar == clust,], title = clust,
                           dotColor = dotColor[clusterVar == clust],
                           ...)
            })
        } else {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotBinary(simMat[clusterVar == clust,],
                           dotColor = dotColor[clusterVar == clust], ...)
            })
        }
        names(plotList) <- levels(clusterVar)
        return(plotList)
    }
}

#' @param features For container object methods. Valid row subsetting index that
#' selects features. Default \code{NULL} uses all available features.
#' @param slot For Seurat method, choose from \code{"data"},
#' \code{"scale.data"} or \code{"counts"}. Default \code{"data"}.
#' @param assay For Seurat method, the specific assay to get data from. Default
#' \code{NULL} to the default assay.
#' @rdname plotBinary
#' @export
#' @method plotBinary Seurat
#' @examples
#'
#' # Seurat example
#' if (FALSE) {
#'     srt <- CreateSeuratObject(rnaRaw)
#'     Idents(srt) <- rnaCluster
#'     srt <- colNormalize(srt)
#'     gene <- selectTopFeatures(srt, vertices = c("OS", "RE"))
#'     srt <- colNormalize(srt, scaleFactor = 1e4, log = TRUE)
#'     plotBinary(srt, features = gene, vertices = c("OS", "RE"))
#' }
plotBinary.Seurat <- function(
        object,
        features = NULL,
        slot = "data",
        assay = NULL,
        ...
) {
    values <- .getSeuratData(object, features = features,
                             slot = slot, assay = assay)
    plotBinary(values[[1]], values[[2]], ...)
}

#' @param subset.row For "SingleCellExperiment" methods. Valid row subsetting
#' index that selects features. Default \code{NULL} uses all available features.
#' @param assay.type For "SingleCellExperiment" methods. Which assay to use for
#' calculating the similarity. Default \code{"logcounts"}.
#' @rdname plotBinary
#' @export
#' @method plotBinary SingleCellExperiment
#' @examples
#'
#' # SingleCellExperiment example
#' if (FALSE) {
#'     library(SingleCellExperiment)
#'     sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#'     colLabels(sce) <- rnaCluster
#'     sce <- colNormalize(sce)
#'     gene <- selectTopFeatures(sce, vertices = c("OS", "RE"))
#'     sce <- colNormalize(sce, scaleFactor = 1e4, log = TRUE)
#'     plotBinary(sce, subset.row = gene, vertices = c("OS", "RE"))
#' }
plotBinary.SingleCellExperiment <- function(
        object,
        subset.row = NULL,
        assay.type = "logcounts",
        ...
) {
    values <- .getSCEData(object, subset.row = subset.row,
                          assay.type = assay.type)
    plotBinary(values[[1]], values[[2]], ...)
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
#'         object,
#'         clusterVar,
#'         features = NULL,
#'         useDatasets = NULL,
#'         ...
#' ) {
#'     values <- .ligerPrepare(object, clusterVar, features, useDatasets)
#'     plotBinary(values[[1]], clusterVar = values[[2]], ...)
#' }


