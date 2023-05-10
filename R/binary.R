#' Create binary plots
#' @description
#' Create binary plots that show similarity between single cells and
#' selected terminals in a 2D space. The two vertices are placed at the left and
#' right of a 2D plot where x-axis shows the similarity. Y-axis is jittered for
#' clear view. A density curve is added for indicating the distribution.
#' @param object An object
#' @param ... Arguments passed to other methods.
#' @rdname plotBinary
#' @export plotBinary
#' @return For distMatrix method, a ggplot object. For other methods, a ggplot
#' object when \code{splitCluster = FALSE}, or a list of ggplot objects when
#' \code{splitCluster = TRUE}.
plotBinary <- function(object, ...) {
    UseMethod('plotBinary', object)
}

#' @rdname plotBinary
#' @param axisPortion Whether to show axis value by proportion, so that the
#' coordinates of each dot sum up to 100. Default \code{TRUE}.
#' @param dotSize,dotColor Dot aesthetics passed to
#' \code{\link[ggplot2]{geom_point}}. Default \code{0.6} and \code{"grey60"}.
#' @param densLinewidth Density plot line aesthetic. Default \code{0.8}.
#' @param labelColors Color of the axis lines and vertex labels. Default
#' \code{c("#3B4992FF", "#EE0000FF")} (blue and red).
#' @param distMethod Method name used for calculating the dist.matrix object.
#' @param title Title text of the plot. Default \code{NULL}.
#' @export
#' @method plotBinary distMatrix
plotBinary.distMatrix <- function(
        object,
        axisPortion = TRUE,
        dotSize = 0.6,
        dotColor = "grey60",
        densLinewidth = 0.8,
        labelColors = c("#3B4992FF", "#EE0000FF"),
        distMethod = attributes(object)$method,
        title = NULL,
        ...
) {
    if (is.null(distMethod) && isFALSE(axisPortion)) {
        stop("Cannot identify distance method for visualizing axis.")
    }
    if (ncol(object) != 3) {
        stop("`distMatrix` object must have three columns for binary plot, ",
             "where the first two are for vertices and the last for cluster ",
             "assignment.")
    }
    topLAB <- colnames(object)[1]
    bottomLAB <- colnames(object)[2]
    object$Y <- stats::runif(nrow(object))
    sumProp <- 1
    if (isTRUE(axisPortion)) {
        object[,1] <- 100 * object[,1]
        object[,2] <- 100 * object[,2]
        sumProp <- 100
    }

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
#' column of the matrix object. Or a character scalar selecting a cell metadata
#' variable from container object.
#' @param vertices A vector of TWO values specifying the clusters as the
#' terminals.
#'
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
#' @rdname plotBinary
#' @export
#' @method plotBinary default
plotBinary.default <- function(
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
    vertices <- unique(vertices)
    if (length(vertices) < 2) {
        stop("Must specify 2 different vertices.")
    } else if (length(vertices) > 2) {
        vertices <- vertices[seq(2)]
        warning("More than 2 vertices specified for binary plot. ",
                "Using the first two.", immediate. = TRUE)
    }
    if (!all(vertices %in% clusterVar)) {
        stop("Specified vertex clusters are not all found in the cluster ",
             "variable")
    }
    distMat <- calcDist(object, clusterVar = clusterVar,
                        vertices = vertices,
                        method = method, normCluster = normCluster,
                        scale = scale)
    if (isFALSE(splitCluster)) plotBinary(object = distMat, ...)
    else {
        if (isTRUE(clusterTitle)) {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotBinary(distMat[distMat$Label == clust,], title = clust,
                           ...)
            })
        } else {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotBinary(distMat[distMat$Label == clust,], ...)
            })
        }
        names(plotList) <- levels(clusterVar)
        return(plotList)
    }
}

#' @rdname plotBinary
#' @param useDatasets For liger method, select datasets where the distance
#' calculation should only be limited within this range. Default \code{NULL}
#' uses all datasets.
#' @param features For container object methods. Valid row subsetting index that
#' selects features. Default \code{NULL}.
#' @export
#' @method plotBinary liger
plotBinary.liger <- function(
        object,
        clusterVar,
        features = NULL,
        useDatasets = NULL,
        ...
) {
    values <- .ligerPrepare(object, clusterVar, features, useDatasets)
    plotBinary(values[[1]], clusterVar = values[[2]], ...)
}

