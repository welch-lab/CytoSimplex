#' Create quaternary plots
#' @description
#' Create quaternary plots that show similarity between single cells and
#' selected terminals in a tetrahedron space.
#' @param object An object
#' @param ... Arguments passed to other methods.
#' @rdname plotQuaternary
#' @export plotQuaternary
#' @return For distMatrix method, a plist (plot3D package product) object. For
#' other methods, a plist object when \code{splitCluster = FALSE}, or a list of
#' plist objects when \code{splitCluster = TRUE}. A plist object can be viewed
#' with \code{print()}, \code{show()} or a direct run of the object variable in
#' console.
#' @examples
#' rnaLog <- colNormalize(rnaRaw, 1e4, TRUE)
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' plotQuaternary(rnaLog[gene, ], rnaCluster, c("RE", "OS", "CH", "ORT"))
plotQuaternary <- function(object, ...) {
    UseMethod('plotQuaternary', object)
}

#' @rdname plotQuaternary
#' @param dotSize,dotColor Dot aesthetics. Default \code{0.6} and
#' \code{"grey60"}.
#' @param labelColors Colors of the vertex labels. Default
#' \code{c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF")} (blue, red,
#' green and purple).
#' @param distMethod Method name used for calculating the dist.matrix object.
#' @param title Title text of the plot. Default \code{NULL}.
#' @param theta,phi Numeric scalar. The angles defining the viewing direction.
#' \code{theta} gives the azimuthal direction and \code{phi} the colatitude.
#' Default \code{-40} and \code{-10}.
#' @export
#' @method plotQuaternary simMat
plotQuaternary.simMat <- function(
        object,
        dotSize = 0.6,
        dotColor = "grey60",
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
        distMethod = attributes(object)$method,
        title = NULL,
        theta = -40,
        phi = -10,
        ...
) {
    if (ncol(object) != 5) {
        stop("`simMat` object must have five columns for quaternary plot, ",
             "where the first four are for vertices and the last for cluster ",
             "assignment.")
    }
    distNorm <- t(apply(object[,seq(4)], 1, .normalize))

    # Compute tetrahedron coordinates according to
    # https://mathoverflow.net/a/184585
    tetra <- qr.Q(qr(matrix(1, nrow = 4)), complete = TRUE)[,-1]

    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    df3D <- distNorm %*% tetra

    # Plot data
    grDevices::pdf(nullfile())
    scatter3D(df3D[,1], df3D[,2], df3D[,3],
              xlim = range(tetra[,1]),
              ylim = range(tetra[,2]),
              zlim = range(tetra[,3]), alpha = 0.8,
              col = dotColor, cex = dotSize/2, pch = 16, d = 3,
              colkey = list(plot = FALSE),
              box = FALSE, theta = theta, phi = phi, plot = FALSE)
    lines3D(tetra[c(1,2,3,4,1,3,1,2,4),1],
            tetra[c(1,2,3,4,1,3,1,2,4),2],
            tetra[c(1,2,3,4,1,3,1,2,4),3],
            col = "grey", add = TRUE, plot = FALSE)
    text3D(tetra[,1], tetra[,2], tetra[,3],
           colnames(object)[seq(4)], col = labelColors,
           add = TRUE, plot = FALSE)
    grDevices::dev.off()
    getplist()
}

setOldClass("plist")

#' @rdname show-plot3D
#' @title Show plist object produced with plot3D package
#' @param object,x plist object
#' @param ... Not used.
#' @export
setMethod("show", "plist", function(object) {
    print(object)
}
)

#' @rdname show-plot3D
#' @method print plist
#' @export
print.plist <- function(x, ...) {
    graphics::par(mar = c(1, 0, 1, 0), oma = c(0, 0, 0, 0),
                  mgp = c(0, 0, 0), xpd = NA)
    plot(x, ...)
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
#' @param force Whether to force calculate the distance when more then 500
#' features are detected, which is generally not recommended. Default
#' \code{FALSE}.
#' @param sigma Gaussian kernel parameter that controls the effect of variance.
#' Only effective when using a distance metric (i.e. \code{method} is
#' \code{"euclidian"} or \code{"cosine"}). Larger value tighten the dot
#' spreading on figure. Default \code{4}.
#' @param scale Whether to min-max scale the distance matrix by clusters.
#' Default \code{TRUE}.
#' @param splitCluster Logical, whether to return a list of plots where each
#' contains cells belonging to only one cluster. Default \code{FALSE}.
#' @param clusterTitle If \code{splitCluster = TRUE}, whether to title each
#' subplot by the name of each cluster. Default \code{TRUE}.
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary default
plotQuaternary.default <- function(
        object,
        clusterVar,
        vertices,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        force = FALSE,
        sigma = 4,
        scale = TRUE,
        splitCluster = FALSE,
        clusterTitle = TRUE,
        dotColor = "grey60",
        ...
) {
    method <- match.arg(method)
    vertices <- unique(vertices)
    if (length(vertices) < 4) {
        stop("Must specify 4 different vertices.")
    } else if (length(vertices) > 4) {
        vertices <- vertices[seq(4)]
        warning("More than 4 vertices specified for quaternary plot. ",
                "Using the first 4.", immediate. = TRUE)
    }
    if (!all(vertices %in% clusterVar)) {
        stop("Specified vertex clusters are not all found in the cluster ",
             "variable")
    }
    if (length(dotColor) == 1) dotColor <- rep(dotColor, ncol(object))
    if (length(dotColor) != ncol(object)) {
        stop("`dotColor` need to be either 1 scalar or match the number of ",
             "samples in `object`.")
    }
    distMat <- calcDist2(object, clusterVar = clusterVar,
                         vertices = vertices, method = method,
                         scale = scale, force = force, sigma = sigma)
    if (isFALSE(splitCluster)) plotQuaternary(object = distMat,
                                              dotColor = dotColor, ...)
    else {
        if (isTRUE(clusterTitle)) {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotQuaternary(distMat[distMat$Label == clust,], title = clust,
                               dotColor = dotColor[clusterVar == clust], ...)
            })
        } else {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotQuaternary(distMat[distMat$Label == clust,],
                               dotColor = dotColor[clusterVar == clust], ...)
            })
        }
        names(plotList) <- levels(clusterVar)
        return(plotList)
    }
}

#' #' @rdname plotQuaternary
#' #' @param useDatasets For liger method, select datasets where the distance
#' #' calculation should only be limited within this range. Default \code{NULL}
#' #' uses all datasets.
#' #' @param features For container object methods. Valid row subsetting index that
#' #' selects features. Default \code{NULL}.
#' #' @export
#' #' @method plotQuaternary liger
#' plotQuaternary.liger <- function(
#'         object,
#'         clusterVar,
#'         features = NULL,
#'         useDatasets = NULL,
#'         ...
#' ) {
#'     values <- .ligerPrepare(object, clusterVar, features, useDatasets)
#'     plotQuaternary(values[[1]], clusterVar = values[[2]], ...)
#' }
