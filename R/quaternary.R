#' Create quaternary simplex plots
#' @description
#' Create quaternary plots that show similarity between single cells and
#' selected four terminals in a tetrahedron space.
#'
#' A dynamic rotating view in a GIF image file can be created with
#' \code{\link{writeQuaternaryGIF}}. Package \code{magick} must be installed in
#' advance. Linux users may refer to this
#' \href{https://cran.r-project.org/web/packages/magick/vignettes/intro.html#Build_from_source}{installation guide}.
#' @param object An object
#' @param ... Arguments passed to other methods.
#' @rdname plotQuaternary
#' @export plotQuaternary
#' @return For "simMat" method, a "plist" (plot3D package product) object. For
#' other methods, a "plist" object when \code{splitCluster = FALSE}, or a list
#' of "plist" objects when \code{splitCluster = TRUE}. A "plist" object can be
#' viewed with \code{print()}, \code{show()} or a direct run of the object
#' variable name in interactive console.
#' @examples
#' rnaNorm <- colNormalize(rnaRaw)
#' gene <- selectTopFeatures(rnaNorm, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' rnaLog <- colNormalize(rnaRaw, 1e4, TRUE)
#' plotQuaternary(rnaLog[gene, ], rnaCluster, c("RE", "OS", "CH", "ORT"))
plotQuaternary <- function(object, ...) {
    UseMethod('plotQuaternary', object)
}

#' @rdname plotQuaternary
#' @param veloMat Aggregated velocity matrix. Output of \code{aggrVeloGraph}.
#' @param nGrid Number of grids along the x-axis of the tetrahedron
#' triangle. Default \code{10}.
#' @param radius Arrow length of unit velocity. Lower this when arrows point
#' outside of the tetrahedron. Default \code{0.2}.
#' @param dotSize,dotColor Dot aesthetics. Default \code{0.6} and
#' \code{"grey60"}.
#' @param labelColors Colors of the vertex labels. Default
#' \code{c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF")} (blue, red,
#' green and purple).
#' @param arrowLinewidth Arrow aesthetics. Default \code{0.6}.
#' @param title Title text of the plot. Default \code{NULL}.
#' @param theta,phi Numeric scalar. The angles defining the viewing direction.
#' \code{theta} gives the azimuthal direction and \code{phi} the colatitude.
#' Default \code{-40} and \code{-10}.
#' @export
#' @method plotQuaternary simMat
plotQuaternary.simMat <- function(
        object,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.2,
        dotSize = 0.6,
        dotColor = "grey60",
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
        arrowLinewidth = 0.6,
        title = NULL,
        theta = 0,
        phi = 0,
        ...
) {
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    df3D <- as.matrix(object) %*% tetra


    # Plot data
    grDevices::pdf(nullfile())
    scatter3D(df3D[,1], df3D[,2], df3D[,3], main = title,
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
    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = object, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            subcoord <- arrowCoords[[i]]
            arrows3D(subcoord[,1], subcoord[,2], subcoord[,3],
                     subcoord[,4], subcoord[,5], subcoord[,6],
                     angle = 20, lwd = arrowLinewidth, length = 0.1,
                     col = labelColors[i], add = TRUE, plot = FALSE)
        }
    }
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
#' column of the matrix object. For "Seurat" method, leave \code{NULL} for using
#' default "Idents", or can also be a variable in \code{meta.data} slot. For
#' "SingleCellExperiment" method, leave \code{NULL} for using "colLabels", or
#' can be a variable in \code{colData} slot.
#' @param vertices Vector of four unique cluster names that will be used for
#' plotting. Or a named list that groups clusters as four terminal vertices.
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
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary default
plotQuaternary.default <- function(
        object,
        clusterVar,
        vertices,
        veloGraph = NULL,
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
    vcheck <- .checkVertex(object, clusterVar, vertices, n = 4)
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

    if (isFALSE(splitCluster)) plotQuaternary(object = simMat,
                                              veloMat = veloMat,
                                              dotColor = dotColor, ...)
    else {
        if (isTRUE(clusterTitle)) {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotQuaternary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               title = clust,
                               dotColor = dotColor[clusterVar == clust], ...)
            })
            plotList$allCells <- plotQuaternary(simMat,
                                                veloMat = veloMat,
                                                dotColor = dotColor,
                                                title = "All cells", ...)
        } else {
            plotList <- lapply(levels(clusterVar), function(clust) {
                plotQuaternary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               dotColor = dotColor[clusterVar == clust], ...)
            })
            plotList$allCells <- plotQuaternary(simMat,
                                                veloMat = veloMat,
                                                dotColor = dotColor, ...)
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
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary Seurat
#' @examples
#'
#' # Seurat example
#' if (FALSE) {
#'     srt <- CreateSeuratObject(rnaRaw)
#'     Idents(srt) <- rnaCluster
#'     srt <- colNormalize(srt)
#'     gene <- selectTopFeatures(srt, vertices = c("OS", "RE", "CH", "ORT"))
#'     srt <- colNormalize(srt, scaleFactor = 1e4, log = TRUE)
#'     plotQuaternary(srt, features = geneSel, vertices = c("OS", "RE", "CH", "ORT"))
#' }
plotQuaternary.Seurat <- function(
        object,
        features = NULL,
        slot = "data",
        assay = NULL,
        ...
) {
    values <- .getSeuratData(object, features = features,
                             slot = slot, assay = assay)
    plotQuaternary(values[[1]], values[[2]], ...)
}

#' @param subset.row For "SingleCellExperiment" methods. Valid row subsetting
#' index that selects features. Default \code{NULL} uses all available features.
#' @param assay.type For "SingleCellExperiment" methods. Which assay to use for
#' calculating the similarity. Default \code{"logcounts"}.
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary SingleCellExperiment
#' @examples
#'
#' # SingleCellExperiment example
#' if (FALSE) {
#'     library(SingleCellExperiment)
#'     sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#'     colLabels(sce) <- rnaCluster
#'     sce <- colNormalize(sce)
#'     gene <- selectTopFeatures(sce, vertices = c("OS", "RE", "CH", "ORT"))
#'     sce <- colNormalize(sce, scaleFactor = 1e4, log = TRUE)
#'     plotQuaternary(sce, subset.row = gene,
#'                    vertices = c("OS", "RE", "CH", "ORT"))
#' }
plotQuaternary.SingleCellExperiment <- function(
        object,
        subset.row = NULL,
        assay.type = "logcounts",
        ...
) {
    values <- .getSCEData(object, subset.row = subset.row,
                          assay.type = assay.type)
    plotQuaternary(values[[1]], clusterVar = values[[2]], ...)
}

# #' @rdname plotQuaternary
# #' @param useDatasets For liger method, select datasets where the distance
# #' calculation should only be limited within this range. Default \code{NULL}
# #' uses all datasets.
# #' @param features For container object methods. Valid row subsetting index that
# #' selects features. Default \code{NULL}.
# #' @export
# #' @method plotQuaternary liger
# plotQuaternary.liger <- function(
#         object,
#         clusterVar,
#         features = NULL,
#         useDatasets = NULL,
#         ...
# ) {
#     values <- .ligerPrepare(object, clusterVar, features, useDatasets)
#     plotQuaternary(values[[1]], clusterVar = values[[2]], ...)
# }

#' Create GIF image for dynamic rotating view of 3D quaternary simplex plot
#' @param object Input object that \code{\link{plotQuaternary}} accepts.
#' @param ... All other arguments needed for \code{\link{plotQuaternary}},
#' except \code{theta} which controls the rotation.
#' @param useCluster One cluster that exists in \code{clusterVar}, if users
#' need to view the plot for specific group. Default \code{NULL} plot all cells.
#' @param gifPath Output GIF image file path. Default \code{"quaternary.gif"}
#' @param tmpDir A temprorary directory to store all PNG files for all
#' perspectives created. Default \code{file.path(getwd(),
#' "quaternary_gif_tmp/")}.
#' @param fps Number of frame per second. Default \code{10}.
#' @param degreePerFrame Number of degree that the tetrahedron is rotated per
#' frame. Default \code{10}.
#' @return No object is returned. The \code{tmpDir} folder will be created with
#' \code{360 / degreePerFrame} PNG image files in it. A GIF image file will be
#' created at \code{gifPath}.
#' @export
#' @examples
#' rnaLog <- colNormalize(rnaRaw, 1e4, TRUE)
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' writeQuaternaryGIF(rnaLog[gene, ], rnaCluster, c("RE", "OS", "CH", "ORT"))
writeQuaternaryGIF <- function(
        object,
        ...,
        useCluster = NULL,
        gifPath = "quaternary.gif",
        tmpDir = file.path(getwd(), "quaternary_gif_tmp/"),
        fps = 10,
        degreePerFrame = 10
) {
    if (!requireNamespace("magick", quietly = TRUE)) {
        stop("Package 'magick' must be installed for creating GIF figure.\n",
             "See https://cran.r-project.org/web/packages/magick/vignettes/",
             "intro.html#Installing_magick for detail.")
    }
    if (100 %% fps != 0) {
        stop("FPS must be a factor of 100.")
    }
    allImgPath <- c()
    if (!dir.exists(tmpDir)) dir.create(tmpDir)

    message("Generating quanternary simplex plots...")
    allTheta <- seq(0, 360, degreePerFrame)
    pb <- utils::txtProgressBar(1, length(allTheta), style = 3)
    for (i in seq_along(allTheta)) {
        theta <- allTheta[i]
        filename <- paste0("quaternary_theta", theta, ".png")
        filename <- file.path(tmpDir, filename)
        allImgPath <- c(allImgPath, filename)
        if (is.null(useCluster)) {
            p <- plotQuaternary(object, ..., theta = theta)
        } else {
            pl <- plotQuaternary(object, ..., splitCluster = TRUE,
                                 theta = theta)
            p <- pl[[useCluster]]
        }

        grDevices::png(filename)
        print(p)
        grDevices::dev.off()
        utils::setTxtProgressBar(pb, value = i)
    }
    cat("\n")
    imgList <- lapply(allImgPath, magick::image_read)
    imgJoined <- magick::image_join(imgList)
    imgAnimated <- magick::image_animate(imgJoined, fps = fps)
    message("Generate the GIF: ", gifPath)
    magick::image_write(image = imgAnimated,
                        path = gifPath)
}
