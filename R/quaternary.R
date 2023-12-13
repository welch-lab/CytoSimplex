#' Create quaternary simplex plots
#' @description
#' Create quaternary plots that show similarity between single cells and
#' selected four terminals in a baricentric coordinate.
#'
#' See \code{\link{plotTernary}} for more details on methodologies.
#'
#' A dynamic rotating view in a GIF image file can be created with
#' \code{\link{writeQuaternaryGIF}}. Package \code{magick} must be installed in
#' advance. Linux users may refer to this
#' \href{https://cran.r-project.org/package=magick/vignettes/intro.html#Build_from_source}{installation guide}.
#' @param x Input data. Can be a \code{matrix} or \code{dgCMatrix} object with
#' cells as columns, a \code{Seurat} or \code{SingleCellExperiment} object.
#' "simMat" method takes intermediate values.
#' @param ... Arguments passed to other methods.
#' @rdname plotQuaternary
#' @export plotQuaternary
#' @return For "simMat" method, a "plist" (plot3D package product) object. For
#' other methods, a "plist" object when \code{splitCluster = FALSE}, or a list
#' of "plist" objects when \code{splitCluster = TRUE}. A "plist" object can be
#' viewed with \code{print()}, \code{show()} or a direct run of the object
#' variable name in interactive console.
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' plotQuaternary(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"), gene)
plotQuaternary <- function(x, ...) {
    UseMethod('plotQuaternary', x)
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
#' spreading on figure. Default \code{0.05}.
#' @param scale Whether to min-max scale the distance matrix by clusters.
#' Default \code{TRUE}.
#' @param returnData Logical. Whether to return similarity and aggregated
#' velocity data if applicable instead of generating plot. Default \code{FALSE}.
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary default
plotQuaternary.default <- function(
        x,
        clusterVar,
        vertices,
        features = NULL,
        veloGraph = NULL,
        byCluster = NULL,
        processed = FALSE,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        force = FALSE,
        sigma = 0.05,
        scale = TRUE,
        dotColor = "grey60",
        returnData = FALSE,
        ...
) {
    method <- match.arg(method)
    vcheck <- .checkVertex(x, clusterVar, vertices, n = 4)
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

    if (isTRUE(returnData)) return(list(sim = simMat, velo = veloMat,
                                        originalCluster = clusterVar,
                                        mappedCluster = vertClust))

    if (is.null(byCluster)) {
        return(plotQuaternary(x = simMat, veloMat = veloMat,
               dotColor = dotColor, ...))
    } else {
        if (identical(byCluster, "all")) {
            pl <- lapply(levels(clusterVar), function(clust) {
                plotQuaternary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               title = clust,
                               dotColor = dotColor[clusterVar == clust], ...)
            })
            pl$allCells <- plotQuaternary(simMat,
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
                plotQuaternary(simMat[clusterVar == clust,],
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
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary Seurat
#' @examples
#' \donttest{
#' # Seurat example
#' library(Seurat)
#' srt <- CreateSeuratObject(rnaRaw)
#' Idents(srt) <- rnaCluster
#' gene <- selectTopFeatures(srt, vertices = c("OS", "RE", "CH", "ORT"))
#' plotQuaternary(srt, features = gene,
#'                vertices = c("OS", "RE", "CH", "ORT"))
#' }
plotQuaternary.Seurat <- function(
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
    plotQuaternary(values[[1]], clusterVar = values[[2]], processed = processed,
                   ...)
}

#' @param assay.type For "SingleCellExperiment" methods. Which assay to use for
#' calculating the similarity. Default \code{"counts"}.
#' @rdname plotQuaternary
#' @export
#' @method plotQuaternary SingleCellExperiment
#' @examples
#' \donttest{
#' # SingleCellExperiment example
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#' colLabels(sce) <- rnaCluster
#' gene <- selectTopFeatures(sce, vertices = c("OS", "RE", "CH", "ORT"))
#' plotQuaternary(sce, features = gene,
#'                vertices = c("OS", "RE", "CH", "ORT"))
#' }
plotQuaternary.SingleCellExperiment <- function(
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
    plotQuaternary(values[[1]], clusterVar = values[[2]], processed = processed,
                   ...)
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
#         x,
#         clusterVar,
#         features = NULL,
#         useDatasets = NULL,
#         ...
# ) {
#     values <- .ligerPrepare(x, clusterVar, features, useDatasets)
#     plotQuaternary(values[[1]], clusterVar = values[[2]], ...)
# }

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
#' @param arrowAngle,arrowLen Arrow aesthetics passed to
#' \code{grid::\link[grid]{arrow}}. The length of the arrow will be internally
#' converted to unit onject in inches. Default \code{20} and \code{0.1}.
#' @param edgeLinewidth Controls the linewidth of the edges of the tetrahedron.
#' Default \code{1}.
#' @param vertexLabelSize Numeric, size of vertex text label relative to default
#' size. Default \code{1}.
#' @param title Title text of the plot. Default \code{NULL}.
#' @param titleSize,titleColor Setting on the main title text. Default \code{1},
#' and \code{"black"}.
#' @param theta,phi Numeric scalar. The angles defining the viewing direction.
#' \code{theta} gives the azimuthal direction and \code{phi} the colatitude.
#' Default \code{20} and \code{0}.
#' @param interactive Logical. Whether to use "rgl" library to create
#' interactive device. Default \code{FALSE}.
#' @export
#' @method plotQuaternary simMat
plotQuaternary.simMat <- function(
        x,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.2,
        dotSize = 0.6,
        dotColor = "grey60",
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
        arrowLinewidth = 0.6,
        arrowAngle = 20,
        arrowLen = 0.1,
        vertexLabelSize = 1,
        edgeLinewidth = 1,
        title = NULL,
        titleSize = 1,
        titleColor = "black",
        theta = 20,
        phi = 0,
        interactive = FALSE,
        ...
) {
    if (isTRUE(interactive)) {
        .plotQuatRGL(
            simMat = x, veloMat = veloMat, nGrid = nGrid, radius = radius,
            dotSize = dotSize, dotColor = dotColor, title = title,
            labelColors = labelColors, arrowLinewidth = arrowLinewidth,
            vertexLabelSize = vertexLabelSize
        )
    } else {
        .plotQuat(
            simMat = x, veloMat = veloMat, nGrid = nGrid, radius = radius,
            dotSize = dotSize, dotColor = dotColor, title = title,
            titleSize = titleSize, titleColor = titleColor,
            labelColors = labelColors, arrowLinewidth = arrowLinewidth,
            arrowAngle = arrowAngle, arrowLen = arrowLen,
            theta = theta, phi = phi, vertexLabelSize = vertexLabelSize,
            edgeLinewidth = edgeLinewidth
        )
    }
}

.plotQuat <- function(
        simMat,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.2,
        dotSize = 0.6,
        dotColor = "grey60",
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
        arrowLinewidth = 0.6,
        arrowAngle = 20,
        arrowLen = 0.1,
        vertexLabelSize = 1,
        edgeLinewidth = 1,
        title = NULL,
        titleSize = 1,
        titleColor = "black",
        theta = 20,
        phi = 0
) {
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    cellCart <- as.matrix(simMat) %*% tetra

    tetraVertices <- rotateByZAxis(tetra, theta)
    cellCart <- rotateByZAxis(cellCart, theta)
    # Plot data
    grDevices::pdf(nullfile())
    scatter3D(cellCart[,1], cellCart[,2], cellCart[,3],
              main = list(title, cex = titleSize, col = titleColor), outer = FALSE,
              xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), zlim = c(0, 1.7),
              alpha = 0.8, col = dotColor, cex = dotSize/2, pch = 16, d = 3,
              colkey = list(plot = FALSE), expand = 0.7,
              box = FALSE, theta = 0, phi = phi, plot = FALSE)
    lines3D(tetraVertices[c(1,2,3,4,1,3,1,2,4), 1],
            tetraVertices[c(1,2,3,4,1,3,1,2,4), 2],
            tetraVertices[c(1,2,3,4,1,3,1,2,4), 3],
            col = "grey", add = TRUE, plot = FALSE, lwd = edgeLinewidth)
    text3D(tetraVertices[,1]*1.1, tetraVertices[,2]*1.1, tetraVertices[,3]*1.03,
           colnames(simMat)[seq(4)], col = labelColors, adj = 0.5, font = 2,
           cex = vertexLabelSize, add = TRUE, plot = FALSE)

    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = simMat, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            subcoord <- arrowCoords[[i]]
            subcoord[,seq(3)] <- rotateByZAxis(subcoord[,seq(3)], theta)
            subcoord[,seq(4,6)] <- rotateByZAxis(subcoord[,seq(4,6)], theta)
            arrows3D(subcoord[,1], subcoord[,2], subcoord[,3],
                     subcoord[,4], subcoord[,5], subcoord[,6],
                     lwd = arrowLinewidth, length = arrowLen,
                     angle = arrowAngle, col = labelColors[i], add = TRUE,
                     plot = FALSE)
        }
    }
    grDevices::dev.off()
    getplist()
}

setOldClass("plist")

#' @rdname show-plist
#' @title Show plist object produced with plot3D package
#' @param object,x plist object
#' @param ... Graphic parameters passed to \code{\link{plot}}. \code{mar} is
#' pre-specified.
#' @export
#' @return No return value. It displays the plot described in a 'plist' object
#' returned by \code{\link{plotQuaternary}}, internally created by package
#' 'plot3D'.
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' plistObj <- plotQuaternary(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"), gene)
#' print(plistObj)
#' # equivalent to
#' show(plistObj)
setMethod("show", "plist", function(object) {
    print(object)
}
)

#' @rdname show-plist
#' @method print plist
#' @export
print.plist <- function(x, ...) {
    oldpar <- graphics::par(no.readonly = TRUE) # code line i
    on.exit(graphics::par(oldpar))
    graphics::par(mar = c(0, 0, 0, 0))
    plot(x, ...)
}

# Rotate cartesien coordinate around Z-axis
# coord - N x 3 matrix
# theta - degree to rotate in azimuthal direction
# Return - rotated N x 3 matrix
rotateByZAxis <- function(coord, theta) {
    if (theta %% 360 == 0) return(coord)
    # Convert theta to radians
    theta_rad <- theta * pi / 180

    # Create rotation matrix
    rotation_matrix <- matrix(c(cos(theta_rad), -sin(theta_rad), 0,
                                sin(theta_rad), cos(theta_rad), 0,
                                0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)

    # Perform rotation
    rotated_points <- as.matrix(coord) %*% t(rotation_matrix)

    return(rotated_points)
}

.plotQuatRGL <- function(
    simMat,
    veloMat = NULL,
    nGrid = 10,
    radius = 0.2,
    dotSize = 0.6,
    dotColor = "grey60",
    labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
    arrowLinewidth = 0.1,
    vertexLabelSize = 1,
    title = NULL
) {
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("Package \"rgl\" is required for generating the interactive ",
             "view of a quaternary simplex plot. Please install with running:",
             "\ninstall.packages(\"rgl\")")
    }
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    cellCart <- as.matrix(simMat) %*% tetra

    # Plot the tetrahedron edges with RGL
    rgl::open3d()
    rgl::bg3d("white")
    rgl::points3d(cellCart, color = dotColor, size = dotSize)
    rgl::lines3d(tetra[c(1,2,3,4,1,3,1,2,4),1],
                 tetra[c(1,2,3,4,1,3,1,2,4),2],
                 tetra[c(1,2,3,4,1,3,1,2,4),3])
    rgl::text3d(tetra[,1]*1.1, tetra[,2]*1.1, tetra[,3]*1.1, colnames(simMat),
                color = labelColors, cex = vertexLabelSize)
    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = simMat, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)

        for (i in seq_along(arrowCoords)) {
            subcoord <- arrowCoords[[i]]
            for (j in seq(nrow(subcoord))) {
                rgl::arrow3d(c(subcoord[j,1], subcoord[j,2], subcoord[j,3]),
                             c(subcoord[j,4], subcoord[j,5], subcoord[j,6]),
                             barblen = 0.005, width = 0.01,
                             theta = 70, lwd = arrowLinewidth,
                             type = "line", color = labelColors[i],
                             thickness = 1)
            }
        }
    }
    if (!is.null(title)) rgl::title3d(title)
}


#' Create GIF image for dynamic rotating view of 3D quaternary simplex plot
#' @param x Input object that \code{\link{plotQuaternary}} accepts.
#' @param ... All other arguments needed for \code{\link{plotQuaternary}}. Must
#' be specified with exact argument names instead of a positional manner.
#' @param cluster One cluster that exists in \code{clusterVar}, if users
#' need to view the plot for specific group. Default \code{NULL} plot all cells.
#' @param gifPath Output GIF image file path. Default \code{"quaternary.gif"}
#' @param tmpDir A temprorary directory to store all PNG files for all
#' perspectives created. Default \code{tempdir()}.
#' @param fps Number of frame per second, must be a factor of 100. Default
#' \code{10}.
#' @param degreePerFrame Number of degree that the tetrahedron is rotated per
#' frame. Default \code{10}.
#' @return No object is returned. The \code{tmpDir} folder will be created with
#' \code{360 / degreePerFrame} PNG image files in it. A GIF image file will be
#' created at \code{gifPath}.
#' @export
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' \donttest{
#' writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster, features = gene,
#'                    vertices = c("RE", "OS", "CH", "ORT"),
#'                    gifPath = paste0(tempfile(), ".gif"))
#' }
writeQuaternaryGIF <- function(
        x,
        ...,
        cluster = NULL,
        gifPath = "quaternary.gif",
        tmpDir = tempdir(),
        fps = 10,
        degreePerFrame = 10
) {
    if (!requireNamespace("magick", quietly = TRUE)) {
        stop("Package 'magick' must be installed for creating GIF figure.\n", # nocov
             "See https://cran.r-project.org/web/packages/magick/vignettes/", # nocov
             "intro.html#Installing_magick for detail.") # nocov
    }
    if (100 %% fps != 0) {
        stop("FPS must be a factor of 100.")
    }
    allImgPath <- c()
    if (!dir.exists(tmpDir)) dir.create(tmpDir)

    methodArgs <- list(...)
    if (is.null(names(methodArgs)) || "" %in% names(methodArgs)) {
        stop("Please set `...` arguments with explicit argument names")
    }
    ignore <- c("theta", "byCluster", "returnData", "interactive", "veloMat")
    if (any(names(methodArgs) %in% ignore)) {
        warning("Arguments ignored: ",
                paste(names(methodArgs)[names(methodArgs) %in% ignore],
                      collapse = ", "))
    }
    methodArgs[ignore] <- NULL

    plotData <- do.call(plotQuaternary, c(list(x = x, returnData = TRUE),
                                          methodArgs))
    clusterVar <- plotData$originalCluster
    if (!is.null(cluster)) {
        if (length(cluster) > 1)
            stop("Can only generate GIF for one cluster at a time")
        if (!cluster %in% clusterVar)
            stop("\"", cluster, "\" is not an available cluster.")
        simMat <- plotData$sim[clusterVar == cluster,]
        veloMat <- plotData$velo[clusterVar == cluster,]
        if (!"title" %in% names(methodArgs)) methodArgs$title <- cluster
    } else {
        simMat <- plotData$sim
        veloMat <- plotData$velo
    }

    message("Generating quanternary simplex plots...")
    allTheta <- seq(0, 360, degreePerFrame)
    pb <- utils::txtProgressBar(1, length(allTheta), style = 3)
    for (i in seq_along(allTheta)) {
        theta <- allTheta[i]
        filename <- paste0("quaternary_theta", theta, ".png")
        filename <- file.path(tmpDir, filename)
        allImgPath <- c(allImgPath, filename)
        p <- do.call(plotQuaternary,
                     c(list(x = simMat, veloMat = veloMat, theta = theta,
                            interactive = FALSE),
                       methodArgs))

        grDevices::png(filename)
        print(p)
        grDevices::dev.off()
        utils::setTxtProgressBar(pb, value = i)
    }
    cat("\n")
    imgList <- lapply(allImgPath, magick::image_read)
    imgJoined <- magick::image_join(imgList)
    imgAnimated <- magick::image_animate(imgJoined, fps = fps)
    message("Generated the GIF: ", gifPath)
    magick::image_write(image = imgAnimated,
                        path = gifPath)
    invisible()
}
