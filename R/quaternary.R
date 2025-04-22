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
#' @param ... Arguments passed on to other S3 methods.
#' @return By default, a "plotly" object. When \code{interactive = FALSE}, a
#' "quatPlot" object when \code{byCluster} is not specified, or a "list" of
#' "quatPlot" objects when \code{byCluster} is specified. When
#' \code{returnData = TRUE}, a list of similarity matrix and aggregated velocity
#' matrix is returned.
#' @family plotQuaternary
#' @rdname plotQuaternary
#' @export
plotQuaternary <- function(x, ...) {
    UseMethod('plotQuaternary', x)
}

#' @inherit plotQuaternary title description
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
#' @inheritDotParams plotQuaternary.simMat nGrid radius dotSize labelColors arrowLinewidth arrowAngle arrowLen vertexLabelSize edgeLinewidth title titleSize titleColor theta phi interactive
#' @rdname plotQuaternary
#' @method plotQuaternary default
#' @export
#' @examples
#' gene <- selectTopFeatures(
#'     x = rnaRaw,
#'     clusterVar = rnaCluster,
#'     vertices = c("RE", "OS", "CH", "ORT")
#' )
#' plotQuaternary(
#'     x = rnaRaw,
#'     clusterVar = rnaCluster,
#'     vertices = c("RE", "OS", "CH", "ORT"),
#'     features = gene
#' )
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
    vcheck <- .checkVertex(x, clusterVar, vertices, n = 4)
    vertClust <- vcheck[[1]]
    vertices <- vcheck[[2]]

    colorArg <- .checkDotColor(
        n = ncol(x), cluster = clusterVar, colorBy = dotColorBy,
        colors = dotColor, palette = palette, direction = direction,
        breaks = breaks,
        legendTitle = legendTitle
    )
    if (isFALSE(processed) && !is.rawCounts(x)) {
        cli::cli_alert_warning(
            "Input matrix is not raw counts (integers). Results may be affected."
        )
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
            cli::cli_abort(
                c("x" = "{.var veloGraph} must have dimension of {nrow(simMat)} x {nrow(simMat)} and has dimnames covering all cells in {.var x}.",
                  "i" = "Given dim: {nrow(veloGraph)} x {ncol(veloGraph)}")
            )
        }
        veloGraph <- veloGraph[rownames(simMat), rownames(simMat)]
        veloMat <- aggrVeloGraph(veloGraph, clusterVar = vertClust,
                                 vertices = vertices)
    }

    if (isTRUE(returnData)) {
        return(
            list(
                sim = simMat,
                velo = veloMat,
                originalCluster = clusterVar,
                mappedCluster = vertClust,
                colorArg = colorArg
            )
        )
    }

    if (is.null(byCluster)) {
        return(plotQuaternary(x = simMat, veloMat = veloMat,
               colorArg = colorArg, ...))
    } else {
        # If specified in `...`, use that. Otherwise, always assume the default
        # to be `TRUE`. Note an instance of this parameter does not exist yet,
        # thus created.
        interactive <- list(...)$interactive %||% TRUE
        if (isTRUE(interactive)) {
            cli::cli_alert_warning(
                "Interactive view with {.pkg plotly} can show desired cluster(s) by clicking on the legend. Ignoring {.var byCluster}."
            )
            return(plotQuaternary(
                x = simMat, veloMat = veloMat, colorArg = colorArg, ...
            ))
        }
        if (identical(byCluster, "all")) {
            pl <- lapply(levels(clusterVar), function(clust) {
                plotQuaternary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               title = clust,
                               colorArg = colorArg[clusterVar == clust], ...)
            })
            pl$allCells <- plotQuaternary(simMat,
                                          veloMat = veloMat,
                                          colorArg = colorArg,
                                          title = "All cells", ...)
            names(pl) <- c(levels(clusterVar), "allCells")
        } else {
            if (any(!byCluster %in% levels(clusterVar))) {
                cli::cli_abort(
                    "{.var byCluster} must be either a vector of cluster names or just {.val all}."
                )
            }
            pl <- lapply(byCluster, function(clust) {
                plotQuaternary(simMat[clusterVar == clust,],
                               veloMat = veloMat[clusterVar == clust,],
                               title = clust,
                               colorArg = colorArg[clusterVar == clust], ...)
            })
            names(pl) <- byCluster
        }
        return(pl)
    }
}

#' Craete quaternary simplex plot with Seurat object
#' @inherit plotQuaternary description return
#' @param x A Seurat object
#' @param layer Layer in the specified assay to use. Default \code{"counts"}.
#' @param assay The assay to get data from. Default \code{NULL} uses the default
#' assay.
#' @param clusterVar A variable name in meta.data (\code{x[[]]}). Default
#' \code{NULL} uses \code{Idents(x)}.
#' @inheritParams plotQuaternary.default
#' @inheritDotParams plotQuaternary.default vertices features veloGraph byCluster method force sigma scale dotColor palette direction breaks legendTitle returnData
#' @inheritDotParams plotQuaternary.simMat nGrid radius dotSize labelColors arrowLinewidth arrowAngle arrowLen vertexLabelSize edgeLinewidth title titleSize titleColor theta phi interactive
#' @rdname plotQuaternary.Seurat
#' @family plotQuaternary
#' @export
#' @method plotQuaternary Seurat
#' @examples
#' \donttest{
#' # Seurat example
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   require(Seurat)
#'   srt <- CreateSeuratObject(rnaRaw)
#'   Idents(srt) <- rnaCluster
#'   gene <- selectTopFeatures(srt, vertices = c("OS", "RE", "CH", "ORT"))
#'   plotQuaternary(srt, features = gene,
#'                  vertices = c("OS", "RE", "CH", "ORT"))
#' }
#' }
plotQuaternary.Seurat <- function(
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
    plotQuaternary(values[[1]], clusterVar = values[[2]], processed = processed,
                   dotColorBy = dotColorBy, legendTitle = legendTitle,
                   ...)
}

#' Create quaternary simplex plot with SingleCellExperiment object
#' @family plotQuaternary
#' @inherit plotQuaternary description return
#' @param x A SingleCellExperiment object
#' @param assay.type Assay to use for calculating the similarity. Default
#' \code{"counts"}.
#' @inheritParams plotQuaternary.default
#' @param clusterVar A variable name in \code{colData(x)}. Default \code{NULL}
#' uses \code{colLabels(x)}.
#' @inheritDotParams plotQuaternary.default vertices features veloGraph byCluster method force sigma scale dotColor palette direction breaks legendTitle returnData
#' @inheritDotParams plotQuaternary.simMat nGrid radius dotSize labelColors arrowLinewidth arrowAngle arrowLen vertexLabelSize edgeLinewidth title titleSize titleColor theta phi interactive
#' @rdname plotQuaternary.SingleCellExperiment
#' @export
#' @method plotQuaternary SingleCellExperiment
#' @examples
#' \donttest{
#' # SingleCellExperiment example
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   require(SingleCellExperiment)
#'   sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#'   colLabels(sce) <- rnaCluster
#'   gene <- selectTopFeatures(sce, vertices = c("OS", "RE", "CH", "ORT"))
#'   plotQuaternary(sce, features = gene,
#'                  vertices = c("OS", "RE", "CH", "ORT"))
#' }
#' }
plotQuaternary.SingleCellExperiment <- function(
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
    plotQuaternary(values[[1]], clusterVar = values[[2]], processed = processed,
                   dotColorBy = dotColorBy, legendTitle = legendTitle,
                   ...)
}

#' Create quaternary simplex plot with pre-calculated similarity matrix
#' @rdname plotQuaternary.simMat
#' @param x simMat object, n cells by 4 vertices, each row summing to 1.
#' @param veloMat Aggregated velocity matrix. Output of \code{aggrVeloGraph}.
#' @param nGrid Number of grids along the x-axis of the tetrahedron. Default
#' \code{10}.
#' @param radius Arrow length of unit velocity. Lower this when arrows point
#' outside of the tetrahedron. Default \code{0.2}.
#' @param dotSize Size of each dot. Default \code{0.6} for static figure, and
#' \code{4} for interactive view.
#' @param colorArg A "colorArg" object, internally prepared by
#' \code{\link{plotQuaternary.default}}. Default \code{NULL}.
#' @param labelColors Colors of the vertex labels. Default
#' \code{c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF")} (blue, red,
#' green and purple).
#' @param arrowLinewidth Arrow aesthetics. Default \code{1.6} for interactive
#' view, \code{0.6} for static figure.
#' @param arrowAngle,arrowLen Arrow aesthetics passed to TODOOOO
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
#' @param interactive Logical. Whether to display plotly interactive view.
#' Default \code{TRUE}.
#' @param ... Not used
#' @method plotQuaternary simMat
#' @return A "quatPlot" object, can be displayed by printing.
#' @export
plotQuaternary.simMat <- function(
        x,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.2,
        dotSize = NULL,
        colorArg = NULL,
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
        arrowLinewidth = NULL,
        arrowAngle = 20,
        arrowLen = 0.1,
        vertexLabelSize = NULL,
        edgeLinewidth = 1,
        title = NULL,
        titleSize = 1,
        titleColor = "black",
        theta = 20,
        phi = 0,
        interactive = TRUE,
        ...
) {
    if (isTRUE(interactive)) {
        .plotQuatPlotly(
            simMat = x, veloMat = veloMat, nGrid = nGrid, radius = radius,
            dotSize = dotSize, colorArg = colorArg,
            title = title, labelColors = labelColors,
            arrowLinewidth = arrowLinewidth, vertexLabelSize = vertexLabelSize,
            arrowLen = arrowLen, arrowAngle = arrowAngle
        )
    } else {
        .plotQuat(
            simMat = x, veloMat = veloMat, nGrid = nGrid, radius = radius,
            dotSize = dotSize, colorArg = colorArg,
            title = title, titleSize = titleSize, titleColor = titleColor,
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
        colorArg = NULL,
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
    dotSize <- dotSize %||% 0.6
    vertexLabelSize <- vertexLabelSize %||% 1
    arrowLinewidth <- arrowLinewidth %||% 0.6
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    cellCart <- as.matrix(simMat) %*% tetra

    tetraVertices <- rotateByZAxis(tetra, theta)
    cellCart <- rotateByZAxis(cellCart, theta)

    # Workaround to match up with dumb plot3D coloring schema especially when
    # it comes to categorical
    # - `colvar` has to be numeric leveling: 1, 1, 2, 3
    # - `col` has to be as many color values as the number of levels of a
    #   non-factor object passed to `colvar`
    colvar <- NULL
    col <- NULL
    breaks <- NULL
    legend <- NULL
    legendTitle <- colorArg$palette$legendTitle
    # PH - placeholder for colkey
    colkeyPH <- FALSE
    if (colorArg$type == "continuous") {
        colvar <- colorArg$colorBy
        col <- colorArg$palette$colors
        breaks <- colorArg$palette$breaks
        colkeyPH <- NULL
    } else if (colorArg$type == "categorical") {
        colvar <- as.numeric(droplevels(colorArg$colorBy))
        allCol <- character(nlevels(colorArg$colorBy))
        for (i in seq(nlevels(colorArg$colorBy))) {
            clustCol <- unique(colorArg$colors[colorArg$colorBy == levels(colorArg$colorBy)[i]])
            if (length(clustCol) == 0) clustCol <- NA
            allCol[i] <- clustCol
        }
        col <- allCol[!is.na(allCol)]
        colkeyPH <- list(plot = FALSE)
        # Make arguments for using base legend instead of using plot3D colkey
        legend <- list(
            x = "right",
            legend = levels(colorArg$colorBy),
            col = allCol,
            bty = "n",
            title = legendTitle,
            pch = 16,
            pt.cex = 1.5,
            y.intersp = 1.5,
            title.cex = 1.2,
            title.font = 2 # for bold face
        )
    } else {
        colvar <- NULL
        col <- colorArg$colors
    }
    # Plot data
    grDevices::pdf(nullfile())
    scatter3D(x = cellCart[,1], y = cellCart[,2], z = cellCart[,3],
              main = list(title, cex = titleSize, col = titleColor),
              outer = FALSE,
              xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), zlim = c(0, 1.7),
              alpha = 0.8, colvar = colvar, col = col,
              colkey = colkeyPH,
              clab = if (colorArg$type == "continuous") c(legendTitle, "") else NULL,
              cex = dotSize/2, pch = 16, d = 3, expand = 0.7,
              box = FALSE, theta = 0, phi = phi, plot = FALSE)

    lines3D(tetraVertices[c(1,2,3,4,1,3,2,4), 1],
            tetraVertices[c(1,2,3,4,1,3,2,4), 2],
            tetraVertices[c(1,2,3,4,1,3,2,4), 3],
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
    plist <- getplist()
    class(plist) <- c("quatPlot", class(plist))
    if (!is.null(legend)) plist$legend <- legend
    return(plist)
}

#' @rdname print-quatPlot
#' @title Show plist object produced with plot3D package
#' @param x quatPlot object, returned by \code{\link{plotQuaternary}} when
#' \code{interactive = FALSE}.
#' @param ... Graphic parameters passed to \code{\link{plot}}. \code{mar} is
#' pre-specified.
#' @export
#' @return No return value. It displays the plot described in a 'plist' object
#' returned by \code{\link{plotQuaternary}}, internally created by package
#' 'plot3D'.
#' @method print quatPlot
#' @examples
#' gene <- selectTopFeatures(
#'     x = rnaRaw,
#'     clusterVar = rnaCluster,
#'     vertices = c("RE", "OS", "CH", "ORT")
#' )
#' quat <- plotQuaternary(
#'     x = rnaRaw,
#'     clusterVar = rnaCluster,
#'     vertices = c("RE", "OS", "CH", "ORT"),
#'     features = gene,
#'     interactive = FALSE
#' )
#' quat; print(quat)
print.quatPlot <- function(x, ...) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mar = c(0, 0, 0, 0))
    legend.par <- x$legend
    plot(x, ...)
    if (!is.null(legend.par)) {
        do.call(graphics::legend, legend.par)
    }
    invisible(NULL)
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

# .plotQuatRGL <- function(
#     simMat,
#     veloMat = NULL,
#     nGrid = 10,
#     radius = 0.2,
#     dotSize = 0.6,
#     dotColor = "grey60",
#     labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
#     arrowLinewidth = 0.1,
#     vertexLabelSize = 1,
#     title = NULL
# ) {
#     if (!requireNamespace("rgl", quietly = TRUE)) {
#         stop("Package \"rgl\" is required for generating the interactive ",
#              "view of a quaternary simplex plot. Please install with running:",
#              "\ninstall.packages(\"rgl\")")
#     }
#     # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
#     cellCart <- as.matrix(simMat) %*% tetra
#
#     # Plot the tetrahedron edges with RGL
#     rgl::open3d()
#     rgl::bg3d("white")
#     rgl::points3d(cellCart, color = dotColor, size = dotSize)
#     rgl::lines3d(tetra[c(1,2,3,4,1,3,1,2,4),1],
#                  tetra[c(1,2,3,4,1,3,1,2,4),2],
#                  tetra[c(1,2,3,4,1,3,1,2,4),3])
#     rgl::text3d(tetra[,1]*1.1, tetra[,2]*1.1, tetra[,3]*1.1, colnames(simMat),
#                 color = labelColors, cex = vertexLabelSize)
#     if (!is.null(veloMat)) {
#         arrowCoords <- calcGridVelo(simMat = simMat, veloMat = veloMat,
#                                     nGrid = nGrid, radius = radius)
#
#         for (i in seq_along(arrowCoords)) {
#             subcoord <- arrowCoords[[i]]
#             for (j in seq(nrow(subcoord))) {
#                 rgl::arrow3d(c(subcoord[j,1], subcoord[j,2], subcoord[j,3]),
#                              c(subcoord[j,4], subcoord[j,5], subcoord[j,6]),
#                              barblen = 0.005, width = 0.01,
#                              theta = 70, lwd = arrowLinewidth,
#                              type = "line", color = labelColors[i],
#                              thickness = 1)
#             }
#         }
#     }
#     if (!is.null(title)) rgl::title3d(title)
# }

.plotQuatPlotly <- function(
        simMat,
        veloMat = NULL,
        nGrid = 10,
        radius = 0.2,
        dotSize = 1,
        colorArg = NULL,
        labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF"),
        arrowLinewidth = 1.2,
        arrowAngle = 20,
        arrowLen = 0.1,
        vertexLabelSize = 16,
        title = NULL
) {
    vertexLabelSize <- vertexLabelSize %||% 16
    dotSize <- dotSize %||% 4
    arrowLinewidth <- arrowLinewidth %||% 1.6
    # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
    cellCart <- as.matrix(simMat) %*% tetra
    # Axis settings to disable everything
    ax <- list(
        gridcolor = 'rgba(255, 255, 255, 255)',
        zerolinecolor = 'rgba(255, 255, 255, 255)',
        showticklabels = FALSE,
        title = '',
        showspikes = FALSE
    )
    # Dot plot
    # it still needs a data.frame as input
    dotsCart <- as.data.frame(cellCart)
    colnames(dotsCart) <- c("x", "y", "z")
    dotsCart$ID <- rownames(cellCart)
    dotsCart[[colorArg$palette$legendTitle]] <- colorArg$colorBy

    hoverDF <- as.data.frame(simMat)
    hoverDF <- cbind(ID = rownames(simMat), hoverDF)
    hoverDF[[colorArg$palette$legendTitle]] <- colorArg$colorBy
    # dotColorBy <- colorArg$colorBy
    # dotColor <- colorArg$colors
    dotColor <- colorArg$palette$colors
    lgdTtlFormula <- stats::formula(paste0("~", colorArg$palette$legendTitle))
    # Draw the dots on one plot
    figDots <- plotly::plot_ly(
        dotsCart, x = ~x, y = ~y, z = ~z, text = ~ID,
        size = I(dotSize), color = lgdTtlFormula, colors = dotColor,
        type = "scatter3d", mode = "markers",
        customdata = .plotly_formatCustomData(hoverDF),
        hoverinfo = "text",
        hovertemplate = "%{customdata}"
    )

    # Draw the tetrahedron edges on another plot
    figEdges <- plotly::plot_ly(
        as.data.frame(tetra)[c(1,2,3,4,1,3,4,2),],
        name = "Tetrahedron", x = ~x, y = ~y, z = ~z,
        color = I("black"),
        type = "scatter3d", mode = "lines",
        hoverinfo = "skip", showlegend = FALSE
    )

    subplotList <- list(cells = figDots, tetra = figEdges)

    # Make vertext label annotation
    # Move position ou
    tetraCenter <- colMeans(tetra)
    # Vector that points from tetraCenter to each vertex
    outerVec <- tetra - tetraCenter
    # Normalize to 0.1 length
    outerVec <- outerVec / sqrt(rowSums(outerVec^2)) * 0.2
    # Vertex label position
    vertPos <- tetra + outerVec
    vertDF <- as.data.frame(vertPos)
    vertDF$Vertex <- .htmlbold(colnames(simMat))
    for (i in seq_len(nrow(vertDF))) {
        # Don't know why adding all text label at a time, with setting
        # different colors, will always issue warning messages like
        # "error_x.color ..."
        subplotList[[paste0("vertLabel_", i)]] <- plotly::plot_ly(
            vertDF[i, , drop = FALSE],
            x = ~x, y = ~y, z = ~z, text = ~Vertex, color = I(labelColors[i]),
            type = "scatter3d", mode = "text", hoverinfo = "skip",
            textposition = "middle center", showlegend = FALSE
        )
    }

    if (!is.null(veloMat)) {
        arrowCoords <- calcGridVelo(simMat = simMat, veloMat = veloMat,
                                    nGrid = nGrid, radius = radius)
        for (i in seq_along(arrowCoords)) {
            legendGroup <- paste0("velo_", i)
            subcoord <- arrowCoords[[i]]
            # For plotly, we don't have another way to create a "segment" of
            # arrow type, but we have to draw a line and a cone.
            for (l in seq_len(nrow(subcoord))) {
                # Now it's for each arrow
                line <- data.frame(
                    x = c(subcoord[l, 1], subcoord[l, 4]),
                    y = c(subcoord[l, 2], subcoord[l, 5]),
                    z = c(subcoord[l, 3], subcoord[l, 6])
                )
                length = sqrt(sum((subcoord[l, 4:6] - subcoord[l, 1:3])^2)) / radius
                subplotList[[paste0("vert_", i, "_line_", l)]] <- plotly::plot_ly(
                    name = colnames(simMat)[i],
                    line, x = ~x, y = ~y, z = ~z,
                    type = "scatter3d", mode = "lines",
                    line = list(width = arrowLinewidth),
                    color = I(labelColors[i]),
                    showlegend = FALSE,
                    legendgroup = legendGroup,
                    hovertemplate = sprintf("Potential to %s: %.3f", colnames(simMat)[i], length)
                )
            }
            # Back to adding cones, we can add all cones for one vertex at a time
            # A cone should be anchored by its end ("tip")
            arrowEnd <- data.frame(
                x = subcoord[, 4],
                y = subcoord[, 5],
                z = subcoord[, 6]
            )
            arrowStart <- data.frame(
                x = subcoord[, 1],
                y = subcoord[, 2],
                z = subcoord[, 3]
            )
            EndToStartVec <- as.matrix(arrowStart - arrowEnd)
            scaleFac <- sqrt(rowSums(EndToStartVec^2))
            coneTail <- arrowEnd + EndToStartVec / scaleFac# * arrowLen
            coneDirec <- stats::setNames(arrowEnd - coneTail, c("u", "v", "w"))
            cone <- cbind(arrowEnd, coneDirec)
            hoverDF <- data.frame(
                sqrt(rowSums((subcoord[, 4:6] - subcoord[, 1:3])^2)) / radius
            )
            colnames(hoverDF) <- paste0("Potential to ", colnames(simMat)[i])
            subplotList[[paste0("vert_", i, "_cone")]] <- plotly::plot_ly(
                name = colnames(simMat)[i],
                cone, x = ~x, y = ~y, z = ~z, u = ~u, v = ~v, w = ~w,
                type = "cone", hoverinfo = "skip",
                # Set colorscale so that the cone is of consistent color
                colorscale = list(c(0, 1), c(labelColors[i], labelColors[i])),
                showscale = FALSE, anchor = "tail", sizeref = arrowLen*3,
                customdata = hoverDF,
                hovertemplate = .plotly_formatCustomData(hoverDF),
                hoverlabel = list(bgcolor = labelColors[i]),
                showlegend = TRUE,
                legendgroup = legendGroup,
                legendgrouptitle = if (i == 1) list(text = .htmlbold("Velocity to")) else NULL
            )
        }
    }

    # Show dots and tetra in the same plot
    final <- plotly::subplot(subplotList)
    # Adding global settings
    final <- plotly::layout(
        p = final,
        scene = list(
            xaxis = ax, yaxis = ax, zaxis = ax,
            showaxislabels = FALSE, showbackground = FALSE,
            showgrid = FALSE, showline = FALSE, showticklabels = FALSE
        ),
        hoverlabel = list(align = "left"),
        legend =  list(
            title = list(
                text = .htmlbold(colorArg$palette$legendTitle)
            ),
            itemsizing = "constant",
            y = 0.5
        )
    )
    # Dirty workaround to remove unreasonable NA entries in layout, that always
    # cause a warning
    final$x$layout <- final$x$layout[grep('NA', names(final$x$layout), invert = TRUE)]

    return(final)
}

# Workaround for plotly bug that it can't handle multi-variable customdata
# So directly template what to be shown be for passing to `hovertemplate`
.plotly_formatCustomData <- function(x) {
    args <- list()
    for (i in seq_len(ncol(x))) {
        if (is.numeric(x[,i])) val <- format(round(x[,i], 3), nsmall = 3)
        else val <- x[,i]
        collist <- list(
            colnames(x)[i], ": ", val, "\n"
        )
        args <- c(args, collist)
    }
    do.call(paste0, args)
}

#' Create GIF image for dynamic rotating view of 3D quaternary simplex plot
#' @param x Input object that \code{\link{plotQuaternary}} accepts.
#' @param ... All other arguments needed for \code{\link{plotQuaternary}}. Must
#' be specified with exact argument names instead of a positional manner.
#' @param cluster One cluster that exists in \code{clusterVar} if users
#' need to view the plot for a specific group. Default \code{NULL} plot all
#' cells.
#' @param filename Output GIF image file path. Default \code{NULL} does not
#' write to file.
#' @param fps Number of frame per second, must be a factor of 100. Default
#' \code{10}.
#' @param degreePerFrame Number of degree that the tetrahedron is rotated per
#' frame. Default \code{10}.
#' @param width,height,res \code{grDevices::\link[grDevices]{png}} parameters to
#' set figure size and resolutation. Width and Height are in inches. Default
#' \code{5}, \code{5}, \code{100}.
#' @return A object of class \code{magick-image} that can be shown in the Viewer
#' panel in RStudio or equivalent display device. If \code{filename} is
#' specified, the GIF image will be written to the file path.
#' @export
#' @examples
#' gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS", "CH", "ORT"))
#' \donttest{
#' writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster, features = gene,
#'                    vertices = c("RE", "OS", "CH", "ORT"),
#'                    gifPath = tempfile(fileext = ".gif"))
#' }
writeQuaternaryGIF <- function(
        x,
        ...,
        cluster = NULL,
        filename = NULL,
        fps = 10,
        degreePerFrame = 10,
        width = 5,
        height = 5,
        res = 100
) {
    if (!requireNamespace("magick", quietly = TRUE)) {
        stop("Package 'magick' must be installed for creating GIF figure.\n", # nocov
             "See https://cran.r-project.org/web/packages/magick/vignettes/", # nocov
             "intro.html#Installing_magick for detail.") # nocov
    }
    if (100 %% fps != 0) {
        cli::cli_abort("{.var FPS} must be a factor of 100.")
    }
    allImgPath <- c()

    methodArgs <- list(...)
    if (is.null(names(methodArgs)) || "" %in% names(methodArgs)) {
        cli::cli_abort(
            "Please set {.var ...} arguments with explicit argument names"
        )
    }
    ignore <- c("theta", "byCluster", "returnData", "interactive", "veloMat")
    if (any(names(methodArgs) %in% ignore)) {
        cli::cli_alert_warning(
            "These arguments are ignored in the GIF generation: {.val {names(methodArgs)[names(methodArgs) %in% ignore]}}"
        )
    }
    methodArgs[ignore] <- NULL

    plotData <- do.call(plotQuaternary, c(list(x = x, returnData = TRUE),
                                          methodArgs))
    clusterVar <- plotData$originalCluster
    methodArgs$colorArg <- plotData$colorArg
    if (!is.null(cluster)) {
        if (length(cluster) != 1)
            cli::cli_abort("Can only generate GIF for one cluster at a time")
        if (!cluster %in% clusterVar)
            cli::cli_abort(
                c("x" = "Specified {.var cluster} {.val {cluster}} is not availeble in the data",
                  "i" = "Available ones are: {.val {levels(clusterVar)}}")
            )
        simMat <- plotData$sim[clusterVar == cluster,]
        veloMat <- plotData$velo[clusterVar == cluster,]
        if (!"title" %in% names(methodArgs)) methodArgs$title <- cluster
        if ("colorArg" %in% names(methodArgs)) methodArgs$colorArg <- methodArgs$colorArg[clusterVar == cluster]
    } else {
        simMat <- plotData$sim
        veloMat <- plotData$velo
    }

    allTheta <- seq(0, 360, degreePerFrame)
    cli::cli_progress_bar(
        name = "Rotating",
        total = length(allTheta),
        status = "0 degree",
        type = "iterator"
    )
    for (i in seq_along(allTheta)) {
        theta <- allTheta[i]
        tmpfile <- tempfile(pattern = "quaternary_", fileext = ".png")
        allImgPath <- c(allImgPath, tmpfile)
        p <- do.call(plotQuaternary,
                     c(list(x = simMat, veloMat = veloMat, theta = theta,
                            interactive = FALSE),
                       methodArgs))

        grDevices::png(tmpfile, width = width, height = height, res = res, unit = "in")
        print(p)
        grDevices::dev.off()
        cli::cli_progress_update(set = i, status = paste0(theta, " degree"))
    }
    cli::cli_progress_done()
    cli::cli_process_start(msg = "Assembling frames")
    imgList <- lapply(allImgPath, magick::image_read)
    imgJoined <- magick::image_join(imgList)
    imgAnimated <- magick::image_animate(imgJoined, fps = fps)
    if (!is.null(filename)) {
        filename <- normalizePath(filename, mustWork = FALSE)
        magick::image_write(image = imgAnimated,
                            path = filename)
        cli::cli_process_done(msg_done = "GIF written to file: {.file {filename}}")
    } else {
        cli::cli_process_done()
    }

    return(imgAnimated)
}
