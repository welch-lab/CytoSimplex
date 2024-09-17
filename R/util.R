`%||%` <- function(x, y) {
    if (is.null(x)) y
    else x
}

.scaleMinMax <- function(x) {
    if (all(x == 0)) return(x)
    else {
        x <- (x - min(x, na.rm = TRUE)) /
            (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        return(x)
    }
}

.normalize <- function(x) {
    x <- x + 1e-8
    x / sum(x, na.rm = TRUE)
}

is.rawCounts <- function(x) {
    if (inherits(x, "dgCMatrix")) {
        is_rawCounts_sparse(x)
    } else if (is.matrix(x)) {
        is_rawCounts_dense(x)
    } else {
        FALSE
    }
}

.checkVertex <- function(
        object,
        clusterVar,
        vertices,
        n
) {
    if (!is.null(n)) {
        if (length(vertices) < n) {
            cli::cli_abort(
                c("x" = "Must specify {n} different vertices.",
                  "i" = "Only {length(vertices)} are specified.")
            )
        }
        if (length(vertices) > n) {
            cli::cli_alert_warning(
                "{n} vertices are expected while {length(vertices)} are specified. Using the first {n}."
            )
            vertices <- vertices[seq_len(n)]
        }
    }

    if (length(clusterVar) != ncol(object)) {
        cli::cli_abort(
            c("x" = "Length of {.var clusterVar} must be {ncol(object)}.",
              "i" = "{length(clusterVar)} is provided.")
        )
    }

    allVClust <- vertices
    if (is.list(vertices)) {
        allVClust <- unlist(vertices)
        if (length(allVClust) != length(unique(allVClust))) {
            stop("Overlap found between elements in list vertex specification.")
        }
    }
    if (!all(allVClust %in% clusterVar)) {
        stop("Specified vertex clusters are not all found in the cluster ",
             "variable")
    }

    if (is.list(vertices)) {
        clusterVar <- as.character(clusterVar)
        for (v in names(vertices)) {
            clusterVar[clusterVar %in% vertices[[v]]] <- v
        }
        return(list(factor(clusterVar), names(vertices)))
    } else {
        return(list(clusterVar, vertices))
    }
}

.getSeuratData <- function(
        object,
        assay = NULL,
        layer = "counts",
        clusterVar = NULL
) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Please install package 'Seurat' before interacting with a ", # nocov
             "Seurat object.\ninstall.packages(\"Seurat\")") # nocov
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Please install package 'SeuratObject' before interacting with a ", # nocov
             "Seurat object.\ninstall.packages(\"Seurat\")") # nocov
    }
    mat <- SeuratObject::LayerData(object, layer = layer, assay = assay)
    clusterVar <- clusterVar %||% SeuratObject::Idents(object)
    if (length(clusterVar) == 1) {
        clusterVar <- object[[clusterVar]][[1]]
    }
    if (length(clusterVar) != ncol(object)) {
        stop("Invalid `clusterVar`.")
    }
    return(list(mat, clusterVar))
}

.getSCEData <- function(
        object,
        clusterVar = NULL,
        assay.type = "logcounts"
) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Please install package 'SingleCellExperiment' before ", # nocov
             "interacting with a SingleCellExperiment object.", # nocov
             "\nBiocManager::install(\"SingleCellExperiment\")") # nocov
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Please install package 'SummarizedExperiment' before ", # nocov
             "interacting with a SingleCellExperiment object.", # nocov
             "\nBiocManager::install(\"SummarizedExperiment\")") # nocov
    }
    mat <- SummarizedExperiment::assay(object, assay.type)

    clusterVar <- clusterVar %||% SingleCellExperiment::colLabels(object)
    if (length(clusterVar) == 0) {
        clusterVar <- NULL
    } else if (length(clusterVar) == 1) {
        clusterVar <- SummarizedExperiment::colData(object)[[clusterVar]]
    } else if (length(clusterVar) != ncol(object)) {
        stop("Invalid `clusterVar`.")
    }
    return(list(mat, clusterVar))
}

# .ligerPrepare <- function(
#         object,
#         clusterVar,
#         features = NULL,
#         useDatasets = NULL
# ) {
#     if (!requireNamespace("rliger2", quietly = TRUE)) {
#         stop("Please install package \"rliger\".")
#     }
#     # rliger2:::.checkUseDataset
#     if (is.null(useDatasets)) {
#         useDatasets <- names(object)
#     } else {
#         if (is.numeric(useDatasets)) {
#             if (max(useDatasets) > length(object)) {
#                 stop("Numeric dataset index out of bound. Only ",
#                      length(object), " datasets exist.")
#             }
#             useDatasets <- unique(useDatasets)
#             useDatasets <- names(object)[useDatasets]
#         } else if (is.logical(useDatasets)) {
#             if (length(useDatasets) != length(object)) {
#                 stop("Logical dataset subscription does not match the number ",
#                      "of datasets (", length(object), ").")
#             }
#             useDatasets <- names(object)[useDatasets]
#         } else if (is.character(useDatasets)) {
#             if (any(!useDatasets %in% names(object))) {
#                 notFound <- useDatasets[!useDatasets %in% names(object)]
#                 stop("Specified dataset name(s) not found: ",
#                      paste(notFound, collapse = ", "))
#             }
#         } else {
#             stop("Please use a proper numeric/logical/character vector to ",
#                  "select dataset to use.")
#         }
#     }
#
#     matList <- rliger2::getMatrix(object, slot = "normData",
#                                   dataset = useDatasets, returnList = TRUE)
#     mat <- rliger2::mergeSparseAll(matList)
#     if (!is.null(features)) {
#         if (!all(features %in% rownames(mat))) {
#             nf <- features[!features %in% rownames(mat)]
#             warning("Following specified features not found in the union of ",
#                     "selected datasets: ", paste(nf, collapse = ", "))
#             features <- features[features %in% rownames(mat)]
#         }
#         if (length(features) > 1) mat <- mat[features,]
#         else {
#             stop("Too few specified features available in selected datasets.")
#         }
#     }
#     if (length(clusterVar) == 1) {
#         clusterVar <- rliger2::cellMeta(
#             object, columns = clusterVar,
#             cellIdx = object$dataset %in% useDatasets
#         )
#     }
#     return(list(mat, clusterVar))
# }

# Emulate ggplot categorical colors
ggColor <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

continuousPaletteName <- stats::setNames(
    object = c(rep("viridis", 16),
               rep("RcolorBrewer", sum(RColorBrewer::brewer.pal.info$category != "qual"))),
    nm = c(
        "magma", "A", "inferno", "B", "plasma", "C", "viridis", "D",
        "cividis", "E", "rocket", "F", "mako", "G", "turbo", "H",
        rownames(RColorBrewer::brewer.pal.info)[RColorBrewer::brewer.pal.info$category != "qual"]
    )
)

categoricalPaletteName <- stats::setNames(
    object = rep("RColorBrewer", sum(RColorBrewer::brewer.pal.info$category == "qual")),
    nm = rownames(RColorBrewer::brewer.pal.info)[RColorBrewer::brewer.pal.info$category == "qual"]
)

# This helper function is mainly for checking if the dot color specficication
# is reasonable, and return properly prepared information that can be used by
# either ggplot2, plot3D or plotly.
#
# Possible scenarios
# 1. Setting all dots to one color
#    Return: type: "single", colorBy: NULL, palette: colors: one color, colors: n colors
# 2. colorBy is categorical, colors are specified for each of the levels
#    Return: type: "categorical", colorBy: factor,
#            palette: colors: nlevels(factor) colors, labels: levels(factor)
#            colors: n colors
# 3. colorBy is continuous, colors are to be determined by viridis (opt+dir)
# or customized palette (colors+breaks)
#    Return:
#      - type: "continuous",
#      - colorBy: numeric,
#      - palette:
#        - colors: viridis colors/or given,
#        - breaks: inferred / or given,
#      - colors: n colors
# 4. colorBy is not given, colors are explicitly given for each dot. Can't infer for annotation but need to set each color
#    Return: type: "customize", colorBy: NULL, palette: NULL, colors: n colors
.checkDotColor <- function(
        n,
        cluster,
        colorBy = NULL,
        colors = NULL,
        palette = "D",
        direction = 1,
        breaks = NULL,
        legendTitle = NULL
) {
    colorArg <- list(
        type = NULL,
        colorBy = colorBy,
        palette = list(
            colors = NULL,
            labels = NULL,
            breaks = breaks,
            palette = palette,
            dep = NULL,
            direction = direction,
            legendTitle = legendTitle
        ),
        colors = NULL
    )
    if (is.null(colorBy)) {
        if (is.null(colors)) {
            # Nothing specified, use cluster for coloring and generate ggplot default
            # S2
            legendTitle <- legendTitle %||% "cluster"
            colorArg <- .genCategCol(n = n, palette = palette, colorBy = cluster, legendTitle = legendTitle)
        } else {
            # No grouping with colors but have colors directly specified
            # Can be one color for all dots, or nlevels colors for each cluster,
            # or n colors for explicit setting

            # check lenght of colors
            if (length(colors) != 1 &&
                length(colors) != nlevels(cluster) &&
                length(colors) != n) {
                cli::cli_abort("Length of {.var dotColor} must be 1 for all dots, {nlevels(cluster)} for each cluster, or {n} for explicit setting.")
            }
            if (length(colors) == 1) {
                # S1
                colorArg$type <- "single"
                colorArg$colors <- rep(colors, n)
                colorArg$palette$colors <- colors
            }
            if (length(colors) == nlevels(cluster)) {
                # S2
                colorArg$type <- "categorical"
                colorArg$colorBy <- cluster
                colorArg$colors <- colors[cluster]
                names(colorArg$colors) <- names(cluster)
                legendTitle <- legendTitle %||% "cluster"
                colorArg$palette <- list(colors = colors, labels = levels(cluster), legendTitle = legendTitle)
            }
            if (length(colors) == n) {
                # S4
                colorArg$type <- "customize"
                colorArg$colors <- colors
            }
        }
    } else {
        # ColorBy is given now
        if (length(colorBy) != n) {
            cli::cli_abort(
                c("x" = "Length of {.var dotColorBy} must be {n}.",
                  "i" = "{length(colorBy)} is provided.")
            )
        }
        # Hack for avoiding undefined variable
        dotColorBy <- NULL
        defaultTitle <- deparse(substitute(dotColorBy, rlang::caller_env()))
        legendTitle <- legendTitle %||% defaultTitle
        colorArg$colorBy <- colorBy
        if (is.character(colorBy)) colorBy <- factor(colorBy)
        if (is.factor(colorBy)) {
            colorArg$type <- "categorical"
            colorArg$palette$labels <- levels(colorBy)
            colorArg$palette$legendTitle <- legendTitle %||% "cluster"
            if (is.null(colors)) {
                colorArg <- .genCategCol(n = n, palette = palette, colorBy = colorBy, legendTitle = colorArg$palette$legendTitle)
            } else {
                if (length(colors) < nlevels(colorBy)) {
                    cli::cli_abort("{.var dotColor} must have as least {nlevels(colorBy)} values.")
                }
                if (length(colors) > nlevels(colorBy)) {
                    cli::cli_alert_warning(
                        "{nlevels(colorBy)} colors expected from {.var dotColor} but {length(colors)} are provided. Using the first {nlevels(colorBy)}.",
                    )
                }
                colors <- colors[seq_len(nlevels(colorBy))]
                colorArg$palette$colors <- colors
            }
            colorArg$colors <- colorArg$palette$colors[colorBy]
        } else {
            colorArg$type <- "continuous"
            colorArg$palette$legendTitle <- legendTitle %||% "value"
            if (is.null(colors)) {
                colorArg$palette$breaks <- breaks
                # No customized color palettes given, use palette
                if (!palette %in% names(continuousPaletteName)) {
                    cli::cli_abort(
                        c("x" = "Specified {.var palette} {.val {palette}} is not supported",
                          "i" = "Available ones include options from {.pkg viridis} and {.pkg RColorBrewer} (type of 'div' or 'seq').",
                          "i" = "Namingly: {.val {names(continuousPaletteName)}}.")
                    )
                }
                dep <- continuousPaletteName[palette]
                colorArg$palette$palette <- palette
                colorArg$palette$dep <- dep
                if (dep == "viridis") {
                    colorArg$palette$colors <- viridis::viridis(100, option = palette, direction = direction)
                } else if (dep == "RcolorBrewer") {
                    maxcolor <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
                    baseColors <- RColorBrewer::brewer.pal(n = maxcolor, name = palette)
                    interpolate <- grDevices::colorRampPalette(baseColors)
                    if (direction == 1) {
                        colorArg$palette$colors <- interpolate(100)
                    } else if (direction == -1) {
                        colorArg$palette$colors <- rev(interpolate(100))
                    } else {
                        cli::cli_abort("{.var direction} must be either 1 or -1.")
                    }
                }
                colorArg$colors <- colorArg$palette$colors[round(.scaleMinMax(colorBy)*99) + 1]
            } else {
                # TODO: still need some consideration
                if (length(colors) < 2) {
                    cli::cli_abort(
                        c("x" = "Length of {.var dotColor} must be at least 2.",
                          "i" = "{length(colors)} is provided.")
                    )
                }
                colorArg$palette$colors <- colors
            }
        } # !is.factor(colorBy) ends
    } # !is.null(colorBy) ends
    class(colorArg) <- "colorArg"
    return(colorArg)
}

#' Subset color setting arguments. Internal use only
#' @method `[` colorArg
#' @param x colorArg object
#' @param i subscriber
#' @param j,... Ignored
#' @param drop Whether to drop factor levels
#' @export
#' @return Subset colorArg
#' @noRd
`[.colorArg` <- function(x, i, j, ..., drop = FALSE) {
    if (x$type == "categorical") {
        x$colorBy <- x$colorBy[i, drop = drop]
    }
    if (x$type == "continuous") {
        x$colorBy <- x$colorBy[i]
    }
    x$colors <- x$colors[i]
    x
}

.genCategCol <- function(n, palette, colorBy, legendTitle) {
    colorArg <- list(
        type = "categorical",
        colorBy = colorBy,
        palette = list(
            colors = NULL,
            labels = levels(colorBy),
            palette = palette,
            dep = NULL,
            legendTitle = legendTitle
        ),
        colors = NULL
    )
    if (palette %in% names(categoricalPaletteName)) {
        dep <- categoricalPaletteName[palette]
        colorArg$palette$dep <- dep
        if (dep != "RColorBrewer") { # nocov start
            stop("Only RColorBrewer palettes are supported for categorical coloring.")
        } # nocov end
        maxcolor <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
        if (nlevels(colorArg$colorBy) > maxcolor) {
            cli::cli_alert_danger(
                "Palette {.val {palette}} has only {maxcolor} colors, but {nlevels(colorBy)} levels are present."
            )
            cli::cli_alert_info(
                "Using ggplot default categorical colors instead."
            )
            colorArg$palette$colors <- ggColor(nlevels(colorArg$colorBy))
        } else {
            baseColors <- RColorBrewer::brewer.pal(n = maxcolor, name = palette)
            colorArg$palette$colors <- baseColors[seq_len(nlevels(colorArg$colorBy))]
        }
    } else {
        colorArg$palette$palette <- NULL
        colorArg$palette$dep <- NULL
        colorArg$palette$colors <- ggColor(nlevels(colorArg$colorBy))
    }
    colorArg$colors <- colorArg$palette$colors[colorArg$colorBy]
    return(colorArg)
}


.htmlbold <- function(x) {
    paste0("<b>", x, "</b>")
}

.colAlpha <- function(color, alpha = 1) {
    if (alpha < 0 || alpha > 1) {
        cli::cli_abort("{.var alpha} must be between 0 and 1.")
    }
    rgba <- grDevices::col2rgb(color, alpha = TRUE)
    grDevices::rgb(rgba[1], rgba[2], rgba[3], alpha = 255*alpha, maxColorValue = 255)
}

# Build plotly required colorscale used for "markers" attribute
# colorBy - is the numeric vector used for coloring
# colorUse - is a vector of colors to be used on the colorbar
# return format
# list(
#     list(value_0, "rgba(r, g, b, a)"),
#     list(value_1, "rgba(r, g, b, a)"),
#     ...
# )
.plotlyColorscale <- function(colorUse) {
    steps <- seq(0, 1, length.out = length(colorUse))
    colorRGB <- grDevices::col2rgb(colorUse)
    colorHTMLRGB <- paste0("rgb(", colorRGB[1,], ",", colorRGB[2,], ",", colorRGB[3,], ")")
    lapply(seq_along(colorUse), function(i) list(steps[i], colorHTMLRGB[i]))
}
