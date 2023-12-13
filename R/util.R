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
            stop("Must specify ", n, " different vertices.")
        }
        if (length(vertices) > n) {
            warning(n, " vertices are expected while ", length(vertices),
                    " are specified. Using the first ", n, '.')
            vertices <- vertices[seq_len(n)]
        }
    }

    if (length(clusterVar) != ncol(object)) {
        stop("Length of `clusterVar` must match `ncol(object)`.")
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

    if (is.null(clusterVar)) {
        if (inherits(object, "SingleCellExperiment")) {
            clusterVar <- SingleCellExperiment::colLabels(object)
        }
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
