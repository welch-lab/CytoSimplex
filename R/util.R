.scaleMinMax <- function(x) {
    if (all(x == 0)) return(x)
    else {
        (x - min(x, na.rm = TRUE)) /
            (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    }
}

.normalize <- function(x) {
    x <- x + 1e-8
    x / sum(x, na.rm = TRUE)
}

#' Normalize each column of the input matrix by the column sum
#' @param x Feature by observation matrix.
#' @param scaleFactor Multiplier on normalized data. Default \code{NULL}.
#' @param log Logical. Whether to take log1p transformation after scaling.
#' Default \code{FALSE}
#' @return Normalized matrix of the same size
#' @export
#' @examples
#' rnaNorm <- colNormalize(rnaRaw)
colNormalize <- function(x, scaleFactor = NULL, log = FALSE) {
    if (inherits(x, "dgCMatrix")) {
        x@x <- x@x / rep.int(Matrix::colSums(x), diff(x@p))
    } else if (is.matrix(x)) {
        x <- colNormalize_dense(x, base::colSums(x))
    } else {
        stop("Input matrix of class ", class(x)[1], " is not yet supported.")
    }
    if (!is.null(scaleFactor)) x <- x * scaleFactor
    if (isTRUE(log)) x <- log1p(x)
    return(x)
}

.ligerPrepare <- function(
        object,
        clusterVar,
        features = NULL,
        useDatasets = NULL
) {
    if (!requireNamespace("rliger2", quietly = TRUE)) {
        stop("Please install package \"rliger\".")
    }
    # rliger2:::.checkUseDataset
    if (is.null(useDatasets)) {
        useDatasets <- names(object)
    } else {
        if (is.numeric(useDatasets)) {
            if (max(useDatasets) > length(object)) {
                stop("Numeric dataset index out of bound. Only ",
                     length(object), " datasets exist.")
            }
            useDatasets <- unique(useDatasets)
            useDatasets <- names(object)[useDatasets]
        } else if (is.logical(useDatasets)) {
            if (length(useDatasets) != length(object)) {
                stop("Logical dataset subscription does not match the number ",
                     "of datasets (", length(object), ").")
            }
            useDatasets <- names(object)[useDatasets]
        } else if (is.character(useDatasets)) {
            if (any(!useDatasets %in% names(object))) {
                notFound <- useDatasets[!useDatasets %in% names(object)]
                stop("Specified dataset name(s) not found: ",
                     paste(notFound, collapse = ", "))
            }
        } else {
            stop("Please use a proper numeric/logical/character vector to ",
                 "select dataset to use.")
        }
    }

    matList <- rliger2::getMatrix(object, slot = "normData",
                                  dataset = useDatasets, returnList = TRUE)
    mat <- rliger2::mergeSparseAll(matList)
    if (!is.null(features)) {
        if (!all(features %in% rownames(mat))) {
            nf <- features[!features %in% rownames(mat)]
            warning("Following specified features not found in the union of ",
                    "selected datasets: ", paste(nf, collapse = ", "))
            features <- features[features %in% rownames(mat)]
        }
        if (length(features) > 1) mat <- mat[features,]
        else {
            stop("Too few specified features available in selected datasets.")
        }
    }
    if (length(clusterVar) == 1) {
        clusterVar <- rliger2::cellMeta(
            object, columns = clusterVar,
            cellIdx = object$dataset %in% useDatasets
        )
    }
    return(list(mat, clusterVar))
}
