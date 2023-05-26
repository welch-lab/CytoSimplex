#' Normalize each column of the input matrix by the column sum
#' @param x Feature by observation matrix.
#' @param scaleFactor Multiplier on normalized data. Default \code{NULL}.
#' @param log Logical. Whether to take log1p transformation after scaling.
#' Default \code{FALSE}
#' @param ... Additional arguments passed to methods
#' @return Normalized matrix of the same size
#' @export
#' @rdname colNormalize
#' @examples
#' rnaNorm <- colNormalize(rnaRaw)
colNormalize <- function(x, scaleFactor = NULL, log = FALSE, ...) {
    UseMethod("colNormalize", x)
}

#' @rdname colNormalize
#' @export
#' @method colNormalize default
colNormalize.default <- function(x, scaleFactor = NULL, log = FALSE, ...) {
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

#' @rdname colNormalize
#' @param slot For "Seurat" method, choose from \code{"counts"},
#' \code{"data"} or \code{"scale.data"}. Default \code{"counts"}.
#' @param assay For "Seurat" method, the specific assay to get data from.
#' Default \code{NULL} to the default assay.
#' @return A Seurat object with normalized data in the specified slot of the
#' specified assay.
#' @export
#' @method colNormalize Seurat
#' @examples
#'
#' # Seurat example
#' if (FALSE) {
#'     library(Seurat)
#'     srt <- CreateSeuratObject(rnaRaw)
#'     srt <- colNormalize(srt)
#' }
colNormalize.Seurat <- function(
    x,
    scaleFactor = NULL,
    log = FALSE,
    assay = NULL,
    slot = "counts",
    ...
) {
    value <- .getSeuratData(x, features = NULL, assay = NULL, slot = slot,
                            clusterVar = NULL)
    mat <- value[[1]]
    norm <- colNormalize.default(mat, scaleFactor = scaleFactor, log = log)
    x <- Seurat::SetAssayData(x, slot = "data", assay = assay, new.data = norm)
    return(x)
}

#' @rdname colNormalize
#' @param assay.type For "SingleCellExperiment" method, the assay type to get
#' data from. Default \code{"counts"}.
#' @return A SingleCellExperiment object with normalized data in the specified
#' assay. \code{"normcounts"} if \code{log = FALSE} and \code{"logcounts"} if
#' \code{log = TRUE}.
#' @export
#' @method colNormalize SingleCellExperiment
#' @examples
#'
#' # SingleCellExperiment example
#' if (FALSE) {
#'     library(SingleCellExperiment)
#'     sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#'     sce <- colNormalize(sce)
#' }
colNormalize.SingleCellExperiment <- function(
    x,
    scaleFactor = NULL,
    log = FALSE,
    assay.type = "counts",
    ...
) {
    value <- .getSCEData(x, subset.row = NULL, clusterVar = NULL,
                         assay.type = assay.type)
    mat <- value[[1]]
    norm <- colNormalize.default(mat, scaleFactor = scaleFactor, log = log)
    if (isTRUE(log)) SingleCellExperiment::logcounts(x) <- norm
    else SingleCellExperiment::normcounts(x) <- norm

    return(x)
}
