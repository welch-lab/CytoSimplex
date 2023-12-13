#' Normalize each column of the input matrix by the column sum
#' @param x Feature by observation matrix. Alternatively, \code{Seurat} object
#' or \code{SingleCellExperiment} object with raw counts available are also
#' supported.
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
    dn <- dimnames(x)
    x <- colNormalize_dense(x, base::colSums(x))
    dimnames(x) <- dn
    if (!is.null(scaleFactor)) x <- x * scaleFactor
    if (isTRUE(log)) x <- log1p(x)
    return(x)
}

#' @rdname colNormalize
#' @export
#' @method colNormalize dgCMatrix
colNormalize.dgCMatrix <- function(x, scaleFactor = NULL, log = FALSE, ...) {
    x@x <- x@x / rep.int(Matrix::colSums(x), diff(x@p))
    if (!is.null(scaleFactor)) x <- x * scaleFactor
    if (isTRUE(log)) x <- log1p(x)
    return(x)
}

#' @rdname colNormalize
#' @param assay For "Seurat" method, the specific assay to get data from.
#' Default \code{NULL} to the default assay.
#' @param layer For "Seurat" method, which layer of the assay to be used.
#' Default \code{"counts"}.
#' @return A Seurat object with normalized data in the specified slot of the
#' specified assay.
#' @export
#' @method colNormalize Seurat
#' @examples
#' \donttest{
#' # Seurat example
#' library(Seurat)
#' srt <- CreateSeuratObject(rnaRaw)
#' srt <- colNormalize(srt)
#' }
colNormalize.Seurat <- function(
    x,
    scaleFactor = NULL,
    log = FALSE,
    assay = NULL,
    layer = "counts",
    ...
) {
    value <- .getSeuratData(x, assay = assay, layer = layer, clusterVar = NULL)
    mat <- value[[1]]
    norm <- colNormalize(mat, scaleFactor = scaleFactor, log = log)
    SeuratObject::LayerData(x, layer = "data", assay = assay) <- norm
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
#' \donttest{
#' # SingleCellExperiment example
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#' sce <- colNormalize(sce)
#' }
colNormalize.SingleCellExperiment <- function(
    x,
    scaleFactor = NULL,
    log = FALSE,
    assay.type = "counts",
    ...
) {
    value <- .getSCEData(x, clusterVar = NULL, assay.type = assay.type)
    mat <- value[[1]]
    norm <- colNormalize(mat, scaleFactor = scaleFactor, log = log)
    if (isTRUE(log)) SingleCellExperiment::logcounts(x) <- norm
    else SingleCellExperiment::normcounts(x) <- norm

    return(x)
}
