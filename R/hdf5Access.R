#' Extract a variable from adata.obs stored in an H5AD file
#' @description
#' Primarily designed for fetching the annotation used for visualization.
#' @param filename File path to the H5AD file.
#' @param obsKey The variable name to extract, must use only one character
#' string.
#' @param named Logical, whether to name the vector with cell IDs that came from
#' \code{adata.obs_names}. Default \code{TRUE}.
#' @param categoricalAsFactor Logical, whether to convert categorical variables
#' to factors. Default \code{TRUE}.
#' @return A vector of the extracted variable, or a factor if the variable is
#' encoded to be categorical and \code{categoricalAsFactor = TRUE}.
#' @export
#' @family H5AD-reader
#' @examples
#' \dontrun{
#' h5adFile <- "path/to/analysis.h5ad"
#' cluster <- readH5ADObsVar(h5adFile, "leiden")
#' }
readH5ADObsVar <- function(
        filename,
        obsKey,
        named = TRUE,
        categoricalAsFactor = TRUE
) {
    adataObj <- .openH5AD(filename, ".h5ad")
    obs <- adataObj[['obs']]
    idsCol <- hdf5r::h5attr(obs, "_index")
    ids <- obs[[idsCol]][]
    if (length(obsKey) != 1 || !is.character(obsKey)) {
        cli::cli_abort("{.field obsKey} must be a single character string.")
    }
    if (!obsKey %in% names(obs)) {
        cli::cli_abort("obs key {.val {obsKey}} not found in {.code adata.obs.columns}.")
    }
    variable <- obs[[obsKey]]

    encodingType <- hdf5r::h5attr(variable, "encoding-type")
    if (encodingType == "array") {
        value <- variable[]
    } else if (encodingType == "string-array") {
        value <- variable[]
    } else if (encodingType == "categorical") {
        categories <- variable[['categories']][]
        # +1 to address Python 0-based indexing
        codes <- variable[['codes']][] + 1
        value <- categories[codes]
        if (isTRUE(categoricalAsFactor)) {
            value <- factor(value, levels = categories)
        }
    } else {
        cli::cli_abort("Unsupported encoding type {.val {encodingType}}.")
    }
    if (isTRUE(named)) names(value) <- ids
    adataObj$close_all()
    return(value)
}

#' Extract a sparse matrix from adata.uns stored in an H5AD file
#' @description
#' Primarily designed for fetching the velocity data presented as a cell-cell
#' transition graph.
#' @inheritParams readH5ADObsVar
#' @param unsKey The \code{adata.uns} key to extract, must use only one
#' character string.
#' @return A CSC-matrix of "dgCMatrix" class
#' @export
#' @family H5AD-reader
#' @examples
#' \dontrun{
#' h5adFile <- "path/to/analysis.h5ad"
#' velo <- readH5ADUnsSpMat(h5adFile, "velo_s_norm_graph")
#' }
readH5ADUnsSpMat <- function(
        filename,
        unsKey
) {
    adataObj <- .openH5AD(filename, ".h5ad")
    uns <- adataObj[['uns']]
    if (length(unsKey) != 1 || !is.character(unsKey)) {
        cli::cli_abort("{.field unsKey} must be a single character string.")
    }
    if (!unsKey %in% names(uns)) {
        cli::cli_abort("uns key {.val {unsKey}} not found in {.code adata.uns.keys()}.")
    }
    spMat <- uns[[unsKey]]
    encodingType <- hdf5r::h5attr(spMat, "encoding-type")
    if (encodingType != "csr_matrix") {
        cli::cli_abort("This function only extracts {.val csr-matrix} encoded data.")
    }
    i <- spMat[['indices']][] + 1
    p <- spMat[['indptr']][]
    x <- spMat[['data']][]
    # Returning csc-matrix from csr-matrix transposes the matrix
    dims <- rev(hdf5r::h5attr(spMat, "shape"))
    spMat <- Matrix::sparseMatrix(i = i, p = p, x = x, dims = dims)
    adataObj$close_all()
    return(spMat)
}

#' Extract `adata.obs_names` from an H5AD file
#' @description
#' It frequently happens that velocity analyses stored in H5AD files do not
#' contain the full raw count data suggested for CytoSimplex visualization.
#' Extracting the cell IDs (e.g. barcodes) helps matching the velocity data to
#' raw count data imported from other sources.
#' @inheritParams readH5ADObsVar
#' @return A character vector of cell IDs.
#' @export
#' @family H5AD-reader
#' @examples
#' \dontrun{
#' h5adFile <- "path/to/analysis.h5ad"
#' cellIDs <- readH5ADObsNames(h5adFile)
#' }
readH5ADObsNames <- function(
        filename
) {
    adataObj <- .openH5AD(filename, '.h5ad')
    obs <- adataObj[['obs']]
    idsCol <- hdf5r::h5attr(obs, "_index")
    ids <- obs[[idsCol]][]
    adataObj$close_all()
    return(ids)
}

#' Extract the raw counts from a LOOM file
#' @description
#' This function is primarily designed for fetching the raw count data from a
#' LOOM file, output by \href{https://velocyto.org/}{Velocyto}. We by default
#' use the spliced counts.
#' @details
#' The velocyto output LOOM file is HDF5 based and is roughly organized as
#' follows:
#' \itemize{
#' \item{\code{"matrix"}: The whole raw counts, which is the sum of spliced, unspliced
#' and ambiguous counts.}
#' \item{layers: A group like a folder
#'   \itemize{
#'   \item{\code{"layers/spliced"}: The spliced counts.}
#'   \item{\code{"layers/unspliced"}: The unspliced counts.}
#'   \item{\code{"layers/ambiguous"}: The ambiguous counts.}
#'   }
#' }
#' }
#'
#' An AnnData object created with Scanpy by default loads the data with a
#' different structure, so that all the four matrices are accessible in
#' \code{adata.layers} and set one of them (by default \code{"layers/spliced"})
#' to \code{adata.X}.
#'
#' @param filename File path to the LOOM file.
#' @param matrixPath A path in the LOOM file to the matrix to extract, following
#' the inner HDF5 structure. Default \code{"layers/spliced"}. See Details.
#' @param cellID The name of the cell ID column in the LOOM column-attributes.
#' The same thing as argument \code{obs_names} of \code{scanpy.read_loom}.
#' Default \code{"CellID"}.
#' @param featureID The name of the feature ID column in the LOOM
#' row-attributes. The same thing as argument \code{var_names} of
#' \code{scanpy.read_loom}. Default \code{"Gene"}.
#' @param chunkSize The maximum size of the chunk to load the matrix. Default
#' 1000.
#' @return A sparse matrix of class "dgCMatrix", with cells as columns and genes
#' as rows.
#' @export
#' @family H5AD-reader
#' @examples
#' \dontrun{
#' loomFile <- "velocyto/out/analysis.loom"
#' rawCounts <- readVelocytoLoom(loomFile)
#' }
readVelocytoLoom <- function(
        filename,
        matrixPath = "layers/spliced",
        cellID = "CellID",
        featureID = "Gene",
        chunkSize = 1000
) {
    loom <- .openH5AD(filename, ".loom")
    cellIDs <- loom[[file.path("col_attrs", cellID)]][]
    featureIDs <- loom[[file.path("row_attrs", featureID)]][]
    matH5D <- loom[[matrixPath]]
    chunkDims <- matH5D$chunk_dims
    cellChunkDims <- chunkDims[1]
    chunkSize <- chunkSize - chunkSize%%cellChunkDims
    nChunks <- ceiling(length(cellIDs) / chunkSize)
    spMat <- NULL
    cli::cli_progress_bar(name = "Loading from LOOM", total = nChunks)
    for (i in seq_len(nChunks)) {
        start <- (i - 1) * chunkSize + 1
        end <- min(i * chunkSize, length(cellIDs))
        chunkMat <- matH5D[start:end, ]
        chunkMat <- t(chunkMat)
        chunkMat <- methods::as(chunkMat, "CsparseMatrix")
        spMat <- cbind(spMat, chunkMat)
        cli::cli_progress_update(set = i)
    }
    cli::cli_process_done()
    dimnames(spMat) <- list(featureIDs, cellIDs)
    return(spMat)
}


.openH5AD <- function(filename, checkExt) {
    if (!file.exists(filename)) {
        cli::cli_abort("{checkExt} File not found: {.file {filename}}")
    }
    if (!endsWith(filename, toupper(checkExt)) &&
        !endsWith(filename, tolower(checkExt))) {
        cli::cli_alert_warning("File extension is not {.var {checkExt}}, opening anyway.")
    }
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg hdf5r} is required for extracting data from H5AD files.")
    }
    return(hdf5r::H5File$new(filename, mode = "r"))
}


