#' Pick top differentially presented features for similarity calculation
#' @description
#' Performs wilcoxon rank-sum test on input matrix. While \code{clusterVar} and
#' \code{vertices} together defines the groups of cells to be set as terminals
#' of the simplex, this function will test each of these groups against the rest
#' of the cells. The U-Statistics (\code{statistic}), p-value (\code{pval}) and
#' adjusted p-value (\code{padj}), together with average presence in group
#' (\code{avgExpr}), log fold-change (\code{logFC}), AUC (\code{auc}),
#' percentage in group (\code{pct_in}) and percentage out of group
#' (\code{pct_out}) will be calculated. Set \code{returnStats = TRUE} to return
#' the full statistics table.
#'
#' Top features are selected by sorting primarily on adjusted p-value, and
#' secondarily on log fold-change, after filtering for up-regulated features.
#' @param x Dense or sparse matrix, observation per column. Preferrably a raw
#' count matrix. Alternatively, a \code{Seurat} object or a
#' \code{SingleCellExperiment} object.
#' @param clusterVar A vector/factor assigning the cluster variable to each
#' column of the matrix object. For "Seurat" method, \code{NULL} (default) for
#' \code{Idents(x)}, or a variable name in \code{meta.data} slot. For
#' "SingleCellExperiment" method, \code{NULL} (default) for \code{colLabels(x)},
#' or a variable name in \code{colData} slot.
#' @param vertices Vector of cluster names that will be used for plotting. Or a
#' named list that groups clusters as a terminal vertex. There must not be any
#' overlap between groups.
#' @param ... Arguments passed to methods.
#' @export
#' @return When \code{returnStats = FALSE} (default), a character vector of at
#' most \code{length(unique(vertices))*nTop} feature names. When
#' \code{returnStats = TRUE}, a data.frame of wilcoxon rank sum test statistics.
#' @rdname selectTopFeatures
#' @examples
#' selectTopFeatures(rnaRaw, rnaCluster, c("OS", "RE"))
selectTopFeatures <- function(
    x,
    clusterVar,
    vertices,
    ...
) {
    UseMethod("selectTopFeatures", x)
}

#' @rdname selectTopFeatures
#' @export
#' @method selectTopFeatures default
#' @param nTop Number of top differentially presented features per terminal.
#' Default \code{30}.
#' @param processed Logical. Whether the input matrix is already processed.
#' \code{TRUE} will bypass internal preprocessing and input matrix will be
#' directly used for rank-sum calculation. Default \code{FALSE} and raw count
#' input is recommended.
#' @param lfcThresh Threshold on log fold-change to identify up-regulated
#' features. Default \code{0.1}.
#' @param returnStats Logical. Whether to return the full statistics table
#' rather then returning the selected genes. Default \code{FALSE}
selectTopFeatures.default <- function(
        x,
        clusterVar,
        vertices,
        nTop = 30,
        processed = FALSE,
        lfcThresh = 0.1,
        returnStats = FALSE,
        ...
) {
    if (ncol(x) != length(clusterVar))
        stop("number of columns of x does not match length of clusterVar")

    if (isFALSE(processed) && !is.rawCounts(x)) {
        warning("Input matrix is not raw counts (integers). ",
                "Results may be affected.")
    }

    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    else clusterVar <- droplevels(clusterVar)

    clusters <- clusterVar[!is.na(clusterVar)]
    x <- x[, !is.na(clusterVar)]

    vcheck <- .checkVertex(x, clusters, vertices, n = NULL)
    clusters <- vcheck[[1]]
    vertices <- vcheck[[2]]

    groupSize <- as.numeric(table(clusters))
    if (sum(groupSize > 0) < 2) {
        stop("Must have at least 2 non-empty groups defined.")
    }

    if (isFALSE(processed)) x <- colNormalize(x, scaleFactor = 1e10, log = TRUE)

    statsTable <- wilcoxauc(x, clusters)

    if (isTRUE(returnStats)) return(statsTable)

    # Take the parts for selected vertices
    statsTable <- statsTable[statsTable$group %in% vertices,]
    # Take only upregulated marker features
    statsTable <- statsTable[statsTable$logFC > lfcThresh,]
    # Take the top N features per group, ordered by logFC
    statsTable <- split(statsTable, droplevels(statsTable$group))
    selected <- unlist(lapply(statsTable, function(tab) {
        tab <- tab[order(tab$padj, -tab$logFC)[seq_len(nTop)],]
        selected <- tab$feature
        message("Selected ", length(selected), " features for \"",
                tab$group[1], "\".")
        return(selected)
    }), use.names = FALSE)
    return(selected[!is.na(selected)])
}

#' @rdname selectTopFeatures
#' @export
#' @method selectTopFeatures Seurat
#' @param layer For "Seurat" method, which layer of the assay to be used.
#' Default \code{"counts"}.
#' @param assay Assay name of the Seurat object to be used. Default \code{NULL}.
#' @examples
#' \donttest{
#' # Seurat example
#' library(Seurat)
#' srt <- CreateSeuratObject(rnaRaw)
#' Idents(srt) <- rnaCluster
#' gene <- selectTopFeatures(srt, vertices = c("OS", "RE"))
#' }
selectTopFeatures.Seurat <- function(
    x,
    clusterVar = NULL,
    vertices,
    assay = NULL,
    layer = "counts",
    processed = FALSE,
    ...
) {
    value <- .getSeuratData(x, assay = assay, layer = layer,
                            clusterVar = clusterVar)
    mat <- value[[1]]
    clusterVar <- value[[2]]

    if (missing(processed)) processed <- layer != "counts"
    selectTopFeatures(mat, clusterVar, vertices, processed = processed, ...)
}

#' @rdname selectTopFeatures
#' @export
#' @method selectTopFeatures SingleCellExperiment
#' @param assay.type Assay name of the SingleCellExperiment object to be used.
#' Default \code{"counts"}.
#' @examples
#' \donttest{
#' # SingleCellExperiment example
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
#' colLabels(sce) <- rnaCluster
#' gene <- selectTopFeatures(sce, vertices = c("OS", "RE"))
#' }
selectTopFeatures.SingleCellExperiment <- function(
    x,
    clusterVar = NULL,
    vertices,
    assay.type = "counts",
    processed = FALSE,
    ...
) {
    value <- .getSCEData(x, assay.type = assay.type, clusterVar = clusterVar)
    mat <- value[[1]]
    clusterVar <- value[[2]]
    if (missing(processed)) processed <- assay.type != "counts"
    selectTopFeatures(mat, clusterVar, vertices, processed = processed, ...)
}

# By default, all group-against-rest tests are conducted all together.
# Features from groups of interests will be selected downstream after this func.
wilcoxauc <- function(x, clusterVar) {
    if (methods::is(x, 'dgeMatrix')) x <- as.matrix(x)
    if (methods::is(x, 'data.frame')) x <- as.matrix(x)
    if (methods::is(x, 'dgTMatrix')) x <- methods::as(x, 'CsparseMatrix')
    if (methods::is(x, 'TsparseMatrix')) x <- methods::as(x, 'CsparseMatrix')
    if (is.null(row.names(x))) {
        rownames(x) <- paste0('Feature', seq(nrow(x)))
    }
    groupSize <- as.numeric(table(clusterVar))

    ## Compute primary statistics
    n1n2 <- groupSize * (ncol(x) - groupSize)
    # rankRes - list(X_ranked, ties), where X_ranked is obs x feature

    rankRes <- colRanking(x)
    ustat <- computeUstat(rankRes$X_ranked, clusterVar, n1n2, groupSize)
    auc <- t(ustat / n1n2)
    pvals <- computePval(ustat, rankRes$ties, ncol(x), n1n2)
    fdr <- apply(pvals, 2, function(p) stats::p.adjust(p, 'BH'))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
    groupSums <- colAggregateSum(x, clusterVar)
    group_nnz <- colNNZAggr(x, clusterVar)
    group_pct <- t(sweep(group_nnz, 1, as.numeric(table(clusterVar)), "/"))

    group_pct_out <- sweep(-group_nnz, 2, colSums(group_nnz), "+")
    group_pct_out <- sweep(group_pct_out, 1,
                           as.numeric(length(clusterVar) - table(clusterVar)),
                           "/")
    group_pct_out <- t(group_pct_out)

    groupMeans <- t(sweep(groupSums, 1, as.numeric(table(clusterVar)), "/"))

    cs <- colSums(groupSums)
    gs <- as.numeric(table(clusterVar))
    lfc <- Reduce(cbind, lapply(seq_along(levels(clusterVar)), function(g) {
        groupMeans[, g] - (cs - groupSums[g, ])/(length(clusterVar) - gs[g])
    }))

    data.frame(
        feature = rep(row.names(x), times = length(levels(clusterVar))),
        group = factor(rep(levels(clusterVar), each = nrow(x)),
                       levels = levels(clusterVar)),
        avgExpr = as.numeric(groupMeans),
        logFC = as.numeric(lfc),
        statistic = as.numeric(t(ustat)),
        auc = as.numeric(auc),
        pval = as.numeric(pvals),
        padj = as.numeric(fdr),
        pct_in = as.numeric(100 * group_pct),
        pct_out = as.numeric(100 * group_pct_out)
    )
}

computeUstat <- function(Xr, cols, n1n2, groupSize) {
    grs <- rowAggregateSum(Xr, cols)

    if (inherits(Xr, 'dgCMatrix')) {
        # With the ranking of only non-zero features, here the tie-ranking of
        # zeros need to be added.
        nnz <- rowNNZAggr_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
        gnz <- groupSize - nnz
        zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
        ustat <- t((t(gnz) * zero.ranks)) + grs - groupSize *
            (groupSize + 1) / 2
    } else {
        ustat <- grs - groupSize * (groupSize + 1) / 2
    }
    return(ustat)
}

computePval <- function(ustat, ties, N, n1n2) {
    z <- ustat - .5 * n1n2
    z <- z - sign(z) * .5
    .x1 <- N ^ 3 - N
    .x2 <- 1 / (12 * (N ^ 2 - N))
    rhs <- unlist(lapply(ties, function(tvals) {
        (.x1 - sum(tvals ^ 3 - tvals)) * .x2
    }))
    usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
    z <- t(z / usigma)
    pvals <- matrix(2 * stats::pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
}


# Utility function to rank columns of matrix
# x - feature by observation matrix.
colRanking <- function(x) {
    if (inherits(x, "dgCMatrix")) {
        x <- Matrix::t(x)
        # This computes the ranking of non-zero values and the ties
        ties <- cpp_rank_matrix_dgc(x@x, x@p, nrow(x), ncol(x))
        return(list(X_ranked = x, ties = ties))
    } else {
        # This directlu computes the feature ranking and the ties.
        rankMatrix_dense(x)
    }
}

# arggregateSum functions returns ngroup x nFeature matrix
colAggregateSum <- function(x, var) {
    if (inherits(x, "dgCMatrix")) {
        colAggregateSum_sparse(x, as.integer(var) - 1, length(unique(var)))
    } else {
        colAggregateSum_dense(x, as.integer(var) - 1, length(unique(var)))
    }
}

rowAggregateSum <- function(x, var) {
    if (inherits(x, "dgCMatrix")) {
        rowAggregateSum_sparse(x, as.integer(var) - 1, length(unique(var)))
    } else {
        rowAggregateSum_dense(x, as.integer(var) - 1, length(unique(var)))
    }
}

colNNZAggr <- function(x, var) {
    if (inherits(x, "dgCMatrix")) {
        colNNZAggr_sparse(x, as.integer(var) - 1, length(unique(var)))
    } else {
        colNNZAggr_dense(x, as.integer(var) - 1, length(unique(var)))
    }
}
