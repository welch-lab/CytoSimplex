#' Pick top differentially presented features for similarity calculation
#' @description
#' Performs wilcoxon rank sum test on input matrix. Comparisons are formed by
#' observations in each group, as defined by \code{clusterVar}, against the
#' other observations. The U-Statistics (\code{statistic}), p-value
#' (\code{pval}) and adjusted p-value (\code{padj}), together with average
#' presence in group (\code{avgExpr}), log fold-change (\code{logFC}), AUC
#' (\code{auc}), percentage in group (\code{pct_in}) and percentage out of group
#' (\code{pct_out}) will be calculated.
#'
#' For single-cell gene expression data, usually, log-transformed normalized
#' gene expression value is used for wilcoxon test. While users might have
#' various forms of matrices, fast preprocessing options are provided. Users are
#' recommended to know how the input matrix has been processed, and tweak the
#' options as needed.
#' @param x Dense or sparse matrix, observation per column.
#' @param clusterVar Grouping labels. Length must match with \code{ncol(x)}
#' @param vertices Vector of cluster names that will be used for plotting.
#' @param nTop Number of top differentially presented features per cluster.
#' Default \code{30}.
#' @param normalize,scaleFactor,log Input data preprocessing options. When all
#' specified, they happens in the order of normalization, multiplying scale
#' factor and finally log transformation. Default combination, \code{FALSE},
#' \code{NULL} and \code{FALSE}, respectively, expects a log-transformed
#' normalized input.
#' @param padjThresh Threshold on adjusted p-value to identify significant
#' features. Default \code{0.01}.
#' @param returnStats Logical. Whether to return the whole statistics table
#' rather then returning the selected genes. Default \code{FALSE}
#' @export
#' @return When \code{returnStats = FALSE} (default), a character vector of
#' \code{length(unique(vertices))*nTop} feature names. When \code{returnStats =
#' TRUE}, a data.frame of wilcoxon rank sum test statistics.
selectTopGenes <- function(
        x,
        clusterVar,
        vertices,
        nTop = 30,
        normalize = FALSE,
        scaleFactor = NULL,
        log = FALSE,
        padjThresh = 0.01,
        returnStats = FALSE
) {
    if (methods::is(x, 'dgeMatrix')) x <- as.matrix(x)
    if (methods::is(x, 'data.frame')) x <- as.matrix(x)
    if (methods::is(x, 'dgTMatrix')) x <- methods::as(x, 'CsparseMatrix')
    if (methods::is(x, 'TsparseMatrix')) x <- methods::as(x, 'CsparseMatrix')
    if (ncol(x) != length(clusterVar))
        stop("number of columns of x does not match length of clusterVar")

    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    else clusterVar <- droplevels(clusterVar)

    clusterVar <- clusterVar[!is.na(clusterVar)]
    x <- x[, !is.na(clusterVar)]

    groupSize <- as.numeric(table(clusterVar))
    if (length(groupSize[groupSize > 0]) < 2) {
        stop('Must have at least 2 non-empty groups defined.')
    }

    # Preprocessing as needed
    if (isTRUE(normalize)) x <- colNormalize(x)
    if (!is.null(scaleFactor)) x <- x * scaleFactor
    if (isTRUE(log)) x <- log1p(x)

    statsTable <- wilcoxauc(x, clusterVar)

    if (isTRUE(returnStats)) return(statsTable)

    # Take the parts for selected vertices
    statsTable <- statsTable[statsTable$group %in% vertices,]
    # Take only significant features
    statsTable <- statsTable[statsTable$padj < padjThresh,]
    # Remove features identified as marker for multiple groups but has less
    # logFC in some of groups
    statsTable <- statsTable[order(statsTable$logFC, decreasing = TRUE),]
    statsTable <- statsTable[!duplicated(statsTable$feature),]
    # Take the top N features per group, ordered by logFC
    statsTable <- split(statsTable, droplevels(statsTable$group))
    unlist(lapply(statsTable, function(tab) {
        tab <- tab[order(tab$logFC, decreasing = TRUE)[seq_len(nTop)],]
        tab$feature
    }), use.names = FALSE)
}

# By default, all group-against-rest tests are conducted all together.
# Features from groups of interests will be selected downstream.
wilcoxauc <- function(x, clusterVar) {
    if (is.null(row.names(x))) {
        rownames(x) <- paste0('Feature', seq(nrow(x)))
    }
    groupSize <- as.numeric(table(clusterVar))

    ## Compute primary statistics
    n1n2 <- groupSize * (ncol(x) - groupSize)
    # rankRes - list(X_ranked, ties), where X_ranked is obs x feature
    rankRes <- rankMatrix(x)
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

    group_means <- t(sweep(groupSums, 1, as.numeric(table(clusterVar)), "/"))

    cs <- colSums(groupSums)
    gs <- as.numeric(table(clusterVar))
    lfc <- Reduce(cbind, lapply(seq_along(levels(clusterVar)), function(g) {
        group_means[, g] -
            (cs - groupSums[g,]) / (length(clusterVar) - gs[g])
    }))

    data.frame(
        feature = rep(row.names(x), times = length(levels(clusterVar))),
        group = factor(rep(levels(clusterVar), each = nrow(x)),
                       levels = levels(clusterVar)),
        avgExpr = as.numeric(group_means),
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
        gnz <- (groupSize - rowNNZAggr(Xr, cols))
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
rankMatrix <- function(x) {
    if (inherits(x, "dgCMatrix")) {
        x <- Matrix::t(x)
        ties <- cpp_rank_matrix_dgc(x@x, x@p, nrow(x), ncol(x))
        return(list(X_ranked = x, ties = ties))
    } else {
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

rowNNZAggr <- function(x, var) {
    if (inherits(x, "dgCMatrix")) {
        rowNNZAggr_sparse(x, as.integer(var) - 1, length(unique(var)))
    } else {
        rowNNZAggr_dense(x, as.integer(var) - 1, length(unique(var)))
    }
}
