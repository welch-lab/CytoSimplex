# X - must be feature x cell matrix
# clusterVar - must match its length to ncol(X)
calcDist <- function(
        X,
        clusterVar,
        vertices,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        normCluster = FALSE,
        scale = TRUE
){
    method <- match.arg(method)
    if (nrow(X) > 500) {
        warning("Detected more than 500 (", nrow(X),
                ") features in input matrix. ",
                "Calculation will be slow and result will be affected. ",
                "Selection on features is recommended.", immediate. = TRUE)
    }
    if (length(clusterVar) != ncol(X)) {
        stop("Length of `clusterVar` has to match `ncol(X)`")
    }
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    vertices <- unique(vertices)
    if (length(vertices) < 2 || length(vertices) > 4) {
        stop("Can only use 2 - 4 clusters.")
    }

    clusterSub <- clusterVar[clusterVar %in% vertices, drop = TRUE]
    Xsub <- X[, clusterVar %in% vertices, drop = FALSE]
    if (method %in% c("pearson", "spearman")) {
        distMat <- stats::cor(as.matrix(X), as.matrix(Xsub),
                              method = method) + 1
    } else if (method == "euclidean") {
        # distMat <- euclidean_dense(as.matrix(X), as.matrix(Xsub))
        distMat <- euclideanDist(X, Xsub)
    } else {
        distMat <- cosineDist(X, Xsub)
    }
    # Now distMat ncol(X) x ncol(Xsel)
    # gc()
    distMat <- sapply(vertices, function(clust) {
        rowMeans(distMat[, clusterSub == clust], na.rm = TRUE)
    })
    # Now distMat ncol(X) x length(vertices)

    if (method %in% c("euclidean", "cosine")) distMat <- -log10(distMat)

    if (isTRUE(normCluster)) {
        distMat <- apply(distMat, 2, .normalize)
    }

    distMat <- t(apply(distMat, 1, .normalize))

    if (isTRUE(scale)) {
        distMat <- apply(distMat, 2, .scaleMinMax)
    }

    distMat <- t(apply(distMat, 1, .normalize))

    distDF <- as.data.frame(distMat)
    rownames(distDF) <- colnames(X)
    colnames(distDF) <- vertices
    distDF$Label <- clusterVar
    attributes(distDF)$method <- method
    attributes(distDF)$class <- c("distMatrix", class(distDF))
    return(distDF)
}

cosineDist <- function(query, target) {
    if (inherits(query, "matrix")) cosine_dense(query, target)
    else if (inherits(query, "dgCMatrix")) cosine_sparse(query, target) # TODO not implemented yet
    else {
        stop("Distance calculation not supported for matrix of class ",
             class(query)[1])
    }
}

euclideanDist <- function(query, target) {
    if (inherits(query, "matrix")) euclidean_dense(query, target)
    else if (inherits(query, "dgCMatrix")) euclidean_sparse(query, target) # TODO not implemented yet
    else {
        stop("Distance calculation not supported for matrix of class ",
             class(query)[1])
    }
}
