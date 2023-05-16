# X - must be feature x cell matrix
# clusterVar - must match its length to ncol(X)
calcDist2 <- function(
        X,
        clusterVar,
        vertices,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        scale = TRUE,
        force = FALSE,
        sigma = 100
) {
    method <- match.arg(method)
    if (nrow(X) > 500 && !isTRUE(force)) {
        stop("Detected more than 500 (", nrow(X),
             ") features in input matrix. Calculation will be slow and ",
             "result will be affected. Selection on features is recommended.",
             "Use `force = TRUE` to continue.")
    }

    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    vertices <- unique(vertices)

    distMat <- matrix(0, nrow = ncol(X), ncol = length(vertices))
    for (i in seq_along(vertices)) {
        centroid <- rowMeans(X[,clusterVar == vertices[i]])
        if (method %in% c("pearson", "spearman")) {
            distMat[,i] <- stats::cor(as.matrix(X), matrix(centroid, ncol = 1))
        } else if (method == "euclidean") {
            distMat[,i] <- euclideanDist(X, centroid)
        } else if (method == "cosine") {
            distMat[,i] <- cosineDist(X, centroid)
        }
    }
    if (method %in% c("euclidean", "cosine")) {
        # Apply Gaussian Kernel to convert distance to similarity
        distMat <- exp(-distMat^2 / sigma)
    }

    if (isTRUE(scale)) {
        distMat <- apply(distMat, 2, .scaleMinMax)
    }

    distMat <- t(apply(distMat, 1, .normalize))

    distDF <- as.data.frame(distMat)
    rownames(distDF) <- colnames(X)
    colnames(distDF) <- vertices
    attributes(distDF)$class <- c("simMat", class(distDF))
    return(distDF)
}

# calcDist <- function(
#         X,
#         clusterVar,
#         vertices,
#         method = c("euclidean", "cosine", "pearson", "spearman"),
#         normCluster = FALSE,
#         scale = TRUE,
#         force = FALSE
# ) {
#     method <- match.arg(method)
#     if (nrow(X) > 500 && !isTRUE(force)) {
#         stop("Detected more than 500 (", nrow(X),
#              ") features in input matrix. Calculation will be slow and ",
#              "result will be affected. Selection on features is recommended.",
#              "Use `force = TRUE` to continue.")
#     }
#
#     if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
#     vertices <- unique(vertices)
#
#     if (length(vertices) < 2 || length(vertices) > 4) {
#         stop("Can only use 2 - 4 clusters.")
#     }
#     if (!all(vertices %in% clusterVar)) {
#         stop("Specified vertices not found in `clusterVar`.")
#     }
#
#     # Normalize raw count input by column (cell) library size
#     clusterSub <- clusterVar[clusterVar %in% vertices, drop = TRUE]
#     Xsub <- X[, clusterVar %in% vertices, drop = FALSE]
#     if (method %in% c("pearson", "spearman")) {
#         distMat <- stats::cor(as.matrix(X), as.matrix(Xsub),
#                               method = method) + 1
#     } else if (method == "euclidean") {
#         distMat <- euclideanDist(X, Xsub)
#     } else {
#         distMat <- cosineDist(X, Xsub)
#     }
#     # Now distMat ncol(X) x ncol(Xsel)
#     distMat <- sapply(vertices, function(clust) {
#         rowMeans(distMat[, clusterSub == clust], na.rm = TRUE)
#     })
#     # Now distMat ncol(X) x length(vertices)
#     if (method %in% c("euclidean", "cosine")) {
#         distMat <- -log10(distMat)
#     }
#
#     if (isTRUE(normCluster)) {
#         distMat <- apply(distMat, 2, .normalize)
#     }
#
#     distMat <- t(apply(distMat, 1, .normalize))
#
#     if (isTRUE(scale)) {
#         distMat <- apply(distMat, 2, .scaleMinMax)
#     }
#
#     distMat <- t(apply(distMat, 1, .normalize))
#
#     distDF <- as.data.frame(distMat)
#     rownames(distDF) <- colnames(X)
#     colnames(distDF) <- vertices
#     attributes(distDF)$class <- c("simMat", class(distDF))
#     return(distDF)
# }

cosineDist <- function(query, target) {
    if (inherits(query, "dgCMatrix")) {
        if (is.null(dim(target))) {
            target <- Matrix::Matrix(target, ncol = 1, nrow = length(target),
                                     sparse = TRUE)
        }
        cosine_sparse(query, target)
    } else {
        if (is.null(dim(target))) target <- matrix(target, ncol = 1)
        cosine_dense(query, target)
    }
}

euclideanDist <- function(query, target) {
    if (inherits(query, "dgCMatrix")) {
        if (is.null(dim(target))) {
            target <- Matrix::Matrix(target, ncol = 1, nrow = length(target),
                                     sparse = TRUE)
        }
        euclidean_sparse(query, target)
    } else {
        if (is.null(dim(target))) target <- matrix(target, ncol = 1)
        euclidean_dense(query, target)
    }
}
