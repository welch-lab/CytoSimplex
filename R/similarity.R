# X - must be feature x cell matrix
# clusterVar - must match its length to ncol(X)
calcSim <- function(
        X,
        clusterVar,
        vertices,
        method = c("euclidean", "cosine", "pearson", "spearman"),
        distKernel = c("gaussian", "log"),
        scale = TRUE,
        force = FALSE,
        sigma = 0.08
) {
    method <- match.arg(method)
    distKernel <- match.arg(distKernel)
    if (nrow(X) > 500 && !isTRUE(force)) {
        stop("Detected more than 500 (", nrow(X),
             ") features in input matrix. Calculation will be slow and ",
             "result will be affected. Selection on features is recommended.",
             "Use `force = TRUE` to continue.")
    }

    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    vertices <- unique(vertices)

    simMat <- matrix(0, nrow = ncol(X), ncol = length(vertices))
    for (i in seq_along(vertices)) {
        centroid <- rowMeans(X[,clusterVar == vertices[i]])
        if (method %in% c("pearson", "spearman")) {
            simMat[,i] <- stats::cor(as.matrix(X), matrix(centroid, ncol = 1), method = method)
        } else if (method == "euclidean") {
            simMat[,i] <- euclideanDist(X, centroid)
        } else if (method == "cosine") {
            simMat[,i] <- cosineDist(X, centroid)
        }
    }
    if (method %in% c("euclidean", "cosine")) {
        simMat <- t(apply(simMat, 1, .normalize))
        simMat <- exp(-simMat^2 / sigma)
    }
    if (isTRUE(scale)) {
        simMat <- apply(simMat, 2, .scaleMinMax)
    }

    simMat <- t(apply(simMat, 1, .normalize))

    simMat <- as.data.frame(simMat)
    rownames(simMat) <- colnames(X)
    colnames(simMat) <- vertices
    attributes(simMat)$class <- c("simMat", class(simMat))
    return(simMat)
}

# calcDist <- function(
#         X,
#         clusterVar,
#         vertices,
#         method = c("euclidean", "cosine", "pearson", "spearman"),
#         normCluster = FALSE,
#         scale = TRUE,
#         force = FALSE
# ){
#     method <- match.arg(method)
#     if (nrow(X) > 500 && !isTRUE(force)) {
#         stop("Detected more than 500 (", nrow(X),
#                 ") features in input matrix. ",
#                 "Calculation will be slow and result will be affected. ",
#                 "Selection on features is recommended.", immediate. = TRUE)
#     }
#     if (length(clusterVar) != ncol(X)) {
#         stop("Length of `clusterVar` has to match `ncol(X)`")
#     }
#     if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
#     vertices <- unique(vertices)
#     if (length(vertices) < 2 || length(vertices) > 4) {
#         stop("Can only use 2 - 4 clusters.")
#     }
#
#     clusterSub <- clusterVar[clusterVar %in% vertices, drop = TRUE]
#     Xsub <- X[, clusterVar %in% vertices, drop = FALSE]
#     if (method %in% c("pearson", "spearman")) {
#         distMat <- stats::cor(as.matrix(X), as.matrix(Xsub),
#                               method = method) + 1
#     } else if (method == "euclidean") {
#         # distMat <- euclidean_dense(as.matrix(X), as.matrix(Xsub))
#         distMat <- euclideanDist(X, Xsub)
#     } else {
#         distMat <- cosineDist(X, Xsub)
#     }
#     # Now distMat ncol(X) x ncol(Xsel)
#     # gc()
#     distMat <- sapply(vertices, function(clust) {
#         rowMeans(distMat[, clusterSub == clust], na.rm = TRUE)
#     })
#     # Now distMat ncol(X) x length(vertices)
#
#     if (method %in% c("euclidean", "cosine")) distMat <- -log10(distMat)
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
#     # distDF$Label <- clusterVar
#     # attributes(distDF)$method <- method
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
