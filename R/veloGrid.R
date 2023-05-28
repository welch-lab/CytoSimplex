aggrVeloGraph <- function(
        graph,
        clusterVar,
        vertices
) {
    if (!inherits(graph, "dgCMatrix")) {
        graph <- methods::as(graph, "CsparseMatrix")
    }
    veloMat <- sapply(vertices, function(clust) {
        rowMeans(graph[, clusterVar == clust], na.rm = TRUE)
    })

    veloMat = t(apply(veloMat, 1, .normalize))

    rownames(veloMat) <- rownames(graph)
    colnames(veloMat) <- vertices
    return(veloMat)
}

# simMat --- Output from `calcDist()`, ncell x 4 simMat, data.frame object
#             First three columns are normalized similarity, 4th col cluster
# veloMat --- Output from `aggrVeloGraph()`, ncell x 3 matrix object
# Returns a data.frame with 12 columns, where
#   1:3 arrow start ternary coords
#   4:6 arrow end ternary coords, pointing to left
#   7:9 arrow end ternary coords, pointing to top
#   10:12 arrow end ternary coords, pointing to right
calcGridVelo <- function(
        simMat,
        veloMat,
        nGrid = 10,
        radius = 0.1
) {
    # Determine what simplex we are working with
    n <- ncol(simMat)
    # Get the vertex coordinates in cartesian
    shape <- if (n == 3) triangle else if (n == 4) tetra
    # Draw breaks on each axis and use combination to get the full grid
    range <- as.data.frame(rbind(apply(shape, 2, min), apply(shape, 2, max)))
    ## Determin grid edge length just by x axis
    win <- abs(range[1,1] - range[2,1])/nGrid
    seg <- win/2
    gridCart.axis <- lapply(range, function(v) seq(v[1] + seg, v[2] - seg, win))
    gridCart <- expand.grid(gridCart.axis)

    # Select those that fall into the simplex space with barycentric coord
    gridBary <- cart2bary(shape, gridCart)
    gridSelect <- rowSums(gridBary > 0 & gridBary < 1) == ncol(gridBary)
    gridCart <- gridCart[gridSelect,]
    rownames(gridCart) <- paste0("grid", seq(nrow(gridCart)))

    # Convert simplex barycentric coord of cells to cartesien coord (2D space)
    cellCart <- as.data.frame(as.matrix(simMat) %*% shape)

    # Aggregate velocity by grid, comparing cell cartesian coordinate
    gridVelo <- matrix(0, nrow = nrow(gridCart), ncol = ncol(veloMat),
                       dimnames = list(rownames(gridCart),
                                       colnames(veloMat)))
    for (i in seq(nrow(gridCart))) {
        bcIdx <- rep(TRUE, nrow(cellCart))
        for (j in seq(ncol(cellCart))) {
            bcIdx <- bcIdx &
                cellCart[,j] > (gridCart[i,j] - seg) &
                cellCart[,j] < (gridCart[i,j] + seg)
        }
        if (sum(bcIdx, na.rm = TRUE) > 4) {
            # Get the velocity value presented as arrow length only when more
            # then 5 cells fall into a grid
            subVelo <- veloMat[bcIdx, , drop = FALSE]
            gridVelo[i,] <- colMeans(subVelo, na.rm = TRUE)
        }
    }
    # Remove zero-velo grid
    gridToKeep <- rowSums(gridVelo, na.rm = TRUE) > 0
    gridCart <- gridCart[gridToKeep, , drop = FALSE]
    gridVelo <- gridVelo[gridToKeep, , drop = FALSE]

    # Calculate cartesian coordinates of each arrow ending
    arrow.cart <- lapply(seq_len(ncol(gridVelo)), function(i) {
        ends <- getArrowEndCoord(G = gridCart, target = shape[i,],
                                  len = gridVelo[,i] * radius)
        colnames(ends) <- paste0(colnames(ends), "end")
        arrow.cart.sub <- cbind(gridCart, ends)
        return(arrow.cart.sub)
    })
    return(arrow.cart)
}

# G      - coordinate of grid centroid, N x m matrix, G[,1] 'x' of each
#          centroid, G[,2] 'y' of each centroid, etc
# target - length m vector, cartesian coordinate of the vertex where the arrow
#          is pointing to
# len    - arrow length. N element vector
getArrowEndCoord <- function(G, target, len) {
    if (nrow(G) == 0) return(G)
    V <- matrix(rep(target, nrow(G)), nrow = nrow(G), byrow = TRUE)
    # V <- matrix(c(rep(xv, nrow(G)), rep(yv, nrow(G))), ncol = 2)
    # Directed vector pointing from G to V
    vecGV <- V - G
    # Euclidean distance from G to V
    lenGV <- sqrt(rowSums(vecGV * vecGV))
    # Limit arrow length to be no longer than lenGV
    len[len > lenGV] <- lenGV[len > lenGV]
    # Directed vector pointing from G to A
    vecGA <- (vecGV / lenGV) * len
    A <- G + vecGA
    # colnames(A) <- c("xend", "yend")
    return(A)
}

# Adopted from geometry::cart2bary
cart2bary <- function(X, P) {
    X <- as.matrix(X)
    P <- as.matrix(P)
    M <- nrow(P)
    N <- ncol(P)
    X1 <- X[1:N, ] - (matrix(1, N, 1) %*% X[N + 1, , drop = FALSE])
    Beta <- (P - matrix(X[N + 1, ], M, N, byrow = TRUE)) %*%
        solve(X1)
    Beta <- cbind(Beta, 1 - apply(Beta, 1, sum))
    return(Beta)
}
