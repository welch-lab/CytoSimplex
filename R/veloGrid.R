aggrVeloGraph <- function(
        graph,
        clusterVar,
        vertices
) {
    if (!inherits(graph, "dgCMatrix")) {
        graph <- methods::as(graph, "CsparseMatrix")
    }
    veloMat <- sapply(vertices, function(clust) {
        rowMeans(graph[, clusterVar == clust])
    })

    veloMat = t(apply(veloMat, 1, .normalize))

    rownames(veloMat) <- rownames(graph)
    colnames(veloMat) <- vertices
    return(veloMat)
}

# distMat --- Output from `calcDist()`, ncell x 4 distMatrix, data.frame object
#             First three columns are normalized similarity, 4th col cluster
# veloMat --- Output from `aggrVeloGraph()`, ncell x 3 matrix object
# Returns a data.frame with 12 columns, where
#   1:3 arrow start ternary coords
#   4:6 arrow end ternary coords, pointing to left
#   7:9 arrow end ternary coords, pointing to top
#   10:12 arrow end ternary coords, pointing to right
calcGridVelo <- function(
        distMat,
        veloMat,
        nGrid = 10,
        radius = 0.1
) {
    win <- 1/nGrid
    seg <- win/2

    # Initialize `nGrid`^2 grid centroids coords in a 1x1 square
    gridCentroidCoord <- data.frame(
        x = rep(seq(win/2, 1 - win/2, win), each = nGrid),
        y = rep(seq(win/2, 1 - win/2, win), nGrid)
    )

    # Select those that fall into the equilateral triangle by y/x < tan(pi/3)
    inLeft <- gridCentroidCoord$y/gridCentroidCoord$x < 3^0.5
    inRight <- gridCentroidCoord$y/(1 - gridCentroidCoord$x) < 3^0.5
    gridCentroidCoord <- gridCentroidCoord[inLeft & inRight, ]
    rownames(gridCentroidCoord) <- paste0("grid", seq(nrow(gridCentroidCoord)))

    # Convert the cell similarity matrix from ternary into 2D space
    distMat <- distMat[,seq(3)]
    colnames(distMat) <- c("x", "y", "z")
    cellCoord <- as.data.frame(as.matrix(distMat) %*% triangle)

    # Aggregate velocity by grid, comparing cell 2D coordinate
    gridVelo <- matrix(0, nrow = nrow(gridCentroidCoord), ncol = 3,
                       dimnames = list(rownames(gridCentroidCoord),
                                       colnames(veloMat)))
    for (i in seq(nrow(gridCentroidCoord))) {
        bcIdx <- cellCoord$x > (gridCentroidCoord$x[i] - seg) &
            cellCoord$x < (gridCentroidCoord$x[i] + seg) &
            cellCoord$y > (gridCentroidCoord$y[i] - seg) &
            cellCoord$y < (gridCentroidCoord$y[i] + seg)

        if (sum(bcIdx) > 0) {
            subVelo <- veloMat[bcIdx, , drop = FALSE]
            gridVelo[i,] <- colMeans(subVelo + 1e-8,
                                     na.rm = TRUE)
        }
    }
    # Remove zero-velo grid
    gridToKeep <- rowSums(gridVelo, na.rm = TRUE) > 0
    gridCentroidCoord <- gridCentroidCoord[gridToKeep, , drop = FALSE]
    gridVelo <- gridVelo[gridToKeep, , drop = FALSE]

    # message(nrow(veloMat), " cells and ", nrow(gridVelo), " grids")
    # Calculate arrow ending 2D coords
    leftEnds <- getArrowEndCoord(G = gridCentroidCoord, xv = 0, yv = 0,
                                 len = gridVelo[,1] * radius)
    topEnds <- getArrowEndCoord(G = gridCentroidCoord, xv = 0.5, yv = 3^0.5/2,
                                 len = gridVelo[,2] * radius)
    rightEnds <- getArrowEndCoord(G = gridCentroidCoord, xv = 1, yv = 0,
                                 len = gridVelo[,3] * radius)

    return(list(
        left = cbind(gridCentroidCoord, leftEnds),
        top = cbind(gridCentroidCoord, topEnds),
        right = cbind(gridCentroidCoord, rightEnds)
    ))
    # return(list(
    #     grid = tlr2xy(gridCentroidCoord, coord_tern(), inverse = TRUE),
    #     left = tlr2xy(leftEnds, coord_tern(), inverse = TRUE),
    #     top = tlr2xy(topEnds, coord_tern(), inverse = TRUE),
    #     right = tlr2xy(rightEnds, coord_tern(), inverse = TRUE)
    # ))
}

# G - coordinate of grid centroid, N x 2 matrix, G[,1] x of each centroid,
#     G[,2] y of each centroid
# xv, yv - coordinate of one vertex
# len - arrow length. N element vector
#
# Using "G" to denote the grid centroid point, "V" to denote the vertex point
# "A" to denote the arrow end point.
#
# Need to be aware that none of G should overlap with V
# According to how grids are initialized, we are fine though.
getArrowEndCoord <- function(G, xv, yv, len) {
    if (nrow(G) == 0) return(G)
    V <- matrix(c(rep(xv, nrow(G)), rep(yv, nrow(G))), ncol = 2)
    # Directed vector pointing from G to V
    vecGV <- V - G
    # Euclidean distance from G to V
    lenGV <- sqrt(rowSums(vecGV * vecGV))
    # Limit arrow length to be no longer than lenGV
    len[len > lenGV] <- lenGV[len > lenGV]
    # Directed vector pointing from G to A
    vecGA <- (vecGV / lenGV) * len
    A <- G + vecGA
    colnames(A) <- c("xend", "yend")
    return(A)
}
