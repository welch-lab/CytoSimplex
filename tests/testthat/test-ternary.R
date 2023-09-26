library(testthat)
library(CytoSimplex)
library(Matrix)

data("rnaRaw", package = "CytoSimplex")
data("rnaCluster", package = "CytoSimplex")
data("rnaVelo", package = "CytoSimplex")
vertices <- c("OS", "RE", "CH")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)

test_that("Test ternary - sparse", {
    expect_error(plotTernary(rnaRaw, rnaCluster, "hi"),
                 "Must specify 3 different vertices.")
    expect_error(plotTernary(rnaRaw, rnaCluster, c("hi", "hey", "yo")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotTernary(rnaRaw, rnaCluster, vertices,
                             dotColor = c("a", "b")),
                 "`dotColor` need to be either 1")
    expect_error(plotTernary(rnaRaw, rnaCluster, vertices, gene,
                             veloGraph = rnaVelo[1:10,]),
                 "`veloGraph must be of shape N x N and has dimnames covering ")
    expect_warning(plotTernary(rnaRaw, rnaCluster, c(vertices, "ORT"), gene),
                   "3 vertices are expected while 4 are specified")
    rnaNorm <- colNormalize(rnaRaw)
    expect_warning(plotTernary(rnaNorm, rnaCluster, vertices, gene),
                   "Input matrix is not raw counts")

    p <- plotTernary(rnaRaw, rnaCluster, vertices, gene, veloGraph = rnaVelo)
    expect_s3_class(p, "ggplot")

    pl <- plotTernary(rnaRaw, rnaCluster, vertices, gene, byCluster = "all")
    expect_identical(class(pl), "list")

    pl <- plotTernary(rnaRaw, rnaCluster, vertices, gene, byCluster = "RE")
    expect_identical(class(pl), "list")

    expect_error(plotTernary(rnaRaw, rnaCluster, vertices, gene,
                            byCluster = "Hi"),
                 "`byCluster` must be either a vector of cluster name ")

    plotData <- plotTernary(rnaRaw, rnaCluster, vertices, gene,
                            veloGraph = rnaVelo, returnData = TRUE)
    expect_identical(class(plotData), "list")
    expect_s3_class(plotData[[1]], "simMat")
    expect_identical(class(plotData[[2]]), c("matrix", "array"))
})

test_that("Test ternary - dense", {
    rnaRawSub <- as.matrix(rnaRaw[gene,])
    p <- plotTernary(rnaRawSub, rnaCluster, vertices)
    expect_s3_class(p, "ggplot")
})
