library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
data("rnaVelo", package = "scPlotSimplex")
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

    pl <- plotTernary(rnaRaw, rnaCluster, vertices, gene, splitCluster = TRUE)
    expect_identical(class(pl), "list")

    pl <- plotTernary(rnaRaw, rnaCluster, vertices, gene, splitCluster = TRUE,
                     clusterTitle = FALSE)
    expect_identical(class(pl), "list")
})

test_that("Test ternary - dense", {
    rnaRawSub <- as.matrix(rnaRaw[gene,])
    p <- plotTernary(rnaRawSub, rnaCluster, vertices)
    expect_s3_class(p, "ggplot")
})
