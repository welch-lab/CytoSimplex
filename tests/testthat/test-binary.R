library(testthat)
library(CytoSimplex)
library(Matrix)

data("rnaRaw", package = "CytoSimplex")
data("rnaCluster", package = "CytoSimplex")
vertices <- c("OS", "RE")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)

test_that("Test binary - sparse", {
    expect_error(plotBinary(rnaRaw, rnaCluster[1:100], c("CH", "ORT")),
                 "Length of `clusterVar` must match")
    expect_error(plotBinary(rnaRaw, rnaCluster, "hi"),
                 "Must specify 2 different vertices.")
    expect_error(plotBinary(rnaRaw, rnaCluster, c("hi", "hey")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotBinary(rnaRaw, rnaCluster, vertices, dotColor = c("a", "b")),
                 "`dotColor` need to be either 1")
    expect_error(plotBinary(rnaRaw, rnaCluster, vertices),
                 "Detected more than 500")

    expect_warning(plotBinary(rnaRaw, rnaCluster, c(vertices, "CH"), features = gene),
                   "2 vertices are expected while 3 are specified")
    rnaNorm <- colNormalize(rnaRaw)
    expect_warning(plotBinary(rnaNorm, rnaCluster, vertices, gene),
                   "Input matrix is not raw counts")
    rnaCluster.char <- as.character(rnaCluster)
    p <- plotBinary(rnaRaw, rnaCluster.char, vertices, gene)
    expect_s3_class(p, "ggplot")

    expect_no_error(
        plotBinary(rnaRaw, rnaCluster, vertices, gene, method = "cosine")
    )
    expect_no_error(
        plotBinary(rnaRaw, rnaCluster, vertices, gene, method = "pearson")
    )
    expect_no_error(
        plotBinary(rnaRaw, rnaCluster, vertices, gene, method = "spearman")
    )

    pl <- plotBinary(rnaRaw, rnaCluster, vertices, gene, byCluster = "all")
    expect_identical(class(pl), "list")
    expect_s3_class(pl[[1]], "ggplot")

    pl <- plotBinary(rnaRaw, rnaCluster, vertices, gene, byCluster = c("RE"))
    expect_identical(class(pl), "list")
    expect_s3_class(pl[[1]], "ggplot")

    expect_error(plotBinary(rnaRaw, rnaCluster, vertices, gene,
                            byCluster = "Hi"),
                 "`byCluster` must be either a vector of cluster name ")

    simData <- plotBinary(rnaRaw, rnaCluster, vertices, gene, returnData = TRUE)
    expect_identical(class(simData), "list")
    expect_s3_class(simData[[1]], "simMat")
})

test_that("Test binary - dense", {
    rnaRawSub <- as.matrix(rnaRaw[gene,])
    p <- plotBinary(rnaRawSub, rnaCluster, vertices)
    expect_s3_class(p, "ggplot")

    expect_no_error(
        plotBinary(rnaRawSub, rnaCluster, vertices, method = "cosine")
    )
})
