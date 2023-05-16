library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
vertices <- c("OS", "RE")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)
rnaLog <- colNormalize(rnaRaw, scaleFactor = 1e4, log = TRUE)

test_that("Test binary - sparse", {
    expect_error(plotBinary(rnaLog[gene,], rnaCluster, "hi"),
                 "Must specify 2 different vertices.")
    expect_error(plotBinary(rnaLog[gene,], rnaCluster, c("hi", "hey")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotBinary(rnaLog[gene,], rnaCluster, vertices,
                            dotColor = c("a", "b")),
                 "`dotColor` need to be either 1")
    expect_error(plotBinary(rnaLog, rnaCluster, vertices),
                 "Detected more than 500")

    expect_warning(plotBinary(rnaLog[gene,], rnaCluster, c(vertices, "CH")),
                   "2 vertices are expected while 3 are specified")


    p <- plotBinary(rnaLog[gene,], rnaCluster, vertices)
    expect_s3_class(p, "ggplot")

    expect_no_error(
        plotBinary(rnaLog[gene,], rnaCluster, vertices, method = "cosine")
    )
    expect_no_error(
        plotBinary(rnaLog[gene,], rnaCluster, vertices, method = "pearson")
    )
    expect_no_error(
        plotBinary(rnaLog[gene,], rnaCluster, vertices, method = "spearman")
    )

    pl <- plotBinary(rnaLog[gene,], rnaCluster, vertices, splitCluster = TRUE)
    expect_identical(class(pl), "list")

    pl <- plotBinary(rnaLog[gene,], rnaCluster, vertices, splitCluster = TRUE,
                     clusterTitle = FALSE)
    expect_identical(class(pl), "list")
})

test_that("Test binary - dense", {
    rnaLogSub <- as.matrix(rnaLog[gene,])
    p <- plotBinary(rnaLogSub, rnaCluster, vertices)
    expect_s3_class(p, "ggplot")

    expect_no_error(
        plotBinary(rnaLogSub, rnaCluster, vertices, method = "cosine")
    )
    expect_no_error(
        plotBinary(rnaLogSub, rnaCluster, vertices, method = "pearson")
    )
    expect_no_error(
        plotBinary(rnaLogSub, rnaCluster, vertices, method = "spearman")
    )
})
