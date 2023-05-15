library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
data("rnaVelo", package = "scPlotSimplex")
vertices <- c("OS", "RE", "CH")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)
rnaLog <- colNormalize(rnaRaw, scaleFactor = 1e4, log = TRUE)

test_that("Test ternary - sparse", {
    expect_error(plotTernary(rnaLog[gene,], rnaCluster, "hi"),
                 "Must specify 3 different vertices.")
    expect_error(plotTernary(rnaLog[gene,], rnaCluster, c("hi", "hey", "yo")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotTernary(rnaLog[gene,], rnaCluster, vertices,
                             dotColor = c("a", "b")),
                 "`dotColor` need to be either 1")
    expect_error(plotTernary(rnaLog[gene,], rnaCluster, vertices,
                             veloGraph = rnaVelo[1:10,]),
                 "`veloGraph must be of shape N x N and has dimnames covering ")
    expect_warning(plotTernary(rnaLog[gene,], rnaCluster, c(vertices, "ORT")),
                   "More than 3 vertices specified for ternary plot. ")


    p <- plotTernary(rnaLog[gene,], rnaCluster, vertices, veloGraph = rnaVelo)
    expect_s3_class(p, "ggplot")

    pl <- plotTernary(rnaLog[gene,], rnaCluster, vertices, splitCluster = TRUE)
    expect_identical(class(pl), "list")

    pl <- plotTernary(rnaLog[gene,], rnaCluster, vertices, splitCluster = TRUE,
                     clusterTitle = FALSE)
    expect_identical(class(pl), "list")
})

test_that("Test ternary - dense", {
    rnaLogSub <- as.matrix(rnaLog[gene,])
    p <- plotTernary(rnaLogSub, rnaCluster, vertices)
    expect_s3_class(p, "ggplot")
})
