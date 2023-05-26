library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
data("rnaVelo", package = "scPlotSimplex")
vertices <- c("OS", "RE", "CH", "ORT")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)
rnaLog <- colNormalize(rnaRaw, scaleFactor = 1e4, log = TRUE)

test_that("Test quaternary - sparse", {
    expect_error(plotQuaternary(rnaLog[gene,], rnaCluster, "hi"),
                 "Must specify 4 different vertices.")
    expect_error(plotQuaternary(rnaLog[gene,], rnaCluster, c("hi", "hey", "yo", "nah")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotQuaternary(rnaLog[gene,], rnaCluster, vertices,
                             dotColor = c("a", "b")),
                 "`dotColor` need to be either 1")
    expect_error(plotQuaternary(rnaLog[gene,], rnaCluster, vertices,
                             veloGraph = rnaVelo[1:10,]),
                 "`veloGraph must be of shape N x N and has dimnames covering ")
    expect_warning(plotQuaternary(rnaLog[gene,], rnaCluster, c(vertices, "Stem")),
                   "4 vertices are expected while 5 are specified.")


    p <- plotQuaternary(rnaLog[gene,], rnaCluster, vertices, veloGraph = rnaVelo)
    expect_s3_class(p, "plist")

    pl <- plotQuaternary(rnaLog[gene,], rnaCluster, vertices, splitCluster = TRUE)
    expect_identical(class(pl), "list")

    pl <- plotQuaternary(rnaLog[gene,], rnaCluster, vertices, splitCluster = TRUE,
                      clusterTitle = FALSE)
    expect_identical(class(pl), "list")

    show(p)
    expect_gt(length(dev.list()), 0)
})

test_that("Test quaternary - dense", {
    rnaLogSub <- as.matrix(rnaLog[gene,])
    p <- plotQuaternary(rnaLogSub, rnaCluster, vertices)
    expect_s3_class(p, "plist")
})

test_that("Test quaternary GIF", {
    grouping <- list(A = c("ORT"),
                     B = c("RE", "OS"),
                     C = "CH",
                     D = "Stem")
    expect_error(writeQuaternaryGIF(rnaLog[gene,], clusterVar = rnaCluster,
                                    vertices = grouping, fps = 33),
                 "FPS must be a factor of 100.")
    writeQuaternaryGIF(rnaLog[gene,], clusterVar = rnaCluster,
                       vertices = grouping,
                       gifPath = "test.gif", tmpDir = "testGif/")
    expect_true(dir.exists("testGif"))
    expect_true(file.exists("test.gif"))

    writeQuaternaryGIF(rnaLog[gene,], clusterVar = rnaCluster,
                       vertices = grouping, useCluster = "Stem",
                       gifPath = "test2.gif", tmpDir = "testGif2/")
    expect_true(dir.exists("testGif2"))
    expect_true(file.exists("test2.gif"))
})
