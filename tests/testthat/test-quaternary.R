library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
data("rnaVelo", package = "scPlotSimplex")
vertices <- c("OS", "RE", "CH", "ORT")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)

test_that("Test quaternary - sparse", {
    expect_error(plotQuaternary(rnaRaw, rnaCluster, "hi"),
                 "Must specify 4 different vertices.")
    expect_error(plotQuaternary(rnaRaw, rnaCluster, c("hi", "hey", "yo", "nah")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotQuaternary(rnaRaw, rnaCluster, vertices,
                             dotColor = c("a", "b")),
                 "`dotColor` need to be either 1")
    expect_error(plotQuaternary(rnaRaw, rnaCluster, vertices, gene,
                             veloGraph = rnaVelo[1:10,]),
                 "`veloGraph must be of shape N x N and has dimnames covering ")
    expect_warning(plotQuaternary(rnaRaw, rnaCluster, c(vertices, "Stem"), gene),
                   "4 vertices are expected while 5 are specified.")
    rnaNorm <- colNormalize(rnaRaw)
    expect_warning(plotQuaternary(rnaNorm, rnaCluster, vertices, gene),
                   "Input matrix is not raw counts")

    p <- plotQuaternary(rnaRaw, rnaCluster, vertices, gene, veloGraph = rnaVelo)
    expect_s3_class(p, "plist")

    pl <- plotQuaternary(rnaRaw, rnaCluster, vertices, gene, splitCluster = TRUE)
    expect_identical(class(pl), "list")

    pl <- plotQuaternary(rnaRaw, rnaCluster, vertices, gene, splitCluster = TRUE,
                      clusterTitle = FALSE)
    expect_identical(class(pl), "list")

    show(p)
    expect_gt(length(dev.list()), 0)
})

test_that("Test quaternary - dense", {
    rnaRawSub <- as.matrix(rnaRaw[gene,])
    p <- plotQuaternary(rnaRawSub, rnaCluster, vertices)
    expect_s3_class(p, "plist")
})

test_that("Test quaternary GIF", {
    grouping <- list(A = c("ORT"),
                     B = c("RE", "OS"),
                     C = "CH",
                     D = "Stem")
    expect_error(writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                                    features = gene, vertices = grouping,
                                    fps = 33),
                 "FPS must be a factor of 100.")
    writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                       vertices = grouping, features = gene,
                       gifPath = "test.gif", tmpDir = "testGif/")
    expect_true(dir.exists("testGif"))
    expect_true(file.exists("test.gif"))
    unlink("testGif/", recursive = TRUE)
    unlink("test.gif")

    writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster, features = gene,
                       vertices = grouping, useCluster = "Stem",
                       gifPath = "test2.gif", tmpDir = "testGif2/")
    expect_true(dir.exists("testGif2"))
    expect_true(file.exists("test2.gif"))
    unlink("testGif2/", recursive = TRUE)
    unlink("test2.gif")
})
