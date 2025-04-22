library(testthat)
library(Matrix)

Sys.setenv("OMP_THREAD_LIMIT" = 2)
data("rnaRaw", package = "CytoSimplex")
data("rnaCluster", package = "CytoSimplex")
data("rnaVelo", package = "CytoSimplex")
vertices <- c("OS", "RE", "CH", "ORT")
gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)

test_that("Test quaternary - sparse", {
    expect_error(plotQuaternary(rnaRaw, rnaCluster, "hi"),
                 "Must specify 4 different vertices.")
    expect_error(plotQuaternary(rnaRaw, rnaCluster, c("hi", "hey", "yo", "nah")),
                 "Specified vertex clusters are not all found in the cluster ")
    expect_error(plotQuaternary(rnaRaw, rnaCluster, vertices,
                             dotColor = c("a", "b")),
                 "Length of `dotColor` must be 1 for all dots")
    expect_error(plotQuaternary(rnaRaw, rnaCluster, vertices, gene,
                             veloGraph = rnaVelo[1:10,]),
                 "`veloGraph` must have dimension of 250 x 250 and has dimnames covering all")
    expect_message(plotQuaternary(rnaRaw, rnaCluster, c(vertices, "Stem"), gene),
                   "4 vertices are expected while 5 are specified.")
    rnaNorm <- colNormalize(rnaRaw)
    expect_message(plotQuaternary(rnaNorm, rnaCluster, vertices, gene),
                   "Input matrix is not raw counts")

    p <- plotQuaternary(rnaRaw, rnaCluster, vertices, gene, veloGraph = rnaVelo)
    expect_s3_class(p, "plotly")

    pl <- plotQuaternary(rnaRaw, rnaCluster, vertices, gene, byCluster = "all", interactive = FALSE)
    expect_identical(class(pl), "list")

    pl <- plotQuaternary(rnaRaw, rnaCluster, vertices, gene, byCluster = "RE",
                      clusterTitle = FALSE, interactive = FALSE)
    expect_identical(class(pl), "list")

    expect_no_error(show(p))

    expect_error(plotQuaternary(rnaRaw, rnaCluster, vertices, gene,
                                byCluster = "Hi", interactive = FALSE),
                 "`byCluster` must be either a vector of cluster names or just")
    expect_no_error(plotQuaternary(rnaRaw, rnaCluster, vertices, gene,
                                   veloGraph = rnaVelo, interactive = TRUE,
                                   title = "All cells"))
    expect_message(plotQuaternary(rnaRaw, rnaCluster, vertices, gene, byCluster = TRUE,
                                  interactive = TRUE),
                   "can show desired cluster\\(s\\) by clicking on the legend")
})

test_that("Test quaternary - dense", {
    rnaRawSub <- as.matrix(rnaRaw[gene,])
    p <- plotQuaternary(rnaRawSub, rnaCluster, vertices)
    expect_s3_class(p, "plotly")
})

test_that("Test quaternary GIF", {
    grouping <- list(A = c("ORT"),
                     B = c("RE", "OS"),
                     C = "CH",
                     D = "Stem")
    expect_error(writeQuaternaryGIF(rnaRaw, rnaCluster),
                 "Please set `...` arguments with explicit argument names")
    expect_error(writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                                    vertices = grouping, features = gene,
                                    cluster = letters),
                 "Can only generate GIF for one cluster at a time")
    expect_error(writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                                    vertices = grouping, features = gene,
                                    cluster = "a"),
                 "Specified `cluster`")
    expect_error(writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                                    features = gene, vertices = grouping,
                                    fps = 33),
                 "`FPS` must be a factor of 100")

    skip_on_cran()
    tmpgif <- tempfile(pattern = "test_", fileext = ".gif")
    expect_message(writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                                      vertices = grouping, features = gene,
                                      cluster = "RE", filename = tmpgif,
                                      theta = 10),
                   "These arguments are ignored")
    expect_true(file.exists(tmpgif))
    unlink(tmpgif)

    tmpgif <- tempfile(pattern = "test2_", fileext = ".gif")
    writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster,
                       vertices = grouping, features = gene,
                       filename = tmpgif)
    expect_true(file.exists(tmpgif))
    unlink(tmpgif)

    tmpgif <- tempfile(pattern = "test3_", fileext = ".gif")
    writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster, features = gene,
                       vertices = grouping, useCluster = "Stem",
                       filename = tmpgif)
    expect_true(file.exists(tmpgif))
    unlink(tmpgif)
})
