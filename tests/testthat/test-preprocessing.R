library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
vertices <- c("OS", "RE")

test_that("Raw count detection", {
    m.double <- rnaRaw
    m.double@x[1] <- m.double@x[1] + 0.1
    expect_false(scPlotSimplex:::is.rawCounts(m.double))
    m.double <- as.matrix(m.double)
    expect_false(scPlotSimplex:::is.rawCounts(m.double))
    expect_false(scPlotSimplex:::is.rawCounts("m.double"))
})

test_that("Test preprocessing - sparse", {
    expect_error(colNormalize("hello"), "Input matrix of class")
    rnaNorm <- colNormalize(rnaRaw)
    expect_equal(sum(colSums(rnaNorm)), ncol(rnaRaw))
    rnaLog <- colNormalize(rnaRaw, scaleFactor = 1e4, log = TRUE)
    expect_s4_class(rnaLog, "dgCMatrix")
})

test_that("Test preprocessing - dense", {
    rnaRaw <- as.matrix(rnaRaw)
    rnaNorm <- colNormalize(rnaRaw)
    expect_equal(sum(colSums(rnaNorm)), ncol(rnaRaw))
    expect_identical(class(rnaNorm), c("matrix", "array"))
})

test_that("Test wilcoxon - sparse", {
    expect_error(selectTopFeatures(rnaRaw, "hello", vertices = vertices),
                 "number of columns of x does not match length")
    expect_error(selectTopFeatures(rnaRaw, rep("a", ncol(rnaRaw)),
                                   vertices = vertices),
                 "Specified vertex clusters are not all found in the cluster")
    expect_error(selectTopFeatures(rnaRaw, rnaCluster,
                                   vertices = list(a = c("OS", "RE"),
                                                   b = c("OS", "CH"))),
                 "Overlap found between elements in list")
    rnaNorm <- colNormalize(rnaRaw)
    expect_warning(selectTopFeatures(rnaNorm, rnaCluster, vertices),
                   "Input matrix is not raw counts")
    gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)
    expect_equal(length(gene), 60)

    st <- selectTopFeatures(rnaRaw, rnaCluster, vertices, returnStats = TRUE)
    expect_equal(nrow(st), length(levels(rnaCluster))*nrow(rnaRaw))
})

test_that("Test wilcoxon - dense", {
    rnaRaw <- as.matrix(rnaRaw)

    gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)
    expect_equal(length(gene), 60)

    st <- selectTopFeatures(rnaRaw, rnaCluster, vertices, returnStats = TRUE)
    expect_equal(nrow(st), length(levels(rnaCluster))*nrow(rnaRaw))
})
