library(testthat)
library(scPlotSimplex)
library(Matrix)

data("rnaRaw", package = "scPlotSimplex")
data("rnaCluster", package = "scPlotSimplex")
vertices <- c("OS", "RE")

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
    expect_error(selectTopFeatures(rnaRaw, "hello"), "number of columns of x")
    expect_error(selectTopFeatures(rnaRaw, rep("a", ncol(rnaRaw))),
                 "Must have at least 2")

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
