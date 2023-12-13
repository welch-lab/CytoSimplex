if (requireNamespace("Seurat", quietly = TRUE) &&
    suppressWarnings(requireNamespace("SingleCellExperiment", quietly = TRUE)))
{
    Sys.setenv("OMP_THREAD_LIMIT" = 2)

    library(testthat)
    library(CytoSimplex)
    library(Matrix)

    data("rnaRaw", package = "CytoSimplex")
    data("rnaCluster", package = "CytoSimplex")
    data("rnaVelo", package = "CytoSimplex")
    vertices <- c("OS", "RE", "CH", "ORT")

    library(Seurat)
    suppressWarnings({srt <- CreateSeuratObject(rnaRaw)})
    Idents(srt) <- rnaCluster
    srt$cellType <- rnaCluster
    srt <- colNormalize(srt)

    suppressPackageStartupMessages({
        suppressWarnings({
            library(SingleCellExperiment)
        })
    })
    sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
    colLabels(sce) <- rnaCluster
    sce <- colNormalize(sce)

    test_that("normalization", {
        srt <- colNormalize(srt)
        expect_equal(sum(colSums(LayerData(srt, layer = "data", assay = "RNA"))),
                     ncol(rnaRaw))

        sce <- colNormalize(sce)
        expect_true("normcounts" %in% assayNames(sce))

        sce <- colNormalize(sce, 1e4, TRUE)
        expect_true("logcounts" %in% assayNames(sce))
    })


    test_that("wilcoxon", {
        gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices[1:2])
        gene.srt <- selectTopFeatures(srt, vertices = vertices[1:2])
        gene.srt2 <- selectTopFeatures(srt, clusterVar = "cellType", vertices[1:2])
        gene.sce <- selectTopFeatures(sce, vertices = vertices[1:2])
        gene.sce2 <- selectTopFeatures(sce, clusterVar = "label", vertices[1:2])

        expect_identical(gene, gene.srt)
        expect_identical(gene, gene.srt2)
        expect_identical(gene, gene.sce)
        expect_identical(gene, gene.sce2)
    })

    test_that("binary", {
        gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices[1:2])
        p1 <- plotBinary(srt, vertices = vertices[1:2], features = gene)
        p2 <- plotBinary(sce, vertices = vertices[1:2], features = gene)
        expect_s3_class(p1, "ggplot")
        expect_s3_class(p2, "ggplot")
        expect_no_warning(plotBinary(srt, slot = "data", vertices = vertices[1:2],
                                     features = gene))
        expect_no_warning(plotBinary(sce, assay.type = "normcounts",
                                     vertices = vertices[1:2], features = gene))
    })

    test_that("ternary", {
        gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices[1:3])
        p1 <- plotTernary(srt, vertices = vertices[1:3], features = gene)
        p2 <- plotTernary(sce, vertices = vertices[1:3], features = gene)
        expect_s3_class(p1, "ggplot")
        expect_s3_class(p2, "ggplot")
        expect_no_warning(plotTernary(srt, slot = "data", vertices = vertices[1:3],
                                      features = gene))
        expect_no_warning(plotTernary(sce, assay.type = "normcounts",
                                      vertices = vertices[1:3], features = gene))
    })

    test_that("quaternary", {
        gene <- selectTopFeatures(rnaRaw, rnaCluster, vertices)
        p1 <- plotQuaternary(srt, vertices = vertices, features = gene)
        p2 <- plotQuaternary(sce, vertices = vertices, features = gene)
        expect_s3_class(p1, "plist")
        expect_s3_class(p2, "plist")
        expect_no_warning(plotQuaternary(srt, slot = "data", vertices = vertices,
                                         features = gene))
        expect_no_warning(plotQuaternary(sce, assay.type = "normcounts",
                                         vertices = vertices, features = gene))
    })
}
