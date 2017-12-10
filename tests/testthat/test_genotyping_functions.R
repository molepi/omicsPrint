context("Genotyping helper functions")

test_that("phasing recovers swapped SNPs", {

    ##generate some data
    x <- matrix(c(1,1,1, 2, 2, 2, 3, 3, 3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))

    ##generate another data set with one swapped SNP
    y <- x
    y[2,] <- abs(2 - y[2,])     ##swapping

    ##constructe hashed relations
    rHash <- .hashRelations(.constructRelations(xnames = colnames(x),
                                                ynames = colnames(y)))
    ##do the testing
    expect_equal(.phasing(x, y, rHash), x)
})


test_that("pruning properly drops SNPs", {

    ##generate some data with missing values
    x <- matrix(c(NA,1,1,
                  0,2,2,
                  3,3,3,
                  1,2,2,
                  2,2,3,
                  3,3,2), 6,3, byrow=TRUE,
                dimnames=list(paste0("SNP", 1:6),
                              paste0("sample", 1:3)))
    ##do the testing
    expect_equal(.pruning(x, verbose=FALSE), x[4:6,])

})


test_that("allele sharing algorithm", {

    ##generate some data
    x <- matrix(c(1,1,1, 2, 2, 2, 3, 3, 3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))

    allelesSharedSq <- .square(x, x, verbose = FALSE)
    allelesSharedRect <- .rectangular(x, x, verbose = FALSE)

    ##do the testing
    expect_equal(nrow(allelesSharedSq), 6) ##expected dimensions
    expect_equal(nrow(allelesSharedRect), 9)

    expect_equal(allelesSharedSq$mean[1:3], c(2,1,0)) ##expected mean IBS
    expect_equal(allelesSharedRect$mean[1:3], c(2,1,0))

    expect_equal(allelesSharedSq$var[1:3], c(0,0,0)) ##expected variance IBS
    expect_equal(allelesSharedRect$var[1:3], c(0,0,0))

})

test_that("SummarizedExperiment input alleleSharing", {
    x <- matrix(c(1,1,1, 2, 2, 2, 3, 3, 3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))

    Xse <- SummarizedExperiment::SummarizedExperiment(assays = list(a = x, b=x))
    Yse <- SummarizedExperiment::SummarizedExperiment(list(a=x))

    Xexpected <- alleleSharing(x)
    XYexpected <- alleleSharing(x, x)

    expect_equal(alleleSharing(Xse, assayNameX = "a"), Xexpected)
    expect_equal(alleleSharing(Xse, Yse, assayNameX = "a", assayNameY = "a"), XYexpected)
    expect_equal(alleleSharing(Xse, assayNameX = "a", assayNameY = "b"), XYexpected)
})

test_that("MultiAssayExperiment input alleleSharing", {
    x <- matrix(c(1,1,1, 2, 2, 2, 3, 3, 3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))

    pheno <- data.frame(id = 1:3, type = c("a", "a", "b"),
                        sex = c("M", "F", "M"),
                        row.names = c("sample1", "sample2", "sample3"))
    expList <- list(exp1 = x, exp2 = x)
    map1 <- data.frame(primary = c("sample1", "sample2", "sample3"),
                       colname = c("sample1", "sample2", "sample3"))
    map2 <- data.frame(primary = c("sample1", "sample2", "sample3"),
                       colname = c("sample1", "sample2", "sample3"))
    sampMap <- MultiAssayExperiment::listToMap(list(exp1 = map1, exp2 = map2))
    sampMap$primary <- as.character(sampMap$primary)
    sampMap$colname <- as.character(sampMap$colname)
    maeX <- MultiAssayExperiment::MultiAssayExperiment(expList, pheno, sampMap)
    maeY <- MultiAssayExperiment::MultiAssayExperiment(expList, pheno, sampMap)

    Xexpected <- alleleSharing(x)
    XYexpected <- alleleSharing(x, x)

    expect_equal(alleleSharing(maeX, assayNameX = "exp1"), Xexpected)
    expect_equal(alleleSharing(maeX, maeY, assayNameX = "exp1", assayNameY = "exp2"), XYexpected)
    expect_equal(alleleSharing(maeX, assayNameX = "exp1", assayNameY = "exp2"), XYexpected)
})

test_that("RaggedExperiment input alleleSharing", {
    x <- matrix(c(1,1,1, 2, 2, 2, 3, 3, 3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))
    sample1 <- GenomicRanges::GRanges(
                                  c(GENEA = "chr1:1-10:-", GENEB = "chr2:15-18:+", GENEC = "chr2:11-18:+"),
                                  score = c(1,1,1))
    sample2 <- GenomicRanges::GRanges(
                                  c(GENEA = "chr1:1-10:-", GENEC = "chr2:11-18:+", GENEB = "chr2:15-18:+"),
                                  score = c(2,2,2))
    sample3 <- GenomicRanges::GRanges(
                                  c(GENEA = "chr1:1-10:-", GENEC = "chr2:11-18:+", GENEB = "chr2:15-18:+"),
                                  score =c(3,3,3))
    colDat <- data.frame(id = 1:3)

    Xragexp <- RaggedExperiment::RaggedExperiment(sample1 = sample1,
                                                  sample2 = sample2,
                                                  sample3 = sample3,
                                                  colData = colDat)
    Yragexp <- RaggedExperiment::RaggedExperiment(sample1 = sample1,
                                                  sample2 = sample2,
                                                  sample3 = sample3,
                                                  colData = colDat)

    Xexpected <- alleleSharing(x)
    XYexpected <- alleleSharing(x, x)

    expect_equal(alleleSharing(Xragexp), Xexpected)
    expect_equal(alleleSharing(Xragexp, Yragexp), XYexpected)

})

test_that("SE input beta2genotype", {
    x <- matrix(rnorm(9,0.5,0.3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(a=x))

    expectected <- beta2genotype(x)
    expect_equal(beta2genotype(se, assayName = "a"), expectected)
})

test_that("MAE input beta2genotype", {
    x <- matrix(rnorm(9,0.5,0.3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))

    pheno <- data.frame(id = 1:3, type = c("a", "a", "b"),
                        sex = c("M", "F", "M"),
                        row.names = c("sample1", "sample2", "sample3"))
    expList <- list(exp1 = x)
    map1 <- data.frame(primary = c("sample1", "sample2", "sample3"),
                       colname = c("sample1", "sample2", "sample3"))
    sampMap <- MultiAssayExperiment::listToMap(list(exp1 = map1))
    sampMap$colname <- as.character(sampMap$colname)
    sampMap$primary <- as.character(sampMap$primary)
    maeX <- MultiAssayExperiment::MultiAssayExperiment(expList, pheno, sampMap)

    expectected <- beta2genotype(x)
    expect_equal(beta2genotype(x, assayName = "a"), expectected)
})
