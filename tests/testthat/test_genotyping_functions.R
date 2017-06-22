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
    x <- matrix(c(NA,1,NA, 2, 2, 2, 3, 3, 3), 3,3,
                dimnames=list(paste0("SNP", 1:3),
                              paste0("sample", 1:3)))
    ##do the testing
    expect_equal(.pruning(x, verbose=FALSE),
                 matrix(c(2, 3), 1, 2,
                        dimnames=list("SNP2", c("sample2", "sample3")))
                 )
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
