.rectangular <- function(x, y, verbose = TRUE){

    if( verbose )
        message("Running `rectangular` IBS algorithm!")

    if( nrow(x) != nrow(y) )
        stop("Dimension mismatch!")

    if( is.null(colnames(x)) | is.null(colnames(y)) )
        stop("Colnames should be provided!")

    na.rm <- sum(is.na(x)) > 0 | sum(is.na(y)) > 0

    ##calculates mean and variance of IBS between all pairs of x and y
    N <- ncol(x)
    M <- ncol(y)
    indices <- seq_len(N)
    mn <- s2 <- numeric(length = N*M)
    for(j in seq_len(M)) {
        ibs <- 2 - abs(x - y[,j])
        mn[indices + N*(j-1)] <- colMeans(ibs, na.rm = na.rm)
        s2[indices + N*(j-1)] <- colVars(as.matrix(ibs), na.rm = na.rm)

        if( verbose & (j %% 100 == 0 | j == 1) )
            message(j*N, " of ", N*M, " (", round(100*j/M, 2), "%) ...")

    }
    data.frame(mean = mn, var = s2,
               colnames.x = rep(colnames(x), M), colnames.y = rep(colnames(y), each=N))
}

.square <- function(x, y, verbose = TRUE){

    if( verbose )
        message("Running `square` IBS algorithm!")

    if( (ncol(x) != ncol(y)) | (nrow(x) != nrow(y)) )
        stop("Dimension mismatch!")

    if( is.null(colnames(x)) | is.null(colnames(y)) )
        stop("Colnames should be provided!")

    na.rm <- sum(is.na(x)) > 0 | sum(is.na(y)) > 0

    ##calculates mean and variance of IBS between all unique pairs of x (and y)
    N <- ncol(x)
    k <- 1
    mn <- s2 <- numeric(length = N*(N+1)/2)
    for(j in seq_len(N)) {
        ibs <- 2 - abs(x[, j:N, drop=FALSE] - y[,j])
        mn[k + 0:(N - j)] <- colMeans(ibs, na.rm = na.rm)
        s2[k + 0:(N - j)] <- colVars(as.matrix(ibs), na.rm = na.rm)
        k <- k + 1 +  N - j
        if( verbose & (j %% 100 == 0 | j == 1) )
            message(k, " of ", N*(N+1)/2, " (", round(100*k/(N*(N+1)/2), 2), "%) ...")
    }

    data.frame(mean = mn, var = s2,
               colnames.x = unlist(sapply(1:N, function(k) colnames(x)[k:N])), colnames.y = rep(colnames(y), N:1))
}


.phasing <- function(x, y, rHash) {
    ##relabel those in x according to those in y

    ##relabelling is based on the idea that snp's in x should be
    ##positively correlated with those in y robust against NA's and
    ##outliers

    ##FIX level `identical`

    identical <- names(grep("identical", as.list(rHash), value = TRUE))

    colnames.x <- gsub(":.*$", "", identical)
    colnames.y <- gsub("^.*:", "", identical)
    colnames.x <- colnames.x[!grepl("NA", colnames.y)]
    colnames.y <- colnames.y[!grepl("NA", colnames.y)]

    keep <- colnames.x %in% colnames(x) & colnames.y %in% colnames(y)
    if( sum(keep) == 0 )
        stop("No overlapping samples!")

    colnames.x <- colnames.x[keep]
    colnames.y <- colnames.y[keep]

    midx <- match(colnames.x, colnames(x))
    midy <- match(colnames.y, colnames(y))

    signs <- sign(unlist(lapply(seq_len(nrow(x)), function(i)
        cov(x[i,midx], y[i,midy], use = "complete.obs", method = "spearman"))))

    for(i in seq_len(nrow(x))) {
        if( signs[i] < 0 ) {
            xi <- x[i,]
            x[i, xi==1] <- 3
            x[i, xi==3] <- 1
        }
    }
    x
}

.hashRelations <- function(relations, idx.col = "idx", idy.col = "idy", rel.col = "relation_type"){
    hash <- new.env()
    keys <- paste(relations[, idx.col], relations[, idy.col], sep=":")
    values <- as.character(relations[, rel.col])

    ##are there inconsistent relations
    if( any(tapply(values, keys, function(x) length(unique(x)) > 1)) )
        stop("Nonconsistent relations!")

    ##redundant relations are automatically removed!
    tmp <- mapply(assign, keys, values, MoreArgs = list(envir = hash))
    hash
}

.constructRelations <- function(xnames, ynames, idx.col = "idx", idy.col = "idy", rel.col = "relation_type") {
    if( !is.null(ynames) )
        relations <- expand.grid(idx = xnames, idy = ynames, stringsAsFactors = FALSE)
    else
        relations <- expand.grid(idx = xnames, idy = xnames)
    relations[,rel.col] <- "unrelated"
    identical <- relations[,idx.col] == relations[,idy.col]

    relations[identical, rel.col] <- "identical"
    relations <- relations[relations[,rel.col] != "unrelated",]
    rownames(relations) <- 1:nrow(relations)
    relations
}


.pruning <- function(x, callRate = 0.95, coverageRate = 2/3, verbose = TRUE) {

    nsnps <- nrow(x)
    nsamples <- ncol(x)

    ##drop SNPs difficult to call
    calledSNPs <- apply(x, 1, function(x) sum(!is.na(x))/nsamples)

    if( verbose )
        message("There are ", sum(calledSNPs <= callRate),
                " SNP dropped because of low call rate!")

    ##if the coverage of called SNPs is not larger then coverageRate do not calculate IBS
    coverage <- apply(x, 2, function(x) sum(!is.na(x))/nsnps)

    if( verbose )
        message("There are ", sum(coverage < coverageRate),
                " sample set to NA because too little SNPs called!")

    x[calledSNPs > callRate, coverage >= coverageRate, drop=FALSE] ##drop those
}

##' Run the allele sharing algorithm based on ibs
##'
##' calculate mean and variance of identity by state between
##' genotypes, coded as 1,2,3, of all sample pairs either give one
##' omic inferred set of SNPs or two from different omic types.
##'
##' 'phasing' is required if methylation data inferred genotypes are
##' compared with DNA based genotypes, i.e., DNA based genotype is 3
##' whereas methylation is 1. The phasing step will fix this.
##'
##' Notice that there are two algorithms for calculating allele
##' sharing. One for the case one matrix is provided and one for the
##' case two, x and y, are provided. If one is provided the algorithm
##' takes in account the symmetric relations between pairs i.e. x12 =
##' x21 etc.
##'
##' To improve the lookup of relations, which can be millions of say
##' 1000 samples are provided, a hash-table is created from the
##' provided data.frame with relations.
##'
##' @title allele sharing based on ibs
##' @param x genotype matrix with row and column names
##' @param y genotype matrix with row and column names
##' @param relations 'data.frame' with relations and their mapping
##'     identifiers, provide columns if different from the default and
##'     beware identifiers should match with colnames of x and y.
##' @param idx.col column name containing mapping identifiers
##' @param idy.col column name containing mapping identifiers
##' @param rel.col column name containing the relations,
##'     i.e. identical, parentoffspring, etc.
##' @param callRate default 0.95 SNPs that are called in less then the
##'     threshold are dropped
##' @param coverageRate default 2/3 samples with less then threshold
##'     SNPs called are set to NA
##' @param phasing default FALSE
##' @param verbose show progress default TRUE
##' @return data.frame with mean and variance ibs between all pairs
##' @author mvaniterson
##' @importFrom matrixStats colVars
##' @importFrom stats cov
##' @export
##' @examples
##' set.seed(12345)
##' beta <- matrix(runif(100*10, 0,1), nrow=100)
##' beta[1:5, 1:5]
##' colnames(beta) <- paste0("sample", 1:10)
##' genotype <- beta2genotype(beta)
##' genotype[1:5, 1:5]
##' data <- alleleSharing(genotype)
##' head(data)
alleleSharing <- function(x, y = NULL, relations = NULL, idx.col = "idx", idy.col = "idy", rel.col = "relation_type", callRate = 0.95, coverageRate = 2/3, phasing = FALSE, verbose = TRUE) {

    if(is.null(colnames(x)))
        stop("Colnames should be given!")

    if(!is.null(y) & is.null(colnames(y)))
        stop("Colnames should be given!")

    if( is.null(relations )) {
        relations <- .constructRelations(xnames = colnames(x), ynames = colnames(y),
                                         idx.col = idx.col, idy.col = idy.col, rel.col = rel.col)
        relations <- relations[!duplicated(relations),]
    }

    if( verbose )
        message("Hash relations")
    rHash <- .hashRelations(relations, idx.col = idx.col, idy.col = idy.col, rel.col = rel.col)

    if( is.null(y) ) {

        x <- .pruning(x, callRate = callRate, coverageRate = coverageRate, verbose = verbose)

        if( verbose )
            message("Using ", nrow(x), " polymorphic SNPs to determine allele sharing.")

        data <- .square(x, x, verbose)
    } else {

        x <- .pruning(x, callRate = callRate, coverageRate = coverageRate, verbose = verbose)
        y <- .pruning(y, callRate = callRate, coverageRate = coverageRate, verbose = verbose)

        rows <- intersect(rownames(x), rownames(y))
        rId <- match(rows, rownames(x))
        x <- x[rId,]
        rId <- match(rows, rownames(y))
        y <- y[rId,]

        if( phasing )
            x <- .phasing(x, y, rHash)

        if( verbose )
            message("Using ", nrow(x), " polymophic SNPs to determine allele sharing.")

        data <- .rectangular(x, y, verbose)
        if( !(any(colnames(x) %in% data$colnames.x) & any(colnames(y) %in% data$colnames.y)) )
            stop("rHash and x or y do not match: probably swap 'x' and 'y'!")
    }

    pairs <- paste(data$colnames.x, data$colnames.y, sep = ":")
    mId <- pairs %in% ls(rHash)
    data$relation <- "unrelated"
    data$relation[mId] <- unlist(mget(pairs[mId], rHash, mode = "character", ifnotfound = list("unrelated")), use.names = FALSE)
    data$relation <- factor(data$relation)
    invisible(data)
}

##' predict mismatches
##'
##' based on all data a classifier is build using Linear Discriminant
##' Analysis and on the same data a prediction is performed in order
##' to detect wrong sample relationships. The assumption is that the
##' majority of sample relations is correct otherwise we could not do
##' this!
##' @title predict mismatches
##' @param data output from allelesharing
##' @param n = 100 default interpolation for showing the
##'     classification boundaries
##' @param plot.it = TRUE default plot classification graph and
##'     returing mismatches otherwise return all
##' @return predicted mismatches
##' @author mvaniterson
##' @importFrom graphics contour legend plot points
##' @importFrom MASS lda
##' @importFrom stats predict
##' @export
##' @examples
##' set.seed(12345)
##' beta <- matrix(runif(100*10, 0,1), nrow=100)
##' beta[1:5, 1:5]
##' colnames(beta) <- paste0("sample", 1:10)
##' genotype <- beta2genotype(beta)
##' genotype[1:5, 1:5]
##' data <- alleleSharing(genotype)
##' head(data)
##' inferRelations(data)
inferRelations <- function(data, n = 100, plot.it = TRUE){

    data <- droplevels(data)
    model <- lda(relation~mean+var, data = data)

    predicted <- predict(model, data)

    data$predicted <- predicted$class

    print(table(`Predicted relation` = data$predicted, `Assumed relation` = data$relation), zero.print = ".")


    if( plot.it ) {

        id <- which(data$predicted != data$relation)

        plot(data[, c("mean", "var")], pch = ".", cex = 3,
             col = as.integer(data$relation),
             xlab = "mean (IBS)", ylab = "variance (IBS)")

        xp <- seq(min(data$mean), max(data$mean), length = n)
        yp <- seq(min(data$var), max(data$var), length = n)
        grid <- expand.grid(mean = xp, var = yp)
        predicted <- predict(model, grid)
        posterior <- predicted$posterior

        if( ncol(posterior) > 2 ) {
            for(k in seq_len(ncol(posterior))) {
                zp <- posterior[, k] - apply(posterior[,-k], 1, max)
                contour(xp, yp, matrix(zp, n), add = TRUE, levels = 0,
                        drawlabels = FALSE, lty = 1, lwd = 2, col = "grey")
            }
        } else {
            zp <- posterior[, 2] - pmax(posterior[, 1])
            contour(xp, yp, matrix(zp, n), add = TRUE, levels = 0,
                    drawlabels = FALSE, lty = 1, lwd = 2, col = "grey")
        }

        points(data[id, c("mean", "var")], pch = ".", cex = 3,
               col = as.integer(data$relation[id]))
        legend("topright", paste("assumed", levels(data$relation)),
               col = 1:nlevels(data$relation), pch = 15, bty = "n")

        return(invisible(data[id,]))
    }
    else
        return(invisible(data))

}

##' convert DNA methylation beta-value to inferred genotypes
##'
##' Using kmeans unsupervised clustering to infer genotypes based on
##' idea's from Leonard Schalkwyk; wateRmelon packages.
##'
##' 'minSep' and 'minSize' ensure good clusters are found.
##' This function is similar to the gaphunter approach implemented in minfi.
##' @title converts beta-values to genotypes (1, 2 and 3)
##' @param betas beta matrix of probes possibly affected SNPs
##' @param na.rm TRUE drop cpg for which no clustering was observed
##' @param minSep minimal separation between clusters
##' @param minSize size of smallest cluster (in percentage)
##' @param centers center of clusters, defaults to 0.2, 0.5, 0.8.
##' @return matrix with genotypes
##' @author mvaniterson
##' @importFrom stats kmeans
##' @export
##' @examples
##' set.seed(12345)
##' beta <- matrix(runif(100*10, 0,1), nrow=100)
##' beta[1:5, 1:5]
##' genotype <- beta2genotype(beta)
##' genotype[1:5, 1:5]
beta2genotype <- function (betas, na.rm = TRUE, minSep = 0.25, minSize = 5, centers = c(0.2, 0.5, 0.8))
{
    genotypes <- apply(betas, 1, function(x) {

        km <- try(kmeans(as.numeric(x), centers), silent = TRUE)

        if ( !inherits(km, "try-error") ) {
            if ( all(abs(rep(km$centers, 3) - rep(km$centers,
                                                  each = 3))[-c(1, 5, 9)] > minSep) ) {
                if ( 100 * min(as.numeric(table(km$cluster)))/length(x) >
                     minSize )
                    return(km$cluster)
            }
        }
        return(rep(NA, length(x)))
    })

    genotypes <- t(genotypes)
    colnames(genotypes) <- colnames(betas)
    rownames(genotypes) <- rownames(betas)
    if ( na.rm ) {
        nas <- apply(genotypes, 1, anyNA)
        genotypes <- genotypes[!nas, ]
    }
    genotypes
}
