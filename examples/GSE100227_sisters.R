
d <- read_tsv("~/Downloads/GSE100227_normalized_data.txt", progress = T)
d <- d[, ! grepl("Detection", colnames(d))]
d <- as.data.frame(d)
rownames(d) <- d$ID_REF
d <- d[, !colnames(d) == "ID_REF"]



library(GEOquery)
library(SummarizedExperiment)
##gset <- getGEO("GSE100227", GSEMatrix=TRUE, AnnotGPL=FALSE)
##this is more stable somehow?
file <- tempfile()
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100227/suppl/GSE100227_normalized_data.txt.gz", file)
gse <- getGEO(filename=file)
gse
se <- makeSummarizedExperimentFromExpressionSet(gse)
se

r <- expand.grid(idx=colnames(se), idy=colnames(se))

r$Xfam <- colData(se)[r$idx, "characteristics_ch1.1"]
r$Yfam <- colData(se)[r$idy, "characteristics_ch1.1"]
r$Xrole <- colData(se)[r$idx, "characteristics_ch1.2"]
r$Yrole <- colData(se)[r$idy, "characteristics_ch1.2"]

fun <- function(x){
    if (x["idx"] == x["idy"]) {
        return("identical")
    } else if (x["Xfam"] != x["Yfam"]) {
        return("unrelated")
    } else if (x["Xrole"] == x["Yrole"]) {
        return(x["Xrole"])
    } else {
        return("Sister")
    }
    return("hello")
}

r$relation_type <- apply(r, 1, fun)
head(r)

data(hm450.manifest.pop.GoNL)
cpgs <- names(hm450.manifest.pop.GoNL[
        mcols(hm450.manifest.pop.GoNL)$MASK.snp5.EUR])

se2 <- se[cpgs,]

gt <- beta2genotype(d2, minSize = 7)
data <- alleleSharing(gt, relations = r2, verbose = TRUE)
mismatches <- inferRelations(data)

#------
r2 <- r
r2[r2$relation_type == "DZ", "relation_type"] <- "Sister"
#------

.hardyweinberg <- function(x, alpha=0.05){
    
    obs <- rowTabulates(matrix(as.integer(x), nrow=nrow(x)))
    n <- ncol(x)
    p <- (2*obs[,1] + obs[,2])/(2*n)
    q <- 1 - p
    exp <- n*cbind(p^2, 2*p*q, q^2)
    
    stat <- rowSums(((obs - exp)^2)/exp)
    
    ##pval <- pchisq(stat, df = 1)
    
    ##true in Hardy-Weinberg equilibrium
    inEquilibrium <- stat < qchisq(1 - alpha/nrow(x), df = 1)
    inEquilibrium[is.na(inEquilibrium)] <- FALSE
    inEquilibrium
}


pairs <- r[r$relation_type == "MZ",]

fun2 <- function(x, pairs){
    out <- 0
    for (twins in rownames(pairs)) {
        idx <- pairs[twins, "idx"]
        idy <- pairs[twins, "idy"]
        if (is.na(x[idx]) || is.na(x[idy])){
            message("here")
        } else if (x[idx] == x[idy]  ) {
            out <- out + 1
        }
    }
    out
}

d3 <- apply(gt, 1, fun2, pairs)
d4 <- as.data.frame(d3)
d4$yay <- .hardyweinberg(gt)
