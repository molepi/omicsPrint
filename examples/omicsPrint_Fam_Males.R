
##try to filter probes like in vignette

getGEO <- function(..., tmax = 5, verbose = TRUE){
    require(GEOquery)
    gse <- NULL
    tryMax <- tmax
    tries <- 0
    if(verbose)
        message("Downloading data from GEO ...")
    while(!inherits(gse,'list') | tries >= tryMax) {
        tries <- tries + 1
        gse <- tryCatch(
            GEOquery::getGEO(...),
            error = function(e) {
                Sys.sleep(5)
                if(verbose) {
                    message(e)
                    message(' ... retrying')
                }
            })
    }
    return(gse)
}

##library(GEOquery)
##gse <- getGEO("GSE73412", GSEMatrix = TRUE)




library(GEOquery)
gse <- getGEO("GSE73412", GSEMatrix = T)

library(omicsPrint)
library(SummarizedExperiment)
se <- makeSummarizedExperimentFromExpressionSet(gse$GSE73412_series_matrix.txt.gz)
r <- expand.grid(idx=colnames(se), idy=colnames(se))
r$Xfam <- substr(as.character(colData(se)[r$idx, "characteristics_ch1.3"]),
                 21, 29)
r$Yfam <- substr(as.character(colData(se)[r$idy, "characteristics_ch1.3"]),
                 21, 29)
r$Xrole <- substr(as.character(colData(se)[r$idx, "characteristics_ch1.3"]),
                  31, 50)
r$Yrole <- substr(as.character(colData(se)[r$idy, "characteristics_ch1.3"]),
                  31, 50)
r[r == " N/A"] <- ""

r$relationship <- "unrelated"

fun <- function(x) {
    if (x["idx"] == x["idy"]){
        return("identical")
    }

    if (x["Xfam"] == "" || x["Yfam"] == ""){
        return("unrelated")
    }
    ##print(x)
    if (x["Yfam"] == x["Xfam"]){
        roles <- x[c("Xrole", "Yrole")]
        ##print(roles)
        if ("(son)" %in% roles && "(father)" %in% roles){
            return("father/son")
        } else if ("(son)" %in% roles && "(grandfather)" %in% roles) {
            return("grandfather/grandson")
        } else if ("(son)" %in% roles && "(uncle)" %in% roles){
            return("uncle/nephew")
        } else if ("(father)" %in% roles && "(grandfather)" %in% roles) {
            return("father/son")
        } else if ("(father)" %in% roles && "(uncle)" %in% roles) {
            return("brother")
        } else if ("(grandfather)" %in% roles && "(uncle)" %in% roles) {
            return("father/son")
        } else if (sum(roles == "(uncle)") == 2 || sum(roles == "(son)") == 2) {
            return("brother")
        }
    } else {
        return("unrelated")
    }
}

r$relationship <- apply(r, 1, fun)

data(hm450.manifest.pop.GoNL)
cpgs <- names(hm450.manifest.pop.GoNL[
    mcols(hm450.manifest.pop.GoNL)$MASK.snp5.EAS])

se2 <- se[cpgs,]

se3 <- se[, se$characteristics_ch1.3 != "family relationship: N/A"]
se4 <- se3[cpgs,]

gt <- beta2genotype(se4, assayName = "exprs")
data <- alleleSharing(gt, relations = r, rel.col = "relationship", verbose = TRUE)
mismatches <- inferRelations(data)
mismatches


## Inspect Hardy-Weinberg filtering

alpha <- seq(0, 1, by = .01)
nsnps <- numeric(length(alpha))
for(i in 1:length(alpha)) {
    x <- omicsPrint:::.pruning(gt, callRate=0.95, coverageRate=2/3, alpha = alpha[i], verbose=TRUE)
    nsnps[i] <- nrow(x)
}

n <- nrow(gt)


par(mar=c(5,4,6,2))
plot(alpha/n, n-nsnps, type = "b",
     log = "x", panel.first = grid(col=1, equilogs=FALSE),
     xlab = "Significant threshold (Bonferonni corrected)",
     ylab = "Number of SNPs",
     main = "Number of SNPs violating Hardy-Weinberg principle\n(at different significance levels)")
at <- axTicks(1)
axis(3, at, labels=at*n)


