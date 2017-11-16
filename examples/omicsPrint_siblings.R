library(omicsPrint)
library(GEOquery)
library(SummarizedExperiment)


gse <- getGEO("GSE102177")
se <- makeSummarizedExperimentFromExpressionSet(gse$GSE102177_series_matrix.txt.gz)

r <- expand.grid(idx=colnames(se), idy=colnames(se))
r$Xfam <- colData(se)[r$idx, "characteristics_ch1.3"]
r$Yfam <- colData(se)[r$idy, "characteristics_ch1.3"]

fun <- function(x){
    if (x["idx"] == x["idy"]){
        return("identical")
    } else if (x["Xfam"] == x["Yfam"]){
        return("sibling")
    } else {
        return("unrelated")
    }
}

r$relation_type <- apply(r, 1, fun)
r

gt <- beta2genotype(se, assayName = "exprs")
data <- alleleSharing(gt, relations = r)
mismatches <- inferRelations(data)
mismatches