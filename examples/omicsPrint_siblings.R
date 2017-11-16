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

## Inspect Hardy-Weinberg filtering

alpha <- seq(0, 1, by = .01)
nsnps <- numeric(length(alpha))
for(i in 1:length(alpha)) {
    x <- omicsPrint:::.pruning(gt, callRate=0.95, coverageRate=2/3, alpha = alpha[i], verbose=TRUE)
    nsnps[i] <- nrow(x)
}

nrow(gt)


pdf("InspectHWfiltering.pdf")
par(mar=c(5,4,6,2))
plot(alpha/1000, 1002-nsnps, type = "b",
     log = "x", panel.first=grid(col=1, equilogs=FALSE),
     xlab = "Significant threshold (Bonferonni corrected)",
     ylab = "Number of SNPs",
     main = "Number of SNPs violating Hardy-Weinberg principle\n(at different significance levels)")
at <- axTicks(1)
axis(3, at, labels=at*1000)
dev.off()

