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

set.seed(12345)
sizes <- c(5, 10, 25, 50, 100, 250, 500, 1000)

nrep <- 25
miss <- matrix(0, length(sizes), nrep)
cntr = 0
for(i in 1:length(sizes)){
    for(j in 1:nrep)     {
        cntr = cntr +1
        indices <- sample(1:nrow(gt), sizes[i])
        data <- alleleSharing(gt[indices,], relations = r)
        predictions <- inferRelations(data, plot.it=FALSE)
        miss[i, j] <- nrow(predictions) - sum(diag(table(predictions$predicted, predictions$relation)))
    }
}

pdf("missclassification.pdf")
rownames(miss) <- sizes
boxplot(t(miss),
        ylab="Number of missclassified relations (out of 666)", xlab="Number of SNPs")
dev.off()
