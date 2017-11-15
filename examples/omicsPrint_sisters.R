library(readr)
library(GEOquery)
library(SummarizedExperiment)
library(omicsPrint)

d <- read_tsv("~/Downloads/GSE100227_normalized_data.txt", progress = T)
d <- d[, ! grepl("Detection", colnames(d))]
d <- as.data.frame(d)
rownames(d) <- d$ID_REF
d <- d[, !colnames(d) == "ID_REF"]


gse <- getGEO("GSE100227")
#save(gse, file="/home/dcats/GSE100227.Rdata")
#load("/home/dcats/GSE100227.Rdata")
samples <- pData(gse$GSE100227_series_matrix.txt.gz)
rownames(samples) <- samples$description


r <- expand.grid(idx=rownames(samples), idy=rownames(samples))
r$Xfam <- substr(samples[r$idx, "characteristics_ch1.1"], 12, 100)
r$Yfam <- substr(samples[r$idy, "characteristics_ch1.1"], 12, 100)
r$Xrole <- substr(samples[r$idx, "characteristics_ch1.2"], 11, 100)
r$Yrole <- substr(samples[r$idy, "characteristics_ch1.2"], 11, 100)

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

data(hm450.manifest.pop.GoNL)
cpgs <- names(hm450.manifest.pop.GoNL[
        mcols(hm450.manifest.pop.GoNL)$MASK.snp5.EUR])

d2 <- d[cpgs,]

gt <- beta2genotype(d)
data <- alleleSharing(gt, relations = r, verbose = TRUE)
mismatches <- inferRelations(data)


#if MZ + MZ -> MZ
#if DZ + DZ -> DZ
#else -> Sisters