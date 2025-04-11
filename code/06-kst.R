# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
set.seed(250326)
pbs <- readRDS(args[[1]])
sce <- readRDS(args[[2]])
ist <- readRDS(args[[3]])

# wrangling
ids <- "EE|entero|goblet"
idx <- grep(ids, ist$clust, value=TRUE)
ncol(sce <- sce[, names(idx)])

# clustering
length(gs <- intersect(rownames(sce), rownames(pbs)))
table(ex <- colSums(counts(sce[gs, ]) > 0) < 10)
ist <- .ist(sce[, !ex], gs=gs, pbs=pbs, nk=0)
table(ist$clust <- factor(ist$clust))

# saving
saveRDS(ist, args[[4]])