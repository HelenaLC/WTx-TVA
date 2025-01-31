# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
set.seed(241101)
sce <- readRDS(args[[1]])
pbs <- readRDS(args[[2]])

# clustering
gs <- intersect(rownames(sce), rownames(pbs))
table(ex <- colSums(counts(sce[gs, ]) > 0) < 10)
ist <- .ist(sce[, !ex], gs=gs, pbs=pbs, nk=0, ns=c(2e4, 4e4, 2e5))
ist$clust <- factor(ist$clust, sort(unique(colnames(pbs))))

# saving
saveRDS(ist, args[[3]])