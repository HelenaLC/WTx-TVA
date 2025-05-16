# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpmisc)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
set.seed(240512)
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])
pbs <- readRDS(args[[3]])$profiles

# clustering
sce <- sce[, names(which(ist$clust == "epi"))]
gs <- intersect(rownames(sce), rownames(pbs))
table(ex <- colSums(counts(sce[gs, ]) > 0) < 20)
jst <- .ist(sce[, !ex], gs=gs, pbs=pbs, nk=seq(1, 5))
jst$clust <- factor(jst$clust, colnames(pbs))

# saving
saveRDS(jst, args[[4]])
