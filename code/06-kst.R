# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
set.seed(250419)
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])
pbs <- readRDS(args[[3]])

# subsetting
gs <- intersect(rownames(sce), rownames(pbs))
cs <- intersect(
    names(which(ist$clust == "epi")),
    colnames(sce)[grep("REF", sce$typ)])
dim(sce <- sce[gs, cs])

# clustering
table(ex <- colSums(counts(sce) > 0) < 10)
ist <- .ist(sce[, !ex], gs=gs, pbs=pbs, nk=0)
table(ist$clust <- factor(ist$clust))

# saving
saveRDS(ist, args[[4]])