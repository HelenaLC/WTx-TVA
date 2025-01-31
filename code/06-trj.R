# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# setup
set.seed(241202)

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
idx <- match(colnames(sce), names(kid <- ist$clust))
kid <- kid[idx]; kid[is.na(kid)] <- -1

# analysis
sce <- slingshot(sce,
    reducedDim="PCA", 
    clusterLabels=kid, 
    approx_points=100,
    start.clus="entero")

# saving
saveRDS(sce, args[[3]])
