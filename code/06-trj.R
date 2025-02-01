# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# setup
set.seed(250201)

# loading
sce <- readRDS(args[[1]])
roi <- readRDS(args[[2]])

# wrangling
ref <- grepl("REF", roi$typ)
ref <- colnames(roi)[ref]
ids <- rep("-1", ncol(sce))
names(ids) <- colnames(sce)
ids[ref] <- "REF"; table(ids)

# analysis
sce <- slingshot(sce,
    reducedDim="PCA", 
    start.clus="REF",
    clusterLabels=ids, 
    approx_points=100)

# saving
saveRDS(sce, args[[3]])
