# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# setup
set.seed(250209)

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
table(clu <- ifelse(grepl("REF_", kid), "one", "two"))

# analysis
sce <- slingshot(sce, 
    approx_points=100,
    clusterLabels=clu, 
    start.clus="one",
    reducedDim="PCA") 

# saving
base::saveRDS(sce, args[[3]])
