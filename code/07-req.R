# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(scuttle)
    library(slingshot)
    library(BiocParallel)
})

# setup
th <- as.numeric(args$ths)
bp <- MulticoreParam(th)
set.seed(250423)

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# reduction
cs <- names(ks <- ist$clust)
gs <- rownames(ist$profiles)
sce <- logNormCounts(sce[, cs])
sce <- runPCA(sce, subset_row=gs, BPPARAM=bp)
sce <- runUMAP(sce, dimred="PCA", BPPARAM=bp)

# trajectory
sce <- slingshot(sce, 
    clusterLabels=ks, 
    start.clus="stem",
    reducedDim="PCA")

# saving
base::saveRDS(sce, args[[3]])