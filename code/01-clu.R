# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(igraph)
    library(BiocParallel)
    library(BiocSingular)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)
set.seed(250416)

# loading
sce <- readRDS(args[[1]])

# reduction
pcs <- calculatePCA(
    assay(sce), ncomponents=20,
    BSPARAM=RandomParam(), BPPARAM=bp)
reducedDim(sce, "PCA") <- pcs

# clustering
g <- buildSNNGraph(sce, use.dimred="PCA", type="jaccard", k=20, BPPARAM=bp)
k <- cluster_leiden(g, objective_function="modularity", resolution=0.5)
table(sce$clu <- factor(k <- k$membership, seq_along(unique(k))))

# saving
saveRDS(sce, args[[2]])