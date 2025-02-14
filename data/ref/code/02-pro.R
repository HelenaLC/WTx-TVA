# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(scran)
    library(scater)
    library(igraph)
})

# loading
set.seed(240924)
sce <- readRDS(args[[1]])

# selection
sce <- logNormCounts(sce)
tbl <- modelGeneVar(sce)
sel <- getTopHVGs(tbl, n=4e3)
rowData(sce)$sel <- rownames(sce) %in% sel

# reduction
sce <- runPCA(sce, subset_row=sel, ncomponents=20)
sce <- runUMAP(sce, dimred="PCA")

# clustering
nng <- buildSNNGraph(sce, use.dimred="PCA", type="jaccard", k=10)
ids <- cluster_louvain(nng, resolution=0.8)$membership
table(sce$kid_lv1 <- factor(ids, seq_along(unique(ids))))

# saving
saveRDS(sce, args[[2]])