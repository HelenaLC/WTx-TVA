# dependencies
suppressPackageStartupMessages({
    library(AUCell)
    library(scuttle)
    library(jsonlite)
    library(HDF5Array)
    library(BiocParallel)
    library(SingleCellExperiment)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)
set.seed(241008)

# loading
sce <- readRDS(args[[1]])
lys <- readRDS(args[[2]])

# wrangling
sce <- logNormCounts(sce, BPPARAM=bp)
set <- lapply(lys, intersect, rownames(sce))
sapply(lys, length); sapply(set, length)
length(unique(unlist(set))); length(set)

# analysis
idx <- split(seq(ncol(sce)), sce$fov)
res <- bplapply(idx, BPPARAM=bp, \(.) {
    mtx <- as(logcounts(sce)[, ., drop=FALSE], "dgCMatrix")
    rnk <- AUCell_buildRankings(mtx, plotStats=FALSE, verbose=FALSE)
    auc <- AUCell_calcAUC(set, rnk, verbose=FALSE)
})
res <- do.call(cbind, res)
res <- res[, colnames(sce)]
colData(res) <- colData(sce)

# saving
rowData(res)$set <- set
saveRDS(res, args[[3]])
