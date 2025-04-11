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
set.seed(250405)

# loading
sce <- readRDS(args[[1]])
sig <- lapply(args[[2]], fromJSON)
(n <- sapply(sig, length)); sum(n)

# wrangling
gs <- (gs <- Reduce(c, sig))[unique(names(gs))]
gs <- lapply(gs, intersect, rownames(sce))
range(sapply(gs, length)); length(gs)

# scoring
idx <- split(seq(ncol(sce)), sce$fov)
res <- bplapply(idx, BPPARAM=bp, \(.) {
    mtx <- as(assay(sce)[, ., drop=FALSE], "dgCMatrix")
    rnk <- AUCell_buildRankings(mtx, plotStats=FALSE, verbose=FALSE)
    auc <- AUCell_calcAUC(gs, rnk, aucMaxRank=400, verbose=FALSE)
})
res <- do.call(cbind, res)
res <- res[, colnames(sce)]
colData(res) <- colData(sce)

# saving
sub <- gsub(".*(imm|epi|str).*", "\\1", args[[2]])
sub <- rep.int(sub, sapply(sig, length))
sub <- split(sub, unlist(lapply(sig, names)))
sub <- sub[rownames(res)]
rowData(res)$set <- gs
rowData(res)$sub <- sub
saveRDS(res, args[[3]])
