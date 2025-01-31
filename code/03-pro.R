# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(scuttle)
    library(HDF5Array)
    library(BiocParallel)
    library(BiocSingular)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)
set.seed(241130)

# loading
sce <- readRDS(args[[1]])
pbs <- readRDS(args[[2]])

# analysis
sce <- logNormCounts(sce, BPPARAM=bp)
sel <- rownames(sce) %in% rownames(pbs)
sce <- runPCA(sce, 
    ncomponents=30, subset_row=sel,
    BSPARAM=RandomParam(), BPPARAM=bp)

# saving
rowData(sce)$sel <- sel
saveRDS(sce, args[[3]])