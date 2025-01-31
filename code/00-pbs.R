# dependencies
suppressPackageStartupMessages({
    library(scuttle)
    library(BiocParallel)
    library(SingleCellExperiment)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)

# loading
sce <- readRDS(args[[1]])
ids <- readRDS(args[[2]])

# wrangling
sce[[wcs$ids]] <- ids
sizeFactors(sce) <- NULL

# analysis
ids <- colData(sce)[c("sid", wcs$ids)]
pbs <- aggregateAcrossCells(sce, ids, BPPARAM=bp)
pbs <- logNormCounts(pbs, log=FALSE, BPPARAM=bp)

# saving
saveRDS(pbs, args[[3]])
