# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(BiocParallel)
    library(SingleCellExperiment)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)

# loading
sce <- readRDS(args[[1]])
ids <- readRDS(args[[2]])

# analysis
mgs <- findMarkers(sce, 
    groups=ids, direction="up",
    add.summary=TRUE, BPPARAM=bp)

# saving
saveRDS(mgs, args[[3]])
