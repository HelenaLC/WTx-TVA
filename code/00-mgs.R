# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scuttle)
    library(HDF5Array)
    library(BiocParallel)
    library(SingleCellExperiment)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)

# loading
lys <- mapply(
    SIMPLIFY=FALSE,
    sce=args[[1]], 
    ist=args[[2]], 
    \(sce, ist) {
        sce <- readRDS(sce)
        ist <- readRDS(ist)
        idx <- names(ids <- ist$clust)
        idx <- match(colnames(sce), idx)
        sce$kid <- ids[idx]
        metadata(sce) <- list()
        altExps(sce) <- list()
        return(sce)
    })

# wrangling
gs <- lapply(lys, rownames)
gs <- Reduce(intersect, gs)
lys <- lapply(lys, \(sce) sce[gs, ])
dim(sce <- do.call(cbind, lys))

# log-library size normalization
sce <- logNormCounts(sce, BPPARAM=bp)

# analysis
mgs <- findMarkers(sce, 
    block=sce$sid,
    groups=sce$kid, 
    direction="up",
    full.stats=TRUE,
    add.summary=TRUE, 
    BPPARAM=bp)

# saving
saveRDS(mgs, args[[3]])
