# dependencies
suppressPackageStartupMessages({
    library(scuttle)
    library(HDF5Array)
    library(BiocParallel)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# analysis
idx <- match(colnames(sce), names(ids <- ist$clust))
counts(sce) <- sweep(counts(sce), 2, sce$Area.um2, `/`)
pbs <- aggregateAcrossCells(sce, ids[idx], statistics="mean", BPPARAM=bp)

# saving
sel <- rownames(ist$profiles)
sel <- rownames(pbs) %in% sel
rowData(pbs)$sel <- sel
saveRDS(pbs, args[[3]])

