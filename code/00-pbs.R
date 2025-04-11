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

# joint
pbs <- aggregateAcrossCells(sce, 
    colData(sce)[ids <- "kid"],
    use.assay.type="logcounts", 
    statistics="mean", 
    BPPARAM=bp)
colData(pbs) <- colData(pbs)[c(ids, "ncells")]

# split
qbs <- aggregateAcrossCells(sce, 
    colData(sce)[ids <- c("sid", "kid")],
    use.assay.type="logcounts",
    statistics="mean", 
    BPPARAM=bp)
colData(qbs) <- colData(qbs)[c(ids, "ncells")]

# saving
pbs$sid <- "all"
out <- cbind(pbs, qbs)
ist <- readRDS(args[[2]][1])
sel <- gs %in% rownames(ist$profiles)
rowData(out)$sel <- sel
saveRDS(out, args[[3]])
