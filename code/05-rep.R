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
set.seed(250201)

# loading
sce <- readRDS(args[[1]])
pbs <- readRDS(args[[2]])

# analysis
logcounts(sce) <- sweep(counts(sce), 2, sce$Area.um2, `/`)
sel <- rownames(sce) %in% rownames(pbs)
sce <- runPCA(sce, 
    ncomponents=30, subset_row=sel,
    BSPARAM=RandomParam(), BPPARAM=bp)

# saving
rowData(sce)$sel <- sel
base::saveRDS(sce, args[[3]])
