# dependencies
suppressPackageStartupMessages({
    library(AUCell)
    library(msigdbr)
    library(scuttle)
    library(HDF5Array)
    library(BiocParallel)
    library(SingleCellExperiment)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)
set.seed(250131)

# loading
sce <- readRDS(args[[1]])
df <- msigdbr()

# retrieve sets
pat <- c(
    "HALLMARK",
    "DESCARTES_FETAL_INTESTINE")
pat <- paste(pat, collapse="|")
fd <- df[grepl(pat, df$gs_name), ]

# get gene symbols by set
fd <- fd[fd$gene_symbol %in% rownames(sce), ]
gs <- split(fd$gene_symbol, fd$gs_name)
range(sapply(gs, length)); length(gs)

# normalize by area
mtx <- sweep(counts(sce), 2, sce$Area.um2, `/`)

# gene set scoring
idx <- split(seq(ncol(sce)), sce$fov)
res <- bplapply(idx, BPPARAM=bp, \(.) {
    mtx <- as(mtx[, ., drop=FALSE], "dgCMatrix")
    rnk <- AUCell_buildRankings(mtx, plotStats=FALSE, verbose=FALSE)
    auc <- AUCell_calcAUC(gs, rnk, verbose=FALSE)
})
res <- do.call(cbind, res)
res <- res[, colnames(sce)]
colData(res) <- colData(sce)

# saving
rowData(res)$set <- gs
saveRDS(res, args[[2]])
