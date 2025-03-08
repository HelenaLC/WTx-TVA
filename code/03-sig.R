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
set.seed(250308)

# loading
sce <- readRDS(args[[1]])
sig <- lapply(args[[2]], read.delim, header=FALSE)
sapply(sig <- lapply(sig, unlist), length)

# retrieve sets
url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=%s&fileType=json"
names(gs) <- gs <- unlist(sig)
gs <- lapply(gs, \(.) {
    url <- sprintf(url, .)
    tf <- tempfile(fileext=".json")
    download.file(url, tf, quiet=TRUE)
    tryCatch({
        gs <- fromJSON(tf)[[1]]$geneSymbols
        intersect(gs, rownames(sce))
    }, error=\(e) NULL)
})
gs <- gs[!vapply(gs, is.null, logical(1))]
range(sapply(gs, length)); length(gs)

# log-library size normalization
mtx <- normalizeCounts(sce, log=TRUE)

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
sub <- gsub(".*(imm|epi|str).*", "\\1", args[[2]])
sub <- rep.int(sub, sapply(sig, length))
names(sub) <- unlist(sig)
sub <- sub[names(gs)]
rowData(res)$set <- gs
rowData(res)$sub <- sub
saveRDS(res, args[[3]])
