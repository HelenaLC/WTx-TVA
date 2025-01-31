# dependencies
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# splitting
i <- colnames(sce); j <- names(k <- ist$clust)
idx <- split(seq(ncol(sce)), k[match(i, j)])
sapply(sub <- lapply(idx, \(.) sce[, .]), ncol)

# saving
for (. in seq_along(sub)) {
    obj <- sub[[.]]
    rds <- args[[2+.]]
    saveRDS(obj, rds)
}