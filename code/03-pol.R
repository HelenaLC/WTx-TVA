# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(arrow)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
xy <- read.csv(gzfile(args[[2]]))
yx <- filter(xy, cell %in% colnames(sce))

# wrangling
cd <- colData(sce)
i <- match(yx$cell, sce$cell)
j <- setdiff(names(cd), names(yx))
at <- arrow_table(cbind(yx, cd[i, j]))

# saving
write_parquet(at, args[[3]])
