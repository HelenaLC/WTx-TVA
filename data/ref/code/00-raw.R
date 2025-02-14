# dependencies
suppressPackageStartupMessages({
    library(DropletUtils)
    library(SingleCellExperiment)
})

# loading
md <- read.csv(args[[1]])
md <- md[md$sid == wcs$sid, ]
dir <- file.path("data/raw", md$library, md$barcode)
dir <- file.path(dir, "count/sample_filtered_feature_bc_matrix")
sce <- read10xCounts(dir)

# wrangling
cd <- md[i <- match(sce$Sample, unique(sce$Sample)), ]
rownames(cd) <- paste(sep="_", "c", 
    with(cd, gsub("[a-z]+", "", library)), 
    with(cd, gsub("[a-z]+", "", barcode)),
    unlist(sapply(table(i), seq)))
colData(sce) <- DataFrame(cd)

metadata(sce) <- list()
rowData(sce) <- rowData(sce)[-3]
names(rowData(sce)) <- c("gene_id", "gene_symbol")

gs <- rowData(sce)$gene_symbol
rownames(sce) <- make.unique(gs)

# saving
saveRDS(sce, args[[2]])