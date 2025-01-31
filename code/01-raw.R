# dependencies
suppressPackageStartupMessages({
    library(Matrix)
    library(HDF5Array)
    library(SparseArray)
    library(SingleCellExperiment)
})

# setup
did <- basename(args[[1]])
md <- read.delim(args[[2]])
(md <- md[md$did == did, ])

# loading
f <- \(.) file.path(args[[1]], .)
cd <- read.csv(f("metadata.csv.gz"))
y <- readSparseCSV(f("counts.csv.gz"), transpose=TRUE)

# for 1st run, exclude targets
# that overlap with CARD8 barcode
if (md$run[1] == 1) {
    ex <- readRDS(args[[3]])
    y <- y[!rownames(y) %in% ex, ]
}

# add coordinates in mm
xy <- "^Center(X|Y)_(local|global)_px"
px <- grep(xy, names(cd), value=TRUE)
mm <- gsub("_px$", "_mm", px)
cd[mm] <- sapply(cd[px], .px2mm)

# coercion
y <- as(y[-1, ], "dgCMatrix"); colnames(y) <- cd$cell
sce <- SingleCellExperiment(list(counts=y), colData=cd)
names(ex) <- ex <- c("Negative", "SystemControl")
ex <- lapply(ex, grep, rownames(y))
altExps(sce) <- lapply(ex, \(.) sce[., ])
(sce <- sce[-unlist(ex), ])

# saving
for (. in seq(nrow(md))) {
  dm <- md[., setdiff(names(md), "fov")]
  colData(sce)[names(dm)] <- as.list(dm)
  fs <- md$fov[.]; cs <- TRUE
  if (fs != "") {
    ss <- lapply(sapply(strsplit(fs, ",")[[1]], strsplit, "-"), as.integer)
    cs <- sce$fov %in% unlist(lapply(ss, \(.) do.call(seq, as.list(range(.)))))
  }
  saveHDF5SummarizedExperiment(sce[, cs], args[[3+.]], replace=TRUE)
}