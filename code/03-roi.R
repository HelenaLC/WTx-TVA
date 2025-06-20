# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
names(roi) <- basename(roi <- args[[2]])

# split ROIs into ones denoting class(ification) 
# as REF/TVA/CRC and trans(itions), respectively
names(grp) <- grp <- ifelse(
    !grepl("ROI", roi), "typ",
    ifelse(grepl("REF$", roi), 
    "ref", "roi"))
roi <- split(roi, grp)
sapply(roi, length)

# assign ROI labels to cells
for (. in names(roi)) {
    if (length(roi[[.]])) {
        idx <- lapply(roi[[.]], \(.) {
            roi <- .align_shape(sce, .)
            sub <- .subset_shape(sce, roi)
            colnames(sub)
        }) 
        lab <- rep.int(names(idx), vapply(idx, length, numeric(1)))
        sce[[.]] <- lab[match(colnames(sce), unlist(idx))]
    } else {
        sce[[.]] <- NA
    }
}

# tabulate assignments
roi <- as.matrix(colData(sce)[names(roi)])
table(rowAlls(is.na(roi)))
apply(roi, 2, table)

# saving
base::saveRDS(sce, args[[3]])