# loading
sce <- readRDS(args[[1]])
names(roi) <- basename(roi <- args[[2]])

# split ROIs into ones denoting class(ification) 
# as REF/TVA/CRC and trans(itions), respectively
grp <- factor(grepl("ROI", roi), c(FALSE, TRUE))
roi <- split(roi, grp); names(roi) <- c("typ", "roi")

# assign ROI labels to cells
for (. in names(roi)) if (length(roi[[.]])) {
    idx <- lapply(roi[[.]], \(.) {
        roi <- .align_shape(sce, .)
        sub <- .subset_shape(sce, roi)
        colnames(sub)
    }) 
    lab <- rep.int(names(idx), vapply(idx, length, numeric(1)))
    sce[[.]] <- lab[match(colnames(sce), unlist(idx))]
}

# tabulate assignments
roi <- as.matrix(colData(sce)[names(roi)])
table(rowAlls(is.na(roi)))
apply(roi, 2, table)

# saving
saveRDS(sce, args[[3]])