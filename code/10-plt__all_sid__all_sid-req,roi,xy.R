# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# loading
ps <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        trj <- readRDS(x)
        sce <- readRDS(y)
        # restrict to cells in REF region
        sce <- sce[, !is.na(sce$ref)]
        idx <- intersect(colnames(sce), colnames(trj))
        sce <- sce[, idx]; trj <- trj[, idx]$slingshot
        # plotting
        id <- gsub(".*([0-9]{3}).*", "\\1", x)
        pt <- unname(.q(averagePseudotime(trj)))
        .plt_xy(sce, pt, id=id, split=FALSE) +
        theme(legend.position="none")
    })

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- ps[[.]]$data
    dx <- diff(range(df$x))/3
    dy <- diff(range(df$y))/3
    pdf(tf[[.]], 
        width=(3+dx)/2.54, 
        height=(0.5+dy)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
