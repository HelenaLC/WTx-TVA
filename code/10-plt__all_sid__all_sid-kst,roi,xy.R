# setwd("~/projects/cosmx-wtx")
# source("code/_utils.R")
# args <- list(
#     list.files("outs", "kst", full.names=TRUE),
#     list.files("outs", "roi", full.names=TRUE)[-8],
#     "plts/kst,roi,xy.pdf")

# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
ps <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        ist <- readRDS(x)$clust
        sce <- readRDS(y)
        # restrict to cells in REF region
        sce <- sce[, !is.na(sce$ref)]
        idx <- match(colnames(sce), names(ist))
        # plot joint/split cluster assignments
        sid <- gsub(".*([0-9]{3}).*", "\\1", x)
        ps <- .plt_xy(sce, ist[idx], id=sid)
        ps[-1] <- lapply(ps[-1], \(p) {
            n <- sum(p$data$.) 
            s <- p$layers[[1]]$aes_params$size
            p$layers[[1]]$aes_params$size <- s*(1e3/n); p })
        ps
    }) |> Reduce(f=c)

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- ps[[.]]$data
    ps[[.]]$guides$guides$colour$params$override.aes$size <- 1
    dx <- diff(range(df$x))/3
    dy <- diff(range(df$y))/3
    pdf(tf[[.]], 
        width=(3+dx)/2.54, 
        height=(0.5+dy)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
