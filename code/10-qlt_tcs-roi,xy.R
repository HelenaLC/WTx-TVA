# dependencies
suppressPackageStartupMessages({
    library(sf)
    library(sp)
    library(arrow)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

# helper to get spatial coordinates
.xy <- \(.) {
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(.)))
    return(colData(.)[xy])
}

.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
args[1:4] <- lapply(args[1:4], `names<-`, sid <- .sid(args[[1]]))

ps <- lapply(sid, \(s) {
    # loading
    sce <- readRDS(args[[1]][s]); assays(sce) <- list()
    roi <- readRDS(args[[2]][s]); assays(roi) <- list()
    pol <- read_parquet(args[[3]][s], as_data_frame=FALSE)
    idx <- match(colnames(sce), colnames(roi))
    sce$typ <- gsub("^.*_", "", roi$typ)[idx]
    # get regions
    ids <- sort(setdiff(unique(roi$roi), NA))
    ids <- grep("(LI|BV|REF)$", ids, value=TRUE, invert=TRUE)
    lapply(ids, \(i) {
        roj <- roi[, grep(i, roi$roi)] # subset region
        ch <- concaveman(as.matrix(.xy(roj))) # concave hull
        hc <- st_buffer(st_polygon(list(ch)), 0.05) # expansion
        # subset cells
        xy <- .xy(sce); yx <- st_coordinates(hc)
        cs <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
        tmp <- sce[, colnames(sce)[cs != 0]]
        # plotting
        .plt_ps(pol, tmp, "typ", id=i, lw=0.05, a=2/3) +
            scale_fill_manual(values=.pal_roj, na.value="grey") +
            theme(legend.position="none") +
            geom_polygon(fill=NA, col="black", linewidth=0.2,
                aes(V1, V2), data.frame(ch), inherit.aes=FALSE)
    })
}) |> Reduce(f=base::c)

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- ps[[.]]$data
    dx <- diff(range(df$x))*10
    dy <- diff(range(df$y))*10
    pdf(tf[[.]], 
        width=(5+dx)/2.54, 
        height=(0.5+dy)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[5]])
