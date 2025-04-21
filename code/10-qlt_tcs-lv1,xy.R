# dependencies
suppressPackageStartupMessages({
    library(sf)
    library(sp)
    library(arrow)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

# args <- list(
#     list.files("outs", "fil", full.names=TRUE),
#     list.files("outs", "roi", full.names=TRUE),
#     list.files("outs", "pol", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/tcs,lv2,xy.pdf")
# args <- lapply(args, \(.) grep("241", ., value=TRUE, invert=TRUE))

# helper to get spatial coordinates
.xy <- \(.) {
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(.)))
    return(colData(.)[xy])
}

.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
args[1:4] <- lapply(args[1:4], `names<-`, sid <- .sid(args[[1]]))

ps <- lapply(sid, \(sid) {
    # loading
    sce <- readRDS(args[[1]][sid]); assays(sce) <- list()
    roi <- readRDS(args[[2]][sid]); assays(roi) <- list()
    pol <- read_parquet(args[[3]][sid], as_data_frame=FALSE)
    ist <- readRDS(args[[4]][sid])$clust
    idx <- match(colnames(sce), names(ist))
    sce$ist <- factor(ist[idx])
    # get regions
    ids <- sort(setdiff(unique(roi$roi), NA))
    ids <- ids[!grepl("REF$", ids)]
    lapply(ids, \(id) {
        roj <- roi[, grep(id, roi$roi)] # subset region
        ch <- concaveman(as.matrix(.xy(roj))) # concave hull
        hc <- st_buffer(st_polygon(list(ch)), 0.05) # expansion
        # subset cells
        xy <- .xy(sce); yx <- st_coordinates(hc)
        cs <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
        tmp <- sce[, colnames(sce)[cs != 0]]
        # plotting
        .plt_ps(pol, tmp, "ist", id=id, lw=0.05, a=2/3) +
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
