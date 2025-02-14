# dependencies
suppressPackageStartupMessages({
    library(sf)
    library(sp)
    library(arrow)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

# # loading
# args <- list(
#     list.files("outs", "fil", full.names=TRUE),
#     list.files("outs", "roi", full.names=TRUE),
#     list.files("outs", "pol", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE))

# helper to get spatial coordinates
.xy <- \(.) {
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(.)))
    return(colData(.)[xy])
}

.sid <- \(.) split(., gsub(".*([0-9]{3}).*", "\\1", .))
.sub <- \(.) split(., gsub(".*(epi|imm|str).*", "\\1", .))

# get clusters per subset across all samples
kid <- lapply(.sub(args[[4]]), \(.) lapply(.sid(.), \(.) readRDS(.)$clust))
ks <- lapply(kid, \(.) sort(unique(unlist(.))))

sub <- names(kid); sid <- names(kid[[1]])
args[1:3] <- lapply(args[1:3], `names<-`, sid)
df <- lapply(sid, \(sid) {
    # loading
    sce <- readRDS(args[[1]][sid]); assays(sce) <- list()
    roi <- readRDS(args[[2]][sid]); assays(roi) <- list()
    pol <- read_parquet(args[[3]][sid], as_data_frame=FALSE)
    # get regions, except lymphatic invasions
    ids <- sort(setdiff(unique(roi$roi), NA))
    ids <- grep("LI$", ids, invert=TRUE, value=TRUE)
    lapply(ids, \(id) {
        roj <- roi[, grep(id, roi$roi)] # subset region
        ch <- concaveman(as.matrix(.xy(roj))) # concave hull
        hc <- st_buffer(st_polygon(list(ch)), 0.05) # expansion
        # subset cells
        xy <- .xy(sce); yx <- st_coordinates(hc)
        cs <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
        cs <- colnames(sce)[cs != 0]
        lapply(sub, \(sub) {
            # get clusters
            kid <- kid[[sub]][[sid]]
            kid <- kid[match(cs, names(kid))]
            data.frame(sid, roi=id, sub, kid)
        }) |> do.call(what=rbind)
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# plotting
ps <- by(df, df$sub, \(fd) {
    id <- fd$sub[1]
    nc <- sum(!is.na(fd$kid))
    .plt_fq(fd, "roi", "kid") +
    ggtitle(.lab(id, nc))
})
xo <- ps[[1]]$scales$scales[[2]]$limits
ps[-1] <- lapply(ps[-1], `+`, scale_x_discrete(limits=xo))
gg <- wrap_plots(ps, ncol=1) +
    plot_layout(guides="collect") &
    theme(
        axis.title=element_blank(),
        legend.justification="left",
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=4))

# saving
ggsave(args[[5]], gg, units="cm", width=10, height=12)
