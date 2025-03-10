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

df <- mapply(
    SIMPLIFY=FALSE,
    sce=args[[1]], roi=args[[2]], 
    pol=args[[3]], ist=args[[4]],
    \(sce, roi, pol, ist) {
        # loading
        sce <- readRDS(sce); assays(sce) <- list()
        roi <- readRDS(roi); assays(roi) <- list()
        pol <- read_parquet(pol, as_data_frame=FALSE)
        kid <- readRDS(ist)$clust
        # get regions corresponding to transition crypts
        ids <- sort(setdiff(unique(roi$roi), NA))
        ids <- grep("(BV|LI)$", ids, invert=TRUE, value=TRUE)
        df <- lapply(ids, \(id) {
            roj <- roi[, grep(id, roi$roi)] # subset region
            ch <- concaveman(as.matrix(.xy(roj))) # concave hull
            hc <- st_buffer(st_polygon(list(ch)), 0.05) # expansion
            # subset cells
            xy <- .xy(sce); yx <- st_coordinates(hc)
            cs <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
            cs <- colnames(sce)[cs != 0]
            # get clusters
            kid <- kid[match(cs, names(kid))]
            # get region type
            typ <- gsub(".*(REF|TVA|CRC).*", "\\1", roi[, cs]$typ)
            data.frame(row.names=NULL, sid=sce$sid[1], roi=id, typ, kid)
        }) |> do.call(what=rbind)
    }) |> do.call(what=rbind)

# plotting
ps <- by(df, df$roi, \(fd) 
    .plt_fq(fd, "typ", "kid", hc=FALSE) + 
    ggtitle(fd$roi[1]))
gg <- wrap_plots(ps, nrow=6) & theme(
    aspect.ratio=1,
    legend.position="none",
    panel.grid=element_blank(),
    axis.title=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    plot.margin=margin(2,0,0,0),
    plot.title=element_text(size=3),
    panel.border=element_rect(fill=NA)) &
    scale_fill_manual(values=.pal_sub)

# saving
ggsave(args[[5]], gg, units="cm", width=7, height=7.5)
