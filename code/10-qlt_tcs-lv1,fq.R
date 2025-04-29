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
        typ <- gsub("^.*_", "", roi$typ)
        roi$typ <- factor(typ, names(.pal_roj))
        ids <- sort(setdiff(unique(roi$roi), NA))
        ids <- grep("(BV|LI|REF)$", ids, invert=TRUE, value=TRUE)
        df <- lapply(ids, \(id) {
            roj <- roi[, grep(id, roi$roi)] # subset region
            ch <- concaveman(as.matrix(.xy(roj))) # concave hull
            hc <- st_buffer(st_polygon(list(ch)), 0.05) # expansion
            # subset cells
            xy <- .xy(sce)
            yx <- st_coordinates(hc)
            cs <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
            cs <- colnames(sce)[cs != 0]
            # in selected (w/o expansion)?
            yx <- .xy(roj)
            ct <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
            ct <- colnames(sce)[ct != 0]
            sel <- cs %in% ct
            # get clusters & region type
            kid <- kid[match(cs, names(kid))]
            typ <- roi[, cs]$typ
            data.frame(
                row.names=NULL, 
                sid=sce$sid[1], 
                roi=id, typ, kid, sel)
        }) |> do.call(what=rbind)
    }) |> do.call(what=rbind)

# plotting
fd <- df[df$sel & df$kid == "epi", ]
p1 <- .plt_fq(fd, "roi", "typ", id="epi") +
    scale_fill_manual(values=.pal_roj)

p2 <- .plt_fq(df, "roi", "kid", id="all") +
    scale_fill_manual(values=.pal_sub) +
    scale_x_discrete(limits=p1$scales$scales[[1]]$limits)

# aesthetics
ps <- lapply(list(p1, p2), 
    \(p) p + .theme_fig + theme(
        legend.position="none",
        axis.title=element_blank(),
        axis.text.y=element_blank()))

# saving
pdf(args[[5]], onefile=TRUE, width=6/2.54, height=4/2.54)
for (p in ps) print(p); dev.off()
