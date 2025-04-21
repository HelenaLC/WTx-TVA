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

.sid <- \(.) split(., gsub(".*([0-9]{3}).*", "\\1", .))
.sub <- \(.) split(., gsub(".*(epi|imm|str).*", "\\1", .))

# get clusters per subset across all samples
kid <- lapply(.sub(args[[4]]), \(.) lapply(.sid(.), \(.) readRDS(.)$clust))
ks <- lapply(kid, \(.) sort(unique(unlist(.))))

sub <- names(kid); sid <- names(kid[[1]])
args[1:3] <- lapply(args[1:3], `names<-`, sid)
df <- lapply(sid, \(s) {
    # loading
    sce <- readRDS(args[[1]][s]); assays(sce) <- list()
    roi <- readRDS(args[[2]][s]); assays(roi) <- list()
    pol <- read_parquet(args[[3]][s], as_data_frame=FALSE)
    # get regions, except lymphatic invasions
    ids <- sort(setdiff(unique(roi$roi), NA))
    ids <- ids[!grepl("REF$", ids)]
    lapply(ids, \(i) {
        roj <- roi[, grep(i, roi$roi)] # subset region
        ch <- concaveman(as.matrix(.xy(roj))) # concave hull
        hc <- st_buffer(st_polygon(list(ch)), 0.05) # expansion
        # subset cells
        xy <- .xy(sce); yx <- st_coordinates(hc)
        cs <- point.in.polygon(xy[, 1], xy[, 2], yx[, 1], yx[, 2])
        cs <- colnames(sce)[cs != 0]
        lapply(sub, \(sub) {
            # get clusters
            kid <- kid[[sub]][[s]]
            kid <- kid[match(cs, names(kid))]
            data.frame(sid=s, roi=i, sub, kid)
        }) |> do.call(what=rbind)
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# plotting
typ <- grepl("(LI|BV)$", df$roi)
df$typ <- c("crypt", "vessel")[typ+1]
ps <- by(df, df$typ, \(df) {
    ps <- by(df, df$sub, \(fd) {
        id <- fd$sub[1]
        nc <- sum(!is.na(fd$kid))
        fd$kid <- factor(fd$kid, ks[[id]])
        # split by crypts vs. vessels
        gg <- .plt_fq(fd, "roi", "kid", h=TRUE) +
            labs(fill=id) + ggtitle(.lab(id, nc)) 
        gg$guides$guides$fill$params$override.aes$size <- 1
        gg + theme(aspect.ratio=length(unique(fd$roi))/10) +
        (if (id != "imm") theme(axis.title.x=element_blank())) +
        (if (id == "epi") theme(axis.text.y=element_text(size=4, hjust=1)))
    })
    xo <- ps[[1]]$scales$scales[[2]]$limits
    ps[-1] <- lapply(ps[-1], `+`, scale_x_discrete(limits=xo))
    wrap_plots(ps, nrow=1) +
        plot_layout(guides="collect") &
        theme(
            legend.justification="left",
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            legend.spacing=unit(-0.2, "lines"))
})

# saving
pdf(args[[5]], onefile=TRUE, width=10/2.54, height=10/2.54)
for (p in ps) print(p); dev.off()
