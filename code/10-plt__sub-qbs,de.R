# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(scater)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
sub <- sce[, sce$sid == "all"]
ks <- colnames(sub) <- paste(sub$kid)

# utils
.gs <- \(x, n, k=NULL) {
    es <- logcounts(x)
    ks <- colnames(x)
    gs <- lapply(ks, \(i) {
        j <- setdiff(ks, i)
        fc <- rowMeans(es[, i]/es[, j])
        names(tail(sort(fc), n))
    }) |> setNames(ks)
    if (!is.null(k)) 
        return(gs[[k]])
    unique(unlist(gs))
}
.df <- \(x) {
    data.frame(
        t(logcounts(x)), 
        k=colnames(x), 
        check.names=FALSE) |>
        pivot_longer(-k, names_to="g") |>
        group_by(g) |> mutate_at("value", .z)
}
.hm <- \(x, nm, nk, xo, yo) {
    ggplot(x, aes(g, k, fill=value)) +
        scale_x_discrete(limits=xo) +
        scale_y_discrete(limits=rev(yo)) +
        scale_fill_gradientn(
            "z-scaled\nmean expr.",
            colors=hcl.colors(11, "Blue-Red 3"),
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
        ggtitle(.lab(nm, nk)) +
        geom_tile() +
        coord_equal(4/3, expand=FALSE) +
        .theme_fig + theme(
            legend.position="none",
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_text(size=4),
            axis.text.x=element_text(size=3, angle=90, vjust=0.5, hjust=1))
}

# joint
df <- .df(sub[gs <- .gs(sub, 10), ])
p0 <- .hm(df, nm=wcs$sub, nk=length(ks), xo=gs, yo=ks)

# split
qs <- lapply(ks, \(k) {
    n <- colData(sub)[k, "ncells"]
    gs <- .gs(sub, 100, k)
    df <- .df(sub[gs, ])
    y <- assay(sub[gs, ])
    .hm(df, nm=k, nk=n, xo=.xo(y), yo=.yo(y))
})

# saving
pdf(args[[2]], width=12/2.54, height=4/2.54, onefile=TRUE)
for (p in c(list(p0), qs)) print(p); dev.off()
