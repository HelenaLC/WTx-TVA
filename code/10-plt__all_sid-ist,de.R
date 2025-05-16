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
ist <- lapply(args[[1]], readRDS)

# wrangling
se <- lapply(ist, \(ist) {
    es <- normalizeCounts(ist$profiles)
    cd <- data.frame(kid=colnames(es))
    SingleCellExperiment(list(counts=es), colData=cd)
})

# utils
.gs <- \(x, n) {
    y <- assay(x)
    names(k) <- k <- colnames(x)
    lapply(k, \(i) {
        j <- setdiff(k, i)
        fc <- y[, i]/rowMeans(y[, j])
        names(tail(sort(fc), n))
    })
}
.df <- \(x) {
    data.frame(
        t(assay(x)), 
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
        ggtitle(bquote(bold(.(nm))~"(N ="~.(nk)*")")) +
        geom_tile() +
        coord_equal(4/3, expand=FALSE) +
        theme_bw(6) + theme(
            legend.position="bottom",
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            panel.grid=element_blank(),
            legend.key=element_blank(),
            plot.background=element_blank(),
            legend.background=element_blank(),
            legend.key.size=unit(0.4, "lines"),
            legend.key.width=unit(0.8, "lines"),
            legend.title=element_text(vjust=1),
            axis.text.y=element_text(size=4),
            axis.text.x=element_text(size=3, angle=90, vjust=0.5, hjust=1))
}

# aggregation
gs <- Reduce(intersect, lapply(se, rownames))
x <- do.call(cbind, lapply(se, \(.) .[gs, ]))
y <- aggregateAcrossCells(x, x$kid, statistics="mean")

# joint
gs <- unique(unlist(.gs(y, 10)))
p0 <- .hm(.df(y[gs, ]), 
    nm="all", nk=ncol(y), 
    xo=gs, yo=colnames(y))

# split
ps <- lapply(colnames(y), \(k) {
    gs <- .gs(y, 100)[[k]]
    ns <- sapply(ist, \(.) sum(.$clust == k))
    .hm(.df(y[gs, ]), 
        nm=k, nk=sum(ns), 
        xo=gs, yo=colnames(y))
})

# saving
pdf(args[[2]], width=15/2.54, height=5/2.54, onefile=TRUE)
for (p in c(list(p0), ps)) print(p); dev.off()
