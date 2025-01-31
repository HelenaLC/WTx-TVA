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
sid <- gsub(".*-(.*),.*\\.rds", "\\1", basename(args[[1]]))
sub <- gsub(".*-.*,(.*)\\.rds", "\\1", basename(args[[1]]))
# 
se <- mapply(
    ist=ist, sid=sid, sub=sub,
    SIMPLIFY=FALSE, \(ist, sid, sub) {
        es <- normalizeCounts(ist$profiles)
        cd <- data.frame(sid, sub, kid=colnames(es))
        SingleCellExperiment(list(counts=es), colData=cd)
    })

# utils
.gs <- \(x) {
    y <- assay(x)
    gs <- lapply(colnames(x), \(i) {
        j <- setdiff(colnames(x), i)
        fc <- y[, i]/rowMeans(y[, j])
        names(tail(sort(fc), 10))
    }) 
    unique(unlist(gs))
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

# joint
ps <- lapply(split(seq_along(se), sub), \(i) {
    gs <- Reduce(intersect, lapply(se[i], rownames))
    x <- do.call(cbind, lapply(se[i], \(.) .[gs, ]))
    y <- aggregateAcrossCells(x, x$kid, statistics="mean")
    .hm(.df(y[gs <- .gs(y), ]), 
        nm=x$sub[1], nk=ncol(y), 
        xo=gs, yo=colnames(y))
})

# split
qs <- lapply(split(seq_along(se), sid), \(i) {
    lapply(split(seq_along(se[i]), sub[i]), \(j) {
        df <- .df((x <- se[i][[j]])[gs <- .gs(x), ])
        nm <- paste0(x$sid[1], ": ", x$sub[1])
        .hm(df, nm, nk=ncol(x), xo=gs, yo=colnames(x))
    })
}) |> Reduce(f=c)

# saving
pdf(args[[2]], width=15/2.54, height=5/2.54, onefile=TRUE)
for (p in c(ps, qs)) print(p); dev.off()
