# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(AUCell)
    library(scuttle)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)

sid <- sort(unique(.sid(args[[1]])))
sub <- sort(unique(.sub(args[[2]])))

df <- lapply(sub, \(sub) {
    lapply(sid, \(sid) {
        sce <- sce[[which(.sid(args[[1]]) == sid)]]
        ist <- ist[[which(.sid(args[[2]]) == sid & .sub(args[[2]]) == sub)]]
        idx <- match(colnames(sce), names(kid <- ist$clust))
        pbs <- aggregateAcrossCells(sce, kid[idx], use.assay.type="AUC", statistics="mean")
        df <- data.frame(sid, sub, sig=rownames(pbs), assay(pbs))
        fd <- pivot_longer(df, all_of(colnames(pbs)), names_to="kid")
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# plotting
.p <- \(x, i, n, xo=TRUE, yo=TRUE) {
    y <- pivot_wider(x, names_from="kid")
    z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
    nm <- paste0(sid, ": ", sub); nk <- ncol(z)
    ggplot(x, aes(kid, sig, fill=value)) +
        (if (xo) scale_x_discrete(limits=.yo(z))) +
        (if (yo) scale_y_discrete(limits=.xo(z))) +
        scale_fill_gradient2("z-scaled\nmean expr.",
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2),
            low="royalblue3", mid="ivory", high="tomato3") +
        ggtitle(bquote(bold(.(i))~"(N ="~.(n)*")")) +
        geom_tile() +
        coord_equal(expand=FALSE) +
        theme_bw(6) + theme(
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            panel.grid=element_blank(),
            legend.key.size=unit(0.4, "lines"),
            legend.title=element_text(vjust=1),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
}

fd <- df |>
    # average across samples
    group_by(sig, sub, kid) |>
    summarise_at("value", mean) |>
    # scale across clusters
    group_by(sig, sub) |>
    mutate(value=.z(value)) |>
    ungroup()
p0 <- lapply(sub, \(sub) {
    x <- select(fd[fd$sub == sub, ], -sub)
    .p(x, sub, length(unique(x$kid))) +
    if (sub != "epi") theme(axis.text.y=element_blank())
}) |> wrap_plots(nrow=1) + plot_layout(guides="collect")

# split by subset
fd <- df |>
    # average by sample
    group_by(sid, sig, sub, kid) |>
    summarise_at("value", mean) |>
    # scale across clusters
    group_by(sid, sig, sub) |>
    mutate(value=.z(value)) |>
    ungroup()
ps <- lapply(sid, \(sid) { 
    lapply(sub, \(sub) {
        x <- select(fd[fd$sid == sid & fd$sub == sub, ], -c(sid, sub))
        .p(x, paste0(sid, ": ", sub), length(unique(x$kid))) +
            if (sub != "epi") theme(axis.text.y=element_blank())
    }) |> wrap_plots(nrow=1) + plot_layout(guides="collect")
})

# split by sample
fd <- df |>
    # average by sample
    group_by(sid, sig, sub, kid) |>
    summarise_at("value", mean) |>
    # scale across clusters
    group_by(sid, sig) |>
    mutate(value=.z(value)) |>
    ungroup()
qs <- lapply(sid, \(sid) { 
    x <- select(fd[fd$sid == sid, ], -c(sid, sub))
    .p(x, sid, length(unique(x$kid)), yo=TRUE)
})

# saving
pdf(args[[3]], width=15/2.54, height=10/2.54, onefile=TRUE)
for (p in c(list(p0), ps, qs)) print(p); dev.off()
