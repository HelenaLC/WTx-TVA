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
        df <- data.frame(sid, sub, sig=rownames(pbs), assay(pbs), check.names=FALSE)
        fd <- pivot_longer(df, all_of(colnames(pbs)), names_to="kid")
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

df <- df |>
    mutate(sig=gsub("HALLMARK_", "", sig)) |>
    mutate(sig=gsub("DESCARTES_FETAL_INTESTINE_", "", sig)) 

# plotting
.p <- \(x, i, n, xo=TRUE, yo=TRUE) {
    y <- pivot_wider(x, names_from="kid")
    z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
    ggplot(x, aes(kid, sig, fill=value)) +
        (if (xo) scale_x_discrete(limits=.yo(z))) +
        (if (yo) scale_y_discrete(limits=.xo(z))) +
        scale_fill_gradient2(
            "z-scaled\nmean AUCell",
            low="turquoise", high="purple",
            limits=c(-2.5, 2.5), n.breaks=5) +
        ggtitle(bquote(bold(.(i))~"(N ="~.(n)*")")) +
        geom_tile() + coord_equal(expand=FALSE) +
        .thm_fig_c("bw") + theme(
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
}

# joint
fd <- df |>
    # average across samples
    group_by(sig, sub, kid) |>
    summarise_at("value", mean) |>
    # scale across clusters
    group_by(sig, sub) |>
    mutate(value=.z(value)) |>
    ungroup()
p0 <- lapply(sub, \(.) {
    x <- select(fd[fd$sub == ., ], -sub)
    .p(x, ., length(unique(x$kid))) + 
        if (. != "epi") theme(axis.text.y=element_blank())
}) |> wrap_plots(nrow=1) + plot_layout(guides="collect")

# split by sample
fd <- df |>
    # average by sample
    group_by(sid, sig, sub, kid) |>
    summarise_at("value", mean) |>
    # scale across clusters
    group_by(sid, sig) |>
    mutate(value=.z(value)) |>
    ungroup()
ps <- lapply(sid, \(sid) { 
    x <- select(fd[fd$sid == sid, ], -c(sid, sub))
    .p(x, sid, length(unique(x$kid)), yo=TRUE)
})

# split by subset
fd <- df |>
    # average by sample
    group_by(sid, sig, sub, kid) |>
    summarise_at("value", mean) |>
    # scale across clusters
    group_by(sid, sig, sub) |>
    mutate(value=.z(value)) |>
    ungroup()
qs <- lapply(sid, \(sid) { 
    lapply(sub, \(sub) {
        x <- select(fd[fd$sid == sid & fd$sub == sub, ], -c(sid, sub))
        .p(x, paste0(sid, ": ", sub), length(unique(x$kid))) +
            if (sub != "epi") theme(axis.text.y=element_blank())
    }) |> wrap_plots(nrow=1) + plot_layout(guides="collect")
})

# saving
pdf(args[[3]], width=15/2.54, height=12/2.54, onefile=TRUE)
for (p in c(list(p0), ps, qs)) print(p); dev.off()
