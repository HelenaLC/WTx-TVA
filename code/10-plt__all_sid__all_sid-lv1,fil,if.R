# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
ist <- lapply(args[[1]], readRDS)
sce <- lapply(args[[2]], readRDS)

# wrangling
df <- mapply(
    SIMPLIFY=FALSE, sce=sce,
    kid=lapply(ist, \(.) .$clust), 
    \(sce, kid) {
        idx <- match(colnames(sce), names(kid))
        data.frame(colData(sce), kid=kid[idx])
    }) |>
    do.call(what=rbind) |>
    filter(!is.na(kid)) |>
    mutate_at("sid", factor) 

ah <- \(., cf=150) asinh(./cf)
nm <- \(.) gsub("Mean\\.", "", .)
xs <- grep("^Mean\\.", names(df), value=TRUE)
fd <- df |>
    select(all_of(c("sid", "kid", xs))) |>
    rename_with(nm, all_of(xs)) |>
    pivot_longer(-c(sid, kid))
mu <- fd |>
    group_by(name, kid) |>
    summarise_at("value", median)

# plotting
gg <- ggplot(fd, aes(value, sid, fill=kid)) + 
    geom_boxplot(alpha=1/3, linewidth=0.1, outliers=FALSE) +
    scale_color_manual(values=c("gold", "cyan", "magenta")) +
    scale_fill_manual(values=c("gold", "cyan", "magenta")) +
    geom_vline(
        aes(xintercept=value, col=kid), mu,
        show.legend=FALSE, linewidth=0.4) +
    facet_wrap(~name, nrow=1, scales="free_x") +
    scale_y_discrete(limits=\(.) rev(.)) +
    scale_x_log10(breaks=10**seq(0, 6),
        labels=scales::label_log(digits=2)) +
    theme_bw(6) + theme(
        aspect.ratio=1,
        legend.position="none",
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        plot.background=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())

# saving
w <- 3*length(unique(fd$name))
ggsave(args[[3]], gg, units="cm", width=w, height=3)
