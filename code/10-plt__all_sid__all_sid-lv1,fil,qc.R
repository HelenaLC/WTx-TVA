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
args[[3]] <- "plts/lv1,fil,qc.pdf"
# wrangling
df <- mapply(
    SIMPLIFY=FALSE, sub=sce, 
    kid=lapply(ist, \(lys) lys$clust), 
    \(sub, kid) {
        idx <- match(colnames(sub), names(kid))
        data.frame(colData(sub), kid=kid[idx])
    }) |> 
    do.call(what=rbind) |>
    filter(!is.na(kid)) |>
    mutate_at("sid", factor) 

fd <- df |>
    mutate(
        `nCount_RNA/um2`=nCount_RNA/Area.um2,
        `nFeature_RNA/um2`=nFeature_RNA/Area.um2) |>
    pivot_longer(matches("^n.*_RNA")) |>
    mutate(
        name=gsub("nCount_RNA", "counts", name),
        name=gsub("nFeature_RNA", "features", name),
        name=factor(name, c(
            "counts", "features",
            "counts/um2", "features/um2")))
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
