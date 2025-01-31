# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- lapply(args[[2]], readRDS)

# wrangling
sub <- gsub(".*,(.*)\\..*", "\\1", basename(args[[2]]))
names(ist) <- names(sub) <- sub
df <- lapply(sub, \(.) {
    idx <- names(kid <- ist[[.]]$clust)
    data.frame(colData(sce)[idx, ], kid) |>
        filter(!is.na(roi)) |>
        group_by(roi, kid) |>
        tally() |> 
        mutate(p=n/sum(n))
}) 

# plotting
ps <- lapply(sub, \(.) {
    nc <- sum(df[[.]]$n)
    nc <- format(nc, big.mark=",")
    id <- paste0(wcs$sid, ": ", .)
    ggplot(df[[.]], aes(p, roi, fill=kid)) + 
    ggtitle(bquote(bold(.(id))~"(N ="~.(nc)*")")) +
    guides(fill=guide_legend(ncol=1, override.aes=list(shape=21, stroke=0, size=2))) +
    geom_col(col="white", linewidth=0.1, width=1, key_glyph="point") +
    scale_fill_manual(NULL, values=.pal) +
    scale_x_continuous(n.breaks=6) +
    coord_cartesian(expand=FALSE) +
    theme_bw(6) + theme(
        aspect.ratio=2/3,
        legend.margin=margin(),
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.size=unit(0, "lines"),
        plot.title=element_text(hjust=0.5))
})

# saving
gg <- wrap_plots(ps, nrow=1)
ggsave(args[[3]], gg, units="cm", width=15, height=4)
