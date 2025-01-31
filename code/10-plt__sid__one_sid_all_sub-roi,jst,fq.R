# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
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
    df <- data.frame(colData(sce)[idx, ], kid) |>
        filter(!is.na(typ), !is.na(kid)) |>
        mutate(typ=gsub(".*_", "", typ)) |>
        group_by(typ, kid) |>
        tally() |> 
        mutate(p=n/sum(n))
}) 

# plotting
ps <- lapply(sub, \(.) {
    nc <- sum(df[[.]]$n)
    nc <- format(nc, big.mark=",")
    yo <- gsub(".*(REF|TVA|CRC).*", "\\1", unique(df[[.]]$typ))
    yo <- order(factor(yo, c("REF", "TVA", "CRC")))
    ggplot(df[[.]], aes(typ, p, fill=kid)) + 
    ggtitle(bquote(bold(.(.))~"(N ="~.(nc)*")")) +
    guides(fill=guide_legend(ncol=1, override.aes=list(shape=21, stroke=0, size=2))) +
    geom_col(col="white", linewidth=0.1, width=1, key_glyph="point") +
    scale_x_discrete(NULL, limits=\(.) .[yo]) +
    scale_y_continuous(wcs$sid, n.breaks=2) +
    scale_fill_manual(NULL, values=.pal) +
    coord_fixed(10, expand=FALSE) +
    theme_bw(6) + theme(
        legend.key=element_blank(),
        axis.ticks.x=element_blank(),
        plot.background=element_blank(),
        legend.key.size=unit(0, "lines"),
        axis.title.y=element_text(face="bold"),
        plot.title=element_text(hjust=0.5, size=6),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        (if (. != "epi") theme(axis.title.y=element_blank()))
})

# saving
gg <- wrap_plots(ps, nrow=1)
ggsave(args[[3]], gg, units="cm", width=15, height=4)
