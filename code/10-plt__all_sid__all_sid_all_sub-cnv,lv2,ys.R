# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidytext)
    library(patchwork)
    library(zellkonverter)
    library(SingleCellExperiment)
})

# loading
args[[2]] <- grep("epi", args[[2]], value=TRUE)
df <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        cnv <- readRDS(x)
        ist <- readRDS(y)$clust
        ist <- gsub("^epi\\.", "", ist)
        # wrangling
        sid <- gsub(".*([0-9]{3}).*", "\\1", x)
        idx <- intersect(colnames(cnv), names(ist))
        kid <- gsub("^epi\\.", "", ist[idx])
        cnv <- .z(colMeans(abs(assay(cnv)))[idx])
        data.frame(sid, kid, cnv)
    }) |> do.call(what=rbind)

# plotting
gg <- ggplot(df, aes(y=cnv, fill=kid,
    x=reorder_within(kid, cnv, sid, median))) + geom_boxplot(
    outlier.shape=16, outlier.size=0.2, outlier.stroke=0,
    key_glyph="point", show.legend=TRUE, alpha=2/3, linewidth=0.2) +
    facet_wrap(~sid, nrow=2, scales="free_x") +
    scale_fill_manual(NULL, values=.pal_kid) +
    geom_hline(yintercept=0, linewidth=0.2) +
    labs(y="z-scaled CNV score") +
    scale_x_reordered(NULL) +
    .thm_fig_d("bw", "f") + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())

# saving
ggsave(args[[3]], gg, unit="cm", width=12, height=5)
