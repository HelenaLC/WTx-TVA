# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(tidytext)
    library(SingleCellExperiment)
})

# loading
df <- lapply(args[[1]], \(x) {
    cnv <- readRDS(x)
    sco <- .z(colMeans(assay(cnv)))
    sid <- gsub(".*([0-9]{3}).*", "\\1", x)
    data.frame(sid, clu=cnv$clu, cnv=sco)
}) |> do.call(what=rbind)

# plotting
gg <- ggplot(df, aes(y=cnv, fill=clu,
    x=reorder_within(clu, cnv, sid, median))) + geom_boxplot(
    outlier.shape=16, outlier.size=0.2, outlier.stroke=0,
    key_glyph="point", show.legend=TRUE, alpha=2/3, linewidth=0.2) +
    scale_fill_manual("cnv", values=pals::kelly()) +
    facet_wrap(~sid, nrow=2, scales="free_x") +
    geom_hline(yintercept=0, linewidth=0.2) +
    labs(y="z-scaled CNV score") +
    scale_x_reordered(NULL) +
    .thm_fig_d("bw", "f") + theme(
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(size=4, angle=90, hjust=1, vjust=0.5))

# saving
gg$guides$guides$fill$params$override.aes$size <- 1.5
ggsave(args[[2]], gg, unit="cm", width=12, height=6)
