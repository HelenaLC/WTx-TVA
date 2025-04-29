# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(zellkonverter)
    library(SingleCellExperiment)
})

# loading
args <- lapply(args, \(.) .[!grepl("210", .)])
df <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        cnv <- readRDS(x)
        sce <- readRDS(y)
        sid <- sce$sid[1]
        # wrangling
        sce <- sce[, !is.na(sce$typ)]
        idx <- intersect(colnames(sce), colnames(cnv))
        cnv <- .z(colMeans(abs(assay(cnv[, idx]))))
        data.frame(cnv, colData(sce[, idx]))
    }) |> do.call(what=rbind)

# plotting
roi <- gsub("^.*_", "", df$typ)
df$roi <- factor(roi, names(.pal_roj))
gg <- ggplot(df, aes(roi, cnv, fill=roi)) + geom_boxplot(
    outlier.shape=16, outlier.size=0.2, outlier.stroke=0,
    key_glyph="point", show.legend=TRUE, alpha=2/3, linewidth=0.2) +
    facet_grid(~sid, space="free", scales="free_x") +
    scale_fill_manual("region", values=.pal_roj) +
    geom_hline(yintercept=0, linewidth=0.2) +
    labs(x=NULL, y="z-scaled CNV score") +
    .thm_fig_d("bw", "f") + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=10, height=4)
