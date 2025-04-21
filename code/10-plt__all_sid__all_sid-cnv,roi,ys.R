# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(zellkonverter)
    library(SingleCellExperiment)
})

args <- lapply(args, \(.) .[!grepl("210", .)])

# loading
df <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        cnv <- readRDS(x)
        sce <- readRDS(y)
        # wrangling
        sce <- sce[, !is.na(sce$typ)]
        i <- intersect(colnames(sce), colnames(cnv))
        df <- data.frame(
            cnv=.z(colMeans(abs(assay(cnv))[, i])),
            roi=sce[, i]$typ, sid=sce$sid[1])
    }) |> do.call(what=rbind)

# plotting
typ <- gsub("^.*_", "", df$roi)
df$typ <- factor(typ, names(.pal_roj))
gg <- ggplot(df, aes(typ, cnv, fill=typ)) + geom_boxplot(
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
