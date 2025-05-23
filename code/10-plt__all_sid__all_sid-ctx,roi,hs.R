# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# Shannon entropy
.h <- \(.) { . <- .[. > 0]; -sum(.*log(.)) }

df <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        df <- readRDS(x)
        se <- readRDS(y)
        # consider only cells in HP regions
        cs <- colnames(se)[!is.na(se$typ)]
        df <- filter(df, cid %in% cs)
        df$typ <- se$typ[match(df$cid, cs)]
        mx <- as.matrix(select(df, where(is.numeric)))
        data.frame(df, h=apply(mx, 1, .h))
    }) |> do.call(what=rbind)

# wrangling
gg <- ggplot(df, aes(
    reorder(typ, h, median, na.rm=TRUE), 
    h, fill=gsub(".*_", "", typ))) + 
    geom_boxplot(
        key_glyph="point", alpha=2/3, linewidth=0.2,
        outlier.shape=16, outlier.size=0.2, outlier.stroke=0) +
    scale_y_continuous("entropy", limits=c(0, 3)) +
    scale_fill_manual("region", values=.pal_roj) +
    geom_hline(yintercept=0, linewidth=0.2) +
    .thm_fig_d("bw", "f") + theme(
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=3, angle=45, hjust=1, vjust=1))

# saving
ggsave(args[[3]], gg, units="cm", width=5, height=5)
