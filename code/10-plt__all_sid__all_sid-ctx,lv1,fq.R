# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        ctx <- readRDS(x)
        ist <- readRDS(y)
        # wrangling
        lv1 <- ist$clust[ctx$cid]
        data.frame(ctx, lv1)
    }) |> do.call(what=rbind)

# plotting
gg <- .plt_fq(df, "ctx", "lv1", hc=TRUE, h=TRUE) +
    scale_fill_manual(NULL, values=.pal_sub) +
    labs(x="niche", title=NULL) +
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(hjust=0))

# saving
ggsave(args[[3]], gg, units="cm", width=7, height=3.5)
