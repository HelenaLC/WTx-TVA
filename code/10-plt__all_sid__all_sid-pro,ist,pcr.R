# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# analysis
xs <- c(kid="k", area="Area",
    CD45="Mean.CD45", CD68="Mean.CD68_CK8_18",
    DAPI="Mean.DAPI", PanCK="Mean.PanCK", PanMem="Mean.CD298_B2M")

df <- mapply(x=sce, y=ist, \(x, y) {
    x$k <- (k <- y$clust)[match(colnames(x), names(k))]
    df <- data.frame(.pcr(x, xs), sid=x$sid[1])
    df$x <- factor(df$x, xs, names(xs)); df
}, SIMPLIFY=FALSE) |> do.call(what=rbind)

ve <- lapply(sce, \(.) {
    pcs <- reducedDim(., "PCA")
    ve <- attr(pcs, "varExplained")
    ve <- round(sum(ve), digits=2)
    ve <- paste0(format(ve, nsmall=2), "%")
    data.frame(sid=.$sid[1], ve)
}) |> do.call(what=rbind)

# aesthetics
pal <- c(
    kid="red",
    area="grey",
    DAPI="blue",
    CD45="magenta",
    CD68="gold",
    PanCK="cyan2",
    PanMem="limegreen")

# plotting
gg <- ggplot(df, aes(pc, r2, col=x)) + 
    geom_line(linewidth=0.3, show.legend=FALSE) + 
    geom_point(size=0.6) + facet_wrap(~sid, nrow=2) + 
    scale_x_continuous("principal component", breaks=c(1, seq(5, 100, 5))) +
    scale_y_continuous("R-squared from LM", limits=c(0, 1), n.breaks=3) +
    geom_label(size=2, hjust=1, vjust=1, col="grey50", 
        aes(Inf, Inf, label=ve), ve, inherit.aes=FALSE) +
    guides(col=guide_legend(override.aes=list(size=1))) +
    scale_color_manual("predictor", values=pal) +
    coord_cartesian(xlim=c(1, 15)) +
    theme_bw(6) + theme(
        plot.background=element_blank(),
        panel.grid.minor=element_blank(),
        legend.key.size=unit(0, "lines"),
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        panel.background=element_rect(fill=NA))

# saving
ggsave(args[[3]], gg, units="cm", width=15, height=5)
