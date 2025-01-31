# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tidytext)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
# 
# analysis
df <- lapply(sce, \(.) {
    pcs <- reducedDim(., "PCA")
    rot <- attr(pcs, "rotation")
    df <- lapply(colnames(pcs)[seq_len(3)], \(pc) {
        top <- order(-abs(rot[, pc])) <= 20
        data.frame(pc,
            sid=.$sid[1],
            rot=rot[top, pc],
            g=rownames(rot)[top])
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)
nc <- max(nchar(df$g))+1
df$g <- formatC(df$g, width=nc, flag = " ")

# plotting
gg <- lapply(unique(df$pc), \(.) {
    ggplot(filter(df, pc == .), aes(rot, fill=rot > 0,
        y=reorder_within(g, abs(rot), sid))) + 
        geom_col(alpha=1/2, show.legend=FALSE) +
        geom_vline(xintercept=0, linewidth=0.2) +
        facet_wrap(~sid, nrow=1, scales="free") +
        scale_fill_manual(values=c("blue", "red")) +
        labs(x=NULL) + scale_y_reordered(.) +
        coord_cartesian(expand=FALSE) +
        theme_bw(6) + theme(
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.background=element_blank()) +
        if (. != "PC1") theme(strip.text=element_blank()) else 
        theme(strip.text=element_text(face="bold"))
}) |> wrap_plots(ncol=1) +  
    plot_layout(guides="collect") & 
    theme(
        legend.position="none",
        plot.margin=margin(b=4),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        strip.background=element_blank(),
        axis.text.y=element_text(size=3))

# saving
ggsave(args[[3]], gg, units="cm", width=15, height=7)
