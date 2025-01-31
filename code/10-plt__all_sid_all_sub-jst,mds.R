# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(ggplot2)
    library(patchwork)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
sid <- gsub("jst-(.*),.*\\.rds", "\\1", basename(args[[1]]))
sub <- gsub("jst-.*,(.*)\\.rds", "\\1", basename(args[[1]]))
gg <- lapply(names(idx <- split(seq_along(ist), sub)), \(sub) {
    mtx <- lapply(ist[idx[[sub]]], \(lys) lys$profiles)
    gs <- Reduce(intersect, lapply(mtx, rownames))
    mty <- do.call(cbind, lapply(mtx, \(.) .[gs, ]))
    mds <- calculateMDS(normalizeCounts(mty))
    colnames(mds) <- c("MDS1", "MDS2")
    df <- data.frame(mds, sub, 
        kid=unlist(unname(sapply(mtx, colnames))),
        sid=rep.int(sid[idx[[sub]]], sapply(mtx, ncol)))
    nc <- format(nrow(df), big.mark=",")
    ggplot(df, aes(MDS1, MDS2, col=kid, label=sid)) +
        scale_color_manual(NULL, values=.pal) +
        ggtitle(bquote(bold(.(sub))~"(N ="~.(nc)*")")) +
        guides(col=guide_legend(override.aes=list(size=1))) +
        geom_vline(xintercept=0, linewidth=0.2) +
        geom_hline(yintercept=0, linewidth=0.2) +
        geom_point() + coord_equal() + 
        theme_bw(6) + theme(
            legend.justification=c(0, 0.5),
            panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5),
            legend.key.size=unit(0.2, "lines"))
}) |> wrap_plots(ncol=1)

# saving
ggsave(args[[2]], gg, units="cm", width=10, height=18)
