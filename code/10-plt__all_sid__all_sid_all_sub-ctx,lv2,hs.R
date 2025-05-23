# args <- list(
#     list.files("outs", "ctx-", full.names=TRUE),
#     list.files("outs", "lv2-", full.names=TRUE),
#     "plts/ctx,lv2,hs.pdf")
# args <- lapply(args, \(.) grep("241", ., value=TRUE, invert=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# Shannon entropy
.h <- \(.) { . <- .[. > 0]; -sum(.*log(.)) }
.by_id <- \(.) split(., gsub(".*([0-9]{3}).*", "\\1", .))

df <- mapply(
    x=args[[1]], y=.by_id(args[[2]]), 
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        df <- readRDS(x)
        ks <- lapply(y, readRDS)
        # wrangling
        ks <- unlist(lapply(ks, \(.) .$clust))
        cs <- intersect(df$cid, names(ks))
        df <- df |>
            filter(cid %in% cs) |>
            mutate(kid=ks[cs])
        mx <- select(df, where(is.numeric))
        h <- apply(as.matrix(mx), 1, .h)
        data.frame(df, h)
    }) |> do.call(what=rbind)

# downsample to similar number of
# cells per section-subpopulation
df <- by(df, df$kid, \(a) {
    n <- 1e4
    by(a, a$sid, \(b) {
        b[sample(nrow(b), min(nrow(b), 1e4)), ]
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# plotting
df$kid <- gsub("epi\\.", "", df$kid)
ps <- by(df, df$sub, \(fd) {
    ggplot(fd, aes(reorder(kid, h, median, na.rm=TRUE), h, fill=kid)) +
        geom_boxplot(
            key_glyph="point", alpha=2/3, linewidth=0.2,
            outlier.shape=16, outlier.size=0.2, outlier.stroke=0) +
        scale_y_continuous("entropy", limits=c(0, 3)) +
        geom_hline(yintercept=median(fd$h), linewidth=0.2) +
        scale_fill_manual(fd$sub[1], values=.pal_kid) +
        .thm_fig_d("bw", "f") + theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            panel.grid.minor=element_blank())
})

# saving
gg <- wrap_plots(ps, nrow=1)
ggsave(args[[3]], gg, units="cm", width=16, height=5)
