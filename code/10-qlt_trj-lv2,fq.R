# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(slingshot)
    library(SingleCellExperiment)
})

# loading
res <- readRDS(args[[1]])
sce <- readRDS(args[[2]])
ist <- lapply(args[[3]], readRDS)

# wrangling
xy <- "Center(X|Y)_global_mm"
xy <- grep(xy, names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
#t <- grep("time", names(colData(res)))
#t <- as.matrix(colData(res)[t])
#t <- .q(rowMeans(t, na.rm=TRUE))
t <- .q(slingAvgPseudotime(res$slingshot))

# analysis
names(sub) <- names(ist) <- sub <- gsub(".*(epi|imm|str).*", "\\1", args[[3]])
idx <- lapply(ist, \(lys) intersect(colnames(sce), names(lys$clust)))
df <- lapply(sub, \(.) {
    k <- 300; i <- TRUE; if (. == "epi") { k <- k+1; i <- -1 }
    nn <- nn2(xy[idx[[.]], ], xy[idx$epi, ], searchtype="radius", r=0.05, k=k)
    is <- nn$nn.idx[, i]; is[is == 0] <- NA
    ks <- ist[[.]]$clust; ks <- factor(ks, sort(unique(ks)))
    ns <- lapply(seq(nrow(is)), \(.) table(ks[is[., ]]))
    fq <- prop.table(do.call(rbind, ns), 1)
    df <- data.frame(check.names=FALSE, t=t[idx$epi], fq)
    fd <- pivot_longer(df, -t, names_to="k", values_to="p")
}) 

# plotting
xs <- seq(-(dx <- 0.01), 1.01, 0.02)
xs <- seq(-(dx <- 0.005), 1.005, 0.01)
ps <- lapply(sub, \(.) {
    fd <- df[[.]] |>
        mutate(
            t=cut(t, breaks=xs),
            t=xs[as.integer(t)]+dx,
            t=factor(t, sort(unique(t)))) |>
        group_by(k, t) |>
        summarize_at("p", mean, na.rm=TRUE)
    ggplot(fd, aes(t, p, fill=k)) +
        geom_col(width=1, key_glyph="point", position="fill") +
        (if (. == "epi") ggtitle(.lab(wcs$sid, ncol(res)))) +
        labs(x=if (. == "str") "pseudotime", y="frequency") +
        scale_fill_manual(., values=.pal_kid) +
        coord_cartesian(expand=FALSE) 
}) 
gg <- wrap_plots(ps, ncol=1, guides="collect") &
    .thm_fig_d("minimal", "f") & theme(
        aspect.ratio=1/2,
        plot.margin=margin(),
        legend.margin=margin(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.spacing=unit(0.2, "lines"),
        legend.key.size=unit(0, "lines")) &
    guides(fill=guide_legend(override.aes=list(
        shape=21, stroke=0, size=1)))

# saving
ggsave(args[[4]], gg, units="cm", width=9, height=9)
