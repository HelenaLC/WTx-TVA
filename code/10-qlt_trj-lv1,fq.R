# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(slingshot)
    library(SingleCellExperiment)
})

ps <- mapply(x=args[[1]], y=args[[2]], z=args[[3]], \(x, y, z) {
    # loading
    res <- readRDS(x)
    sce <- readRDS(y)
    ist <- readRDS(z)
    
    # wrangling
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(sce)))
    xy <- as.matrix(colData(sce)[xy])
    t <- .q(slingAvgPseudotime(res$slingshot))
    
    # analysis
    nn <- nn2(xy, searchtype="radius", r=0.05, k=301)
    is <- nn$nn.idx[, -1]; is[is == 0] <- NA
    ks <- factor(ist$clust)
    ns <- lapply(seq(nrow(is)), \(.) table(ks[is[., ]]))
    fq <- prop.table(do.call(rbind, ns), 1)
    cs <- match(colnames(sce), colnames(res))
    df <- data.frame(check.names=FALSE, t=t[cs], fq)
    
    # plotting
    xs <- seq(-(dx <- 0.005), 1.005, 0.01)
    fd <- df |>
        filter(!is.na(t)) |>
        pivot_longer(-t, 
            names_to="k",
            values_to="p") |>
        mutate(
            t=cut(t, breaks=xs),
            t=xs[as.integer(t)]+dx,
            t=factor(t, paste((xs+dx)[-length(xs)]))) |>
        group_by(k, t) |>
        summarize_at("p", mean, na.rm=TRUE)
    
    ggplot(fd, aes(t, p, fill=k)) +
        ggtitle(.lab(paste(sce$sid[1]), ncol(sce))) +
        geom_col(width=1, key_glyph="point", position="fill") +
        labs(x="pseudotime", y="frequency") +
        scale_fill_manual(values=.pal_sub) +
        coord_cartesian(expand=FALSE) +
        .thm_fig_d("minimal", "f") + theme(
            legend.position="none",
            aspect.ratio=1/2,
            plot.margin=margin(),
            legend.margin=margin(),
            axis.text=element_blank(),
            axis.ticks=element_blank()) &
        guides(fill=guide_legend(override.aes=list(
            shape=21, stroke=0, size=1)))
})

# saving
pdf(args[[4]], onefile=TRUE, width=6/2.54, height=3.5/2.54)
for (p in ps) print(p); dev.off()
