# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(slingshot)
    library(SingleCellExperiment)
})

# get unique clusters across sections
ks <- sort(unique(unlist(lapply(args[[2]], \(.) readRDS(.)$clust))))

ps <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)
        sid <- paste(sce$sid[1])
        kid <- factor(ist$clust, ks)
        # analysis
        xy <- "Center(X|Y)_global_mm"
        xy <- grep(xy, names(colData(sce)))
        xy <- as.matrix(colData(sce)[xy])
        nn <- nn2(xy, searchtype="radius", r=0.05, k=301)
        is <- nn$nn.idx[, -1]; is[is == 0] <- NA
        ns <- lapply(seq(nrow(is)), \(.) table(kid[is[., ]]))
        fq <- prop.table(do.call(rbind, ns), 1)
        t <- .q(slingAvgPseudotime(sce$slingshot))
        df <- data.frame(check.names=FALSE, fq, t)
        # plotting
        xs <- seq(-(dx <- 0.005), 1+dx, 0.01)
        fd <- df |>
            pivot_longer(-t, 
                names_to="k",
                values_to="p") |>
            mutate(
                t=xs[as.integer(cut(t, breaks=xs))]+dx,
                t=factor(t, paste((xs+dx)[-length(xs)]))) |>
            group_by(k, t) |> summarize_at("p", mean, na.rm=TRUE)
        # plotting
        ggplot(fd, aes(t, p, fill=k)) +
            coord_cartesian(expand=FALSE) +
            ggtitle(.lab(sid, ncol(sce))) +
            labs(x="pseudotime", y="frequency") +
            scale_fill_manual(NULL, values=.pal_kid) +
            geom_col(width=1, key_glyph="point", position="fill") +
            guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=1))) +
            .thm_fig_d("minimal", "f") + theme(
                aspect.ratio=1/2,
                plot.margin=margin(),
                legend.margin=margin(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                panel.grid=element_blank())
    })

# saving
pdf(args[[3]], onefile=TRUE, width=8/2.54, height=3.5/2.54)
for (p in ps) print(p); dev.off()
