# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tidytext)
    library(HDF5Array)
    library(SingleCellExperiment)
})

args <- list(
    list.files("outs", "fil", full.names=TRUE),
    list.files("outs", "lv2", full.names=TRUE),
    "plts/fil,lv2,ds.pdf")

# loading
df <- mapply(SIMPLIFY=FALSE, sce=args[[1]], 
    ist=grep("epi", args[[2]], value=TRUE), 
    \(sce, ist) {
        sce <- readRDS(sce)
        ist <- readRDS(ist)
        # get nearest neighbors
        xy <- grep("global_mm", names(colData(sce)))
        xy <- as.matrix(colData(sce)[xy])
        ns <- nn2(xy, k=301, r=0.05, searchtype="radius")
        is <- ns$nn.idx[, -1]; is[is == 0] <- NA
        # quantify by subpopulation
        ks <- (ks <- ist$clust)[match(rownames(xy), names(ks))]
        ks[is.na(ks)] <- "" # exclude non-epithelia
        nn <- matrix(ks[c(is)], nrow(is), ncol(is))
        ns <- rowSums(nn != "", na.rm=TRUE)
        df <- data.frame(
            sid=factor(sce$sid[1]), 
            kid=ks, n=ns, check.names=FALSE)
    }) |> do.call(what=rbind)

# wrangling
fd <- df |>
    filter(kid != "", !is.na(n))

# downsample to similar number of
# cells per section-subpopulation
fd <- by(fd, fd$kid, \(a) {
    n <- 1e4
    by(a, a$sid, \(b) {
        b[sample(nrow(b), min(nrow(b), 1e4)), ]
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind) |>
    mutate(m=.z(n)) |>
    mutate(kid=gsub("^epi\\.", "", kid))

# plotting
gg <- ggplot(fd,
    aes(reorder(kid, m, median), m, fill=kid)) + 
    geom_boxplot(
        outlier.shape=16, outlier.size=0.2, outlier.stroke=0,
        key_glyph="point", show.legend=TRUE, alpha=2/3, linewidth=0.2) +
    labs(x=NULL, y="cellular density (50um radius)") +
    scale_fill_manual(NULL, drop=FALSE, values=.pal_kid) +
    scale_y_continuous(n.breaks=6, limits=c(-2.5, 2.5)) +
    geom_hline(yintercept=0, linewidth=0.2) +
    scale_x_reordered() +
    .thm_fig_d("bw", "f") + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())

# saving
ggsave(args[[3]], gg, unit="cm", width=6, height=5)
