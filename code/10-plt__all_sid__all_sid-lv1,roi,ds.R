# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
df <- mapply(
    ist=args[[1]], sce=args[[2]], 
    SIMPLIFY=FALSE, \(ist, sce) {
        ist <- readRDS(ist)
        sce <- readRDS(sce)
        # get nearest neighbors
        xy <- grep("global_mm", names(colData(sce)))
        xy <- as.matrix(colData(sce)[xy])
        ns <- nn2(xy, k=301, r=0.05, searchtype="radius")
        is <- ns$nn.idx[, -1]; is[is == 0] <- NA
        ks <- (ks <- ist$clust)[match(rownames(xy), names(ks))]
        nn <- matrix(ks[c(is)], nrow(is), ncol(is))
        # quantify by subset & region
        ks <- sort(unique(ist$clust))
        ns <- sapply(ks, \(.) rowSums(nn == ., na.rm=TRUE))
        data.frame(ns, roi=sce$typ, sid=factor(sce$sid[1]))
}) |> do.call(what=rbind)

# wrangling
fd <- df |>
    pivot_longer(
        where(is.numeric), 
        names_to="sub", 
        values_to="n") |>
    group_by(sid, sub) |> 
    mutate(m=.z(n)) |>
    filter(!is.na(roi)) |>
    mutate(roj=gsub("^.*_", "\\1", roi)) |>
    mutate(roj=factor(roj, names(.pal_roj)))

# plotting
gg <- ggplot(fd, aes(roj, m, fill=roj)) + 
    geom_boxplot(
        outlier.shape=16, outlier.size=0.2, outlier.stroke=0,
        key_glyph="point", show.legend=TRUE, alpha=2/3, linewidth=0.2) +
    scale_fill_manual("region", drop=FALSE, values=.pal_roj) +
    facet_grid(sub~sid, space="free", scales="free_x") +
    labs(x=NULL, y="cellular density (50um radius)") +
    #scale_y_continuous(limits=c(-2, 2.5)) +
    geom_hline(yintercept=0, linewidth=0.2) +
    .thm_fig_d("bw", "f") + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())

# saving
ggsave(args[[3]], gg, unit="cm", width=10, height=6)
