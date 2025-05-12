# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(AUCell)
    library(ggplot2)
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
sig <- readRDS(args[[2]])

# wrangling
sub <- rowData(sig)$sub
epi <- sapply(sub, \(.) any(. == "epi"))
y <- assay(sig)[which(epi), colnames(sce)]
t <- .q(slingAvgPseudotime(sce$slingshot))
df <- pivot_longer(data.frame(t, t(y)), -t)

# average by binned pseudotime
xs <- seq(-(dx <- 0.005), 1.005, 0.01)
mu <- df |>
    mutate(
        t=cut(t, breaks=xs),
        t=xs[as.integer(t)]+dx,
        t=factor(t, sort(unique(t)))) |>
    group_by(name, t) |>
    summarize_at("value", mean, na.rm=TRUE) |>
    mutate_at("value", .z)

# hierarchical clustering
y <- pivot_wider(mu)
z <- as.matrix(y[, -1])
rownames(z) <- y[[1]]

# plotting
gg <- ggplot(mu, aes(t, name, fill=value)) + geom_tile() + 
    scale_fill_gradient2(
        "z-scaled\nmean AUCell",
        low="turquoise", high="purple",
        limits=c(-2.5, 2.5), n.breaks=5) +
    geom_vline(xintercept=paste(seq(0.2, 0.8, 0.2)), linewidth=0.1) +
    scale_y_discrete(position="right", limits=.yo(z)) +
    ggtitle(.lab(wcs$x, ncol(sce))) +
    coord_equal(2, expand=FALSE) +
    labs(x="pseudotime", y=NULL) +
    .thm_fig_c("minimal") + theme(
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=3),
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"))

# saving
ggsave(args[[3]], gg, units="cm", width=10, height=8)
