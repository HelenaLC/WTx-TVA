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

# loading
epi <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# get nearest neighbors
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
ns <- nn2(xy, k=301, r=0.05, searchtype="radius")
is <- ns$nn.idx[, -1]; is[is == 0] <- NA

# count epithelial neighbors
ks <- rep("foo", ncol(sce))
names(ks) <- colnames(sce)
ks[colnames(epi)] <- "epi"
nn <- matrix(ks[c(is)], nrow(is), ncol(is))
ns <- rowSums(nn == "epi", na.rm=TRUE)

# quantify by binned pseudotime
t <- .q(slingshot::slingAvgPseudotime(epi$slingshot))
t <- t[match(colnames(sce), colnames(epi))]
df <- data.frame(n=ns, k=ks, t)
xs <- seq(-(dx <- 0.005), 1.005, 0.01)
mu <- df |>
    filter(!is.na(t)) |>
    mutate(
        t=cut(t, breaks=xs),
        t=xs[as.integer(t)]+dx,
        t=factor(t, sort(unique(t)))) |>
    group_by(t) |>
    summarize_at("n", mean, na.rm=TRUE) |>
    mutate_at("n", .z)

# plotting
gg <- ggplot(mu, aes(t, n, group="")) + 
    geom_hline(yintercept=0, linewidth=0.2) +
    geom_path(col="red", linewidth=0.4) +
    ggtitle(.lab(wcs$x, ncol(epi))) +
    scale_x_discrete(breaks=seq(0, 1, 0.2)) +
    labs(x="pseudotime", y="epithelial density\n(50um radius)") +
    .thm_fig_c("minimal") + theme(
        axis.text.x=element_blank(),
        panel.grid.minor=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=5, height=3)
