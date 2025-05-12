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
ccc <- readRDS(args[[3]])

# wrangling
xy <- "Center(X|Y)_global_mm"
xy <- grep(xy, names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
t <- .q(slingAvgPseudotime(res$slingshot))

# analysis
xs <- seq(-(dx <- 0.005), 1.005, 0.01)
cs <- match(colnames(sce), colnames(res))
nn <- nn2(xy, searchtype="radius", r=0.05, k=301)
is <- nn$nn.idx[, -1]
df <- data.frame(t=t[cs], i=is) |>
    filter(!is.na(t)) |>
    mutate(
        t=cut(t, breaks=xs),
        t=xs[as.integer(t)]+dx,
        t=factor(t, paste((xs+dx)[-length(xs)])))
mu <- by(df, df$t, \(fd) {
    i <- (i <- fd[, -1])[i != 0]
    colMeans((ccc$s[i, ]+ccc$r[i, ])/2)
}, simplify=FALSE) |> array2DF()
names(mu)[1] <- "t"
mu <- mu |>
    pivot_longer(-1) |>
    mutate_at("value", as.numeric)

# selection
nm <- mu |>
    group_by(name) |>
    filter(!is.na(value)) |>
    filter(grepl("-", name)) |>
    filter(!grepl("total", name)) |>
    summarise_at("value", \(.) diff(range(.))) |>
    slice_max(value, n=40) |>
    pull(name)
mv <- mu |>
    filter(!is.na(value)) |>
    filter(name %in% nm) |>
    group_by(name) |>
    mutate_at("value", .z) |>
    ungroup()

# hierarchical ordering
mx <- pivot_wider(mv)
my <- as.matrix(mx[, -1])
rownames(my) <- mx[[1]]

yo <- mv |>
    group_by(name) |>
    summarise_at("value", which.max) |>
    arrange(-value) |> pull(name)

# plotting
gg <- ggplot(mv, aes(t, name, fill=value)) +
    ggtitle(.lab(paste(sce$sid[1]), ncol(sce))) +
    scale_y_discrete(limits=yo) +
    geom_tile() +
    scale_fill_gradient2(
        "z-scaled\nmean CCC",
        limits=c(-2.5, 2.5), n.breaks=6,
        low="cadetblue", mid="ivory", high="firebrick") +
    labs(x="pseudotime", y=NULL) +
    coord_cartesian(expand=FALSE) +
    .thm_fig_c("minimal") + theme(
        plot.margin=margin(),
        legend.margin=margin(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=2))

# saving
ggsave(args[[4]], gg, units="cm", width=6, height=4)
