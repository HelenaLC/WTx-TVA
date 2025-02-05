# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(AUCell)
    library(scuttle)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sig <- lapply(args[[1]], readRDS)
sce <- lapply(args[[2]], readRDS)

# wrangling
df <- bind_rows(mapply(SIMPLIFY=FALSE, x=sce, y=sig, \(x, y)
    data.frame(colData(x), t(assay(y)), check.names=FALSE)))

# average by regions
gs <- rownames(assay(sig[[1]]))
mu <- df |>
    filter(!is.na(typ)) |>
    group_by(typ) |>
    summarize(across(all_of(gs), mean), .groups="drop") |>
    pivot_longer(all_of(gs))

# z-scale across regions
fd <- mu |>
    group_by(name) |>
    mutate_at("value", .z) |>
    ungroup() |>
    mutate(roi=typ, typ=gsub(".*(REF|TVA|CRC).*", "\\1", roi)) |>
    mutate(typ=factor(typ, c("REF", "TVA", "CRC"))) |>
    mutate(name=gsub("HALLMARK_", "", name)) |>
    mutate(name=gsub("DESCARTES_FETAL_INTESTINE_", "", name))

# hierarchical clustering
y <- pivot_wider(select(fd, -typ))
z <- as.matrix(y[, -1])
rownames(z) <- y[[1]]

# plotting
p <- ggplot(distinct(fd, roi, typ), aes(roi, 1, fill=typ)) +
    geom_tile(show.legend=FALSE, col="white", linewidth=0.1) +
    scale_fill_manual(values=c("seagreen", "royalblue", "tomato")) +
    coord_equal(expand=FALSE) + theme_void()
q <- ggplot(fd, aes(roi, name, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(
        "z-scaled\nmean AUCell",
        low="turquoise", high="purple",
        limits=c(-2.5, 2.5), n.breaks=5) +
    labs(x="pathological region", y="gene signature") +
    scale_y_discrete(limits=.yo(z), position="right") +
    coord_equal(3/4, expand=FALSE) +
    .thm_fig_c("minimal") + theme(
        axis.title=element_blank(),
        axis.text.y=element_text(size=3),
        axis.text.x=element_text(size=4, angle=90, hjust=1, vjust=0.5))

# layout
ny <- length(unique(fd$name))
gg <- (p/q) +
    plot_layout(heights=c(1/3, ny/4), guides="collect") &
    scale_x_discrete(limits=.xo(z)) &
    theme(plot.margin=margin())

# saving
ggsave(args[[3]], gg, units="cm", width=10, height=8)
