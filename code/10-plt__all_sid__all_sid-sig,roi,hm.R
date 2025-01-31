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
df <- mapply(x=sce, y=sig, \(x, y) {
    df <- data.frame(colData(x), t(assay(y)), check.names=FALSE)
    fd <- df |>
        group_by(sid, typ) |>
        filter(!is.na(typ)) |>
        summarize(across(all_of(rownames(y)), mean), .groups="drop") |>
        pivot_longer(all_of(rownames(y))) |>
        mutate_at("sid", factor)
}, SIMPLIFY=FALSE) |> do.call(what=rbind)

# z-scale across regions
df <- df |>
    group_by(name) |>
    mutate_at("value", .z) |> ungroup() |>
    mutate(roi=typ, typ=gsub(".*(REF|TVA|CRC).*", "\\1", roi)) |>
    mutate(typ=factor(typ, c("REF", "TVA", "CRC")))

# hierarchical clustering
fd <- pivot_wider(select(df, -c(sid, typ)))
. <- as.matrix(fd[, -1]); rownames(.) <- fd[[1]]

# plotting
p <- ggplot(distinct(df, roi, typ), aes(roi, 1, fill=typ)) +
    geom_tile(show.legend=FALSE, col="white", linewidth=0.1) +
    scale_fill_manual(values=c("seagreen", "royalblue", "tomato")) +
    coord_equal(expand=FALSE) + theme_void()
q <- ggplot(df, aes(roi, name, fill=value)) +
    geom_tile() +
    coord_equal(expand=FALSE) +
    scale_fill_gradientn(
        "z-scaled\nmean AUC",
        colors=rev(hcl.colors(9, "PiYG")),
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
    labs(x="pathological region", y="gene signature") +
    scale_y_discrete(limits=.yo(.), position="right") +
    theme_bw(6) + theme(
        axis.ticks=element_blank(),
        panel.grid=element_blank(),
        axis.text=element_text(size=3.25),
        legend.key.width=unit(0.4, "lines"),
        legend.key.height=unit(0.2, "lines"),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

# layout
ny <- length(unique(df$name))
gg <- (p/q) +
    plot_layout(heights=c(1, ny), guides="collect") &
    scale_x_discrete(limits=.xo(.)) &
    theme(
        legend.position="top",
        plot.margin=margin(),
        legend.margin=margin(b=-8),
        plot.background=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=5, height=5.2)
