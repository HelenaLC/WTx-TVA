# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tidytext)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
df <- mapply(x=sce, y=ist, \(x, y) {
    cs <- match(colnames(x), names(k <- y$clust))
    data.frame(colData(x), k=k[cs]) |>
        filter(!is.na(k), !is.na(typ)) |>
        group_by(sid, typ, k) |> tally() |>
        group_by(typ) |> mutate(p=n/sum(n)) |>
        mutate(roi=typ, 
            typ=gsub(".*(REF|TVA|CRC).*", "\\1", roi),
            typ=factor(typ, c("REF", "TVA", "CRC")))
}, SIMPLIFY=FALSE) |> do.call(what=rbind)

# order by stromal frequency
xo <- df |>
    filter(k == "str") |>
    arrange(-p) |> pull(roi)
df <- mutate(df, roi=factor(roi, xo))

# label total
ns <- df |>
    group_by(typ) |>
    summarise_at("n", sum)
ns <- setNames(ns[[2]], ns[[1]])
ns <- format(ns, big.mark=",")
ns <- sprintf("%s (N = %s)", names(ns), ns)
df$typ <- factor(df$typ, labels=ns)

# cell counts
ns <- df |>
    group_by(typ, roi) |>
    summarise_at("n", sum) |>
    mutate(m=formatC(n, width=6, big.mark=",", flag=""))

# plotting
gg <- ggplot(df, aes(n, roi, fill=k)) + 
    facet_grid(rows="typ", scales="free_y", space="free_y") +
    guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=2))) +
    geom_col(width=1, alpha=2/3, col="white", linewidth=0.1, position="fill") +
    geom_text(size=1.2, hjust=1, aes(0.98, roi, label=m), ns, inherit.aes=FALSE) +
    scale_fill_manual(NULL, values=.pal_sub) +
    labs(x="frequency", y="region") +
    coord_cartesian(expand=FALSE) +
    .thm_fig_d("minimal") + theme(
        legend.position="none",
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

# saving
ggsave(args[[3]], gg, units="cm", width=5, height=7)
