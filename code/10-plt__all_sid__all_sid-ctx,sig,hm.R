# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(AUCell)
    library(scuttle)
})

# loading
ctx <- lapply(args[[1]], readRDS)
sig <- lapply(args[[2]], readRDS)

# wrangling
df <- mapply(x=ctx, y=sig, \(x, y) {
    # average by spatial context
    mu <- aggregateAcrossCells(y[, x$cid], ids=x$ctx,
        use.assay.type=assayNames(y), statistics="mean")
    # correct for number of genes in set
    ns <- sapply(rowData(y)$set, length)
    sid <- mu$sid[1]; mu <- assay(mu)/ns
    data.frame(sid, set=rownames(mu), mu)
}, SIMPLIFY=FALSE) |> bind_rows(.id="sid")

# average across samples, by gene set & spatial context
fd <- df |>
    pivot_longer(
        -c(sid, set), 
        names_to="ctx", 
        values_to="auc") |>
    group_by(set, ctx) |>
    summarise_at("auc", mean, na.rm=TRUE) |>
    group_by(set) |>
    mutate_at("auc", .z) |>
    mutate(set=gsub("HALLMARK_", "", set)) |>
    mutate(set=gsub("DESCARTES_FETAL_INTESTINE_", "", set))

# order according to hierarchical clustering
mx <- pivot_wider(fd, names_from="set", values_from="auc")
my <- as.matrix(mx[, -1]); rownames(my) <- mx[[1]]
xo <- .xo(my); yo <- .yo(my)

# plotting
gg <- ggplot(fd, aes(ctx, set, fill=auc)) +
    geom_tile() +
    scale_fill_gradient2(
        "z-scaled\nmean AUCell",
        low="turquoise", high="purple",
        limits=c(-2.5, 2.5), n.breaks=5) +
    scale_x_discrete(limits=xo) +
    scale_y_discrete(limits=yo) +
    coord_equal(3/4, expand=FALSE) +
    .thm_fig_c("minimal") + theme(
        axis.title=element_blank(),
        axis.text.y=element_text(size=3),
        axis.text.x=element_text(size=4, angle=90, hjust=1, vjust=0.5))

# saving
ggsave(args[[3]], gg, units="cm", width=10, height=8)