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
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
.sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
auc <- mapply(
    SIMPLIFY=FALSE, x=args[[1]], 
    y=split(args[[2]], .sid(args[[2]])), \(x, y) {
        sce <- readRDS(x)
        ist <- lapply(y, readRDS)
        kid <- lapply(ist, \(.) .$clust)
        sub <- rep.int(.sub(y), sapply(kid, length))
        kid <- unlist(kid); names(sub) <- names(kid)
        idx <- intersect(colnames(sce), names(kid))
        sid <- .sid(x); sub <- sub[idx]; kid <- kid[idx]
        SummarizedExperiment(
            list(counts=assay(sce[, idx])),
            colData=data.frame(sid, sub, kid))
    }) |> do.call(what=cbind)

pat <- "^(HALLMARK_|DESCARTES_FETAL_INTESTINE_)"
rownames(auc) <- gsub(pat, "", rownames(auc))

# plotting
.p <- \(x, i, xo=TRUE, yo=TRUE) {
    n <- length(unique(x$kid))
    x <- select(x, kid, name, value)
    y <- pivot_wider(x, names_from="kid")
    z <- as.matrix(y[, -1]); rownames(z) <- y[[1]]
    ggplot(x, aes(kid, name, fill=value)) +
        (if (xo) scale_x_discrete(limits=.yo(z))) +
        (if (yo) scale_y_discrete(limits=.xo(z))) +
        scale_fill_gradient2(
            "z-scaled\nmean AUCell",
            low="turquoise", high="purple",
            limits=c(-2.5, 2.5), n.breaks=5) +
        ggtitle(bquote(bold(.(i))~"(N ="~.(n)*")")) +
        geom_tile() + coord_equal(expand=FALSE) +
        .thm_fig_c("bw") + theme(
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
}

# joint
cd <- colData(auc)[ids <- c("sub", "kid")]
mu <- aggregateAcrossCells(auc, cd, statistics="mean")
df <- data.frame(colData(mu)[ids], t(assay(mu)))
fd <- df |>
    pivot_longer(all_of(rownames(mu))) |>
    group_by(sub, name) |>
    mutate(value=.z(value)) |>
    ungroup()
p0 <- by(fd, fd$sub, \(.) {
    .p(., i <- .$sub[1]) + 
    if (i != "epi") theme(axis.text.y=element_blank())
}) |> wrap_plots(nrow=1) + plot_layout(guides="collect")

# split
cd <- colData(auc)[ids <- c("sub", "sid", "kid")]
mu <- aggregateAcrossCells(auc, cd, statistics="mean")
df <- data.frame(colData(mu)[ids], t(assay(mu)))

# by sample
fd <- df |>
    pivot_longer(all_of(rownames(mu))) |>
    group_by(sid, name) |>
    mutate(value=.z(value)) |>
    select(-sub) |>
    ungroup()
ps <- by(fd, fd$sid, \(.) .p(., .$sid[1]))

# by subset
fd <- df |>
    pivot_longer(all_of(rownames(mu))) |>
    group_by(sid, sub, name) |>
    mutate(value=.z(value)) |>
    ungroup()
qs <- by(fd, fd$sid, \(.) by(., .$sub, \(.) {
    .p(., paste0(.$sid[1], ": ", i <- .$sub[1])) + 
    if (i != "epi") theme(axis.text.y=element_blank())
}) |> wrap_plots(nrow=1) + plot_layout(guides="collect"))

# saving
pdf(args[[3]], width=15/2.54, height=12/2.54, onefile=TRUE)
for (p in c(list(p0), ps, qs)) print(p); dev.off()
