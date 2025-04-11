# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(ggnewscale)
    library(RColorBrewer)
    library(SingleCellExperiment)
})

# loading
get_sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
get_sub <- \(.) gsub(".*(epi|imm|str).*", "\\1", .)
df <- lapply(args[[1]], \(x) {
    sce <- readRDS(x)
    i <- get_sid(args[[2]]) == get_sid(x)
    lapply(args[[2]][i], \(y) {
        kid <- readRDS(y)$clust
        idx <- match(colnames(sce), names(kid))
        data.frame(
            colData(sce)[c("sid", "typ")], 
            sub=get_sub(y), kid=kid[idx])
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind)

# wrangling
df <- df |>
    filter(!is.na(typ), !is.na(kid)) |>
    group_by_all() |> 
    tally() |> 
    mutate(p=n/sum(n)) |> 
    ungroup() |>
    dplyr::rename(roi=typ) |>
    mutate_at("sid", factor) |>
    mutate(typ=gsub(".*(REF|TVA|CRC).*", "\\1", roi)) |>
    mutate(typ=factor(typ, c("REF", "TVA", "CRC")))

# plotting
gd <- guides(fill=guide_legend(ncol=1, 
    override.aes=list(shape=21, size=1, stroke=0)))
ps <- lapply(group_split(df, sub), \(fd) {
    sub <- fd$sub[1]
    x <- pivot_wider(
        select(ungroup(fd), roi, kid, p),
        names_from="kid", values_from="p")
    y <- as.matrix(x[, -1]); rownames(y) <- x[[1]]
    nc <- format(sum(fd$n), big.mark=",")
    p1 <- ggplot(fd, aes(roi, p, fill=kid)) + 
        geom_col(width=1, key_glyph="point") +
        scale_fill_manual(sub, values=.pal_kid) +
        coord_cartesian(expand=FALSE) +
        labs(y="frequency") +
        gd + theme_bw(6) + theme(
            plot.margin=margin(),
            axis.ticks=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank()) +
        if (sub != "str") theme(axis.text=element_blank()) else
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    p2 <- ggplot(distinct(fd, roi, .keep_all=TRUE), aes(roi)) + 
        geom_tile(key_glyph="point", col="white", aes(y=0, fill=typ)) + 
        scale_fill_manual(values=.pal_roi) +
        ggtitle(bquote(bold(.(sub))~"(N ="~.(nc)*")")) +
        gd + new_scale_fill() + gd +
        geom_tile(key_glyph="point", col="white", aes(y=1, fill=sid)) + 
        scale_fill_manual(values=.pal_sid) +
        coord_equal(2/3, expand=FALSE) + theme_void(6) + 
        theme(plot.title=element_text(hjust=0.5)) +
        if (sub != "epi") theme(legend.position="none")
    (p2/p1) + 
        plot_layout(heights=c(1, 7)) & 
        scale_x_discrete(NULL, limits=\(.) .xo(y))
})

# saving
gg <- wrap_plots(ps, ncol=1) +
    plot_layout(guides="collect") & 
    # keep order fixed across subsets
    tail(ps[[1]][[2]]$scales$scales, 1) &
    theme(
        plot.background=element_blank(),
        panel.background=element_blank(),
        legend.background=element_blank(),
        legend.spacing=unit(0.2, "lines"),
        legend.key.size=unit(0, "lines"),
        legend.text=element_text(size=4),
        legend.key=element_blank(),
        panel.grid=element_blank(),
        legend.margin=margin())
ggsave(args[[3]], gg, units="cm", width=8, height=12)
