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
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
df <- mapply(x=rep(args[[1]], each=3), y=args[[2]], \(x, y) {
    sce <- readRDS(x); ist <- readRDS(y)
    sub <- gsub(".*(epi|imm|str).*", "\\1", y)
    kid <- (k <- ist$clust)[match(colnames(sce), names(k))]
    data.frame(colData(sce)[c("sid", "typ")], sub, kid) |>
        filter(!is.na(typ), !is.na(kid)) |>
        group_by_all() |>
        tally() |>
        mutate(p=n/sum(n))
}, SIMPLIFY=FALSE) |>
    do.call(what=rbind) |>
    dplyr::rename(roi=typ) |>
    mutate(typ=gsub(".*(REF|TVA|CRC).*", "\\1", roi)) |>
    mutate(typ=factor(typ, c("REF", "TVA", "CRC"))) |>
    mutate_at("sid", factor)

# aesthetics
gd <- guides(fill=guide_legend(ncol=1, override.aes=list(shape=21, size=1, stroke=0)))

# plotting
ps <- lapply(group_split(ungroup(df), sub), \(fd) {
    sub <- fd$sub[1]
    x <- pivot_wider(
        select(ungroup(fd), roi, kid, p),
        names_from="kid", values_from="p")
    y <- as.matrix(x[, -1]); rownames(y) <- x[[1]]
    nc <- format(sum(fd$n), big.mark=",")
    p1 <- ggplot(fd, aes(roi, p, fill=kid)) + 
        geom_col(width=1, key_glyph="point") +
        scale_y_continuous(
            "proportion of cells", n.breaks=6, 
            labels=scales::label_percent()) +
        scale_fill_manual(sub, values=.pal) +
        coord_cartesian(expand=FALSE) +
        gd + theme_bw(6) + theme(
            plot.margin=margin(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank())
    p2 <- ggplot(distinct(fd, roi, .keep_all=TRUE), aes(roi)) + 
        geom_tile(key_glyph="point", col="white", aes(y=0, fill=typ)) + 
        scale_fill_manual(values=c("limegreen", "royalblue", "tomato")) +
        ggtitle(bquote(bold(.(sub))~"(N ="~.(nc)*")")) +
        gd + new_scale_fill() + gd +
        geom_tile(key_glyph="point", col="white", aes(y=1, fill=sid)) + 
        scale_fill_manual(values=unname(pals::polychrome(8))) +
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
    ps[[1]][[2]]$scales$scales[[3]] &
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
ggsave(args[[3]], gg, units="cm", width=8, height=10)
