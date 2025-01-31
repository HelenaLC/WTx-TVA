# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
sub <- c("epi", "imm", "str")
pal <- c("gold", "cyan", "magenta")
df <- mapply(
    kid=lapply(ist, \(lys) lys$clust), SIMPLIFY=FALSE,
    sid=gsub("lv1-(.*)\\.rds", "\\1", basename(args[[1]])),
    \(sid, kid) data.frame(sid, as.data.frame(table(kid), responseName="n"))) |>
    do.call(what=rbind)
fd <- df |>
    group_by(sid) |>
    summarise_at("n", sum) |>
    mutate(kid="imm")

# plotting
gg <- ggplot(df, aes(kid, log10(n), fill=kid)) + facet_grid(~sid) +
    geom_bar(key_glyph="point", width=1, stat="identity", position="dodge") +
    geom_text(
        aes(label=format(n, big.mark=",")),
        size=2, y=4, angle=90, hjust=1, vjust=0.5) +
    guides(fill=guide_legend(override.aes=list(size=2, shape=21))) +
    geom_text(
        aes(kid, label=format(n, big.mark=",")),
        size=2, y=5.5, hjust=0.5, data=fd, inherit.aes=FALSE) +
    scale_y_continuous("# cells (log10)", limits=c(0, 6)) +
    scale_fill_manual(NULL, breaks=sub, values=pal) +
    coord_cartesian(expand=FALSE, ylim=c(3, 6)) +
    theme_bw(6) + theme(
        aspect.ratio=3,
        plot.margin=margin(),
        legend.margin=margin(),
        legend.position="bottom",
        panel.grid=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        legend.key.size=unit(0.5, "lines"))

# saving
ggsave(args[[2]], gg, units="cm", width=9, height=4.5)
