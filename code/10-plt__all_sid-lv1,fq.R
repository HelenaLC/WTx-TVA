# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(ggplot2)
    library(patchwork)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
sub <- c("epi", "imm", "str")
pal <- c("gold", "cyan", "magenta")
df <- mapply(
    kid=lapply(ist, \(lys) lys$clust), SIMPLIFY=FALSE,
    sid=gsub("lv1-(.*)\\.rds", "\\1", basename(args[[1]])),
    \(sid, kid) data.frame(sid, kid)) |> do.call(what=rbind) |>
    group_by(sid, kid) |> tally() |> mutate(p=n/sum(n))
df$sid <- factor(df$sid, rev(sort(unique(df$sid))))
fd <- df |>
    group_by(sid) |>
    summarise_at("n", sum)

# plotting
p1 <- ggplot(fd, aes(n, sid)) + 
    geom_bar(stat="identity", fill="lightgrey") +
    geom_text(hjust=1.1, size=1.5,
        aes(label=format(n, big.mark=","))) +
    coord_cartesian(expand=FALSE, clip="off") +
    scale_x_reverse("number of cells", limits=c(750e3, 0)) +
    theme_minimal(6) + theme(
        aspect.ratio=a <- 2,
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank(),
        axis.title.y=element_blank())
p2 <- ggplot(df, aes(p, sid, fill=kid)) +
    guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=2))) +
    geom_col(col="white", linewidth=0.1, width=1, key_glyph="point") +
    scale_x_continuous("proportion of cells", n.breaks=6, labels=percent_format()) +
    scale_fill_manual(NULL, breaks=sub, values=pal) +
    coord_cartesian(expand=FALSE) +
    theme_bw(6) + theme(
        aspect.ratio=b <- 2/3,
        legend.position="none",
        panel.grid=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())
gg <- (plot_spacer() | p1 | p2) + 
    plot_layout(widths=c(1, 1/a, 1/b), guides="collect") &
    theme(
        plot.margin=margin(r=2),
        plot.background=element_blank())

# saving
ggsave(args[[2]], gg, units="cm", width=5.5, height=2)
