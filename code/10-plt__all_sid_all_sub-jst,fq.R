# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
df <- mapply(
    sid=gsub("jst-(.*),.*\\.rds", "\\1", basename(args[[1]])),
    sub=gsub("jst-.*,(.*)\\.rds", "\\1", basename(args[[1]])),
    kid=lapply(ist, \(lys) unname(lys$clust)), SIMPLIFY=FALSE,
    \(sid, sub, kid) data.frame(sid, sub, kid)) |>
    do.call(what=rbind) |>
    group_by(sid, sub, kid) |>
    tally() |> mutate(p=n/sum(n))

# plotting
is <- split(seq(nrow(df)), df$sub)
gg <- lapply(names(is), \(.) {
    fd <- df[is[[.]], ]
    nc <- format(sum(fd$n), big.mark=",")
    ggplot(fd, aes(p, sid, fill=kid)) +
        guides(fill=guide_legend(ncol=1, override.aes=list(shape=21, stroke=0, size=2))) +
        geom_col(col="white", linewidth=0.1, width=1, key_glyph="point") +
        scale_x_continuous(limits=c(0, 1), n.breaks=6) +
        ggtitle(bquote(bold(.(.))~"(N ="~.(nc)*")")) +
        scale_y_discrete(limits=\(.) rev(.)) +
        scale_fill_manual(NULL, values=.pal) +
        coord_cartesian(expand=FALSE) +
        theme_bw(6) + theme(
            aspect.ratio=2/3,
            legend.margin=margin(),
            panel.grid=element_blank(),
            axis.title=element_blank(),
            legend.justification=c(0, 0.5),
            plot.title=element_text(hjust=0.5),
            legend.key.size=unit(0.2, "lines"))
}) |> wrap_plots(ncol=1)

# saving
ggsave(args[[2]], gg, units="cm", width=8, height=12)
