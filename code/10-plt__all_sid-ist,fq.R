# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(ggplot2)
})

# loading & wrangling
df <- lapply(args[[1]], \(rds) {
    sid <- gsub(".*-(.*)\\.rds", "\\1", basename(rds))
    kid <- (ist <- readRDS(rds))$clust
    data.frame(sid, table(kid))
}) |> do.call(what=rbind) |>
    dplyr::rename(n=Freq) |>
    group_by(sid) |>
    mutate(p=n/sum(n)) |>
    mutate_at("sid", factor)

# plotting by section
p <- ggplot(df, aes(p, sid, fill=kid)) +
    geom_col(col="white", linewidth=0.1, width=1, key_glyph="point") +
    scale_x_continuous("proportion of cells", 
        n.breaks=6, labels=percent_format()) +
    scale_y_discrete(limits=\(.) rev(.)) +
    coord_cartesian(expand=FALSE, clip="off") +
    theme_bw(6) + theme(
        panel.grid=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

# plotting by cluster
fd <- df |>
    group_by(kid) |>
    summarise_at("n", sum) |>
    arrange(-n) |>
    mutate(kid=factor(kid, kid))
q <- ggplot(fd, aes(kid, n, fill=kid)) +
    geom_col(col="white", linewidth=0.1, width=1, key_glyph="point") +
    scale_y_continuous("number of cells", 
        limits=c(0, 8e5), breaks=seq(0, 8e5, 2e5), 
        labels=scales::label_number(scale=1e-3, suffix="k")) +
    coord_cartesian(expand=FALSE) +
    theme_bw(6) + theme(
        aspect.ratio=2/3,
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.text.x=element_text(size=3.5, angle=90, hjust=1, vjust=0.5))

# aesthetics
ps <- lapply(list(p, q), \(.) . + theme(
    legend.key=element_blank(),
    plot.background=element_blank(),
    legend.key.size=unit(0, "lines"),
    legend.text=element_text(size=4),
    legend.background=element_blank()) +
    scale_fill_manual(NULL, breaks=levels(df$kid), values=.pal) +
    guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=1.25))))

# saving
pdf(args[[2]], onefile=TRUE, width=6/2.54, height=2.5/2.54)
for (p in ps) print(p); dev.off()
