# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
t <- as.matrix(colData(sce)[grep("time", names(colData(sce)))])
t <- (t <- rowMeans(t, na.rm=TRUE))/max(t, na.rm=TRUE)
k <- (k <- ist$clust)[match(colnames(sce), names(k))]
df <- filter(data.frame(t, k), !is.na(t), !is.na(k))

dx <- diff(xs <- seq(-0.01, 1.01, 0.02))[1]
fd <- df |>
    mutate(
        t=cut(t, breaks=xs), 
        t=xs[as.integer(t)],
        t=factor(t, labels=xs[-1]-dx/2)) |>
    group_by(t, k) |> tally() |> mutate(p=n/sum(n))

# plotting
ggplot(fd,  aes(t, p, fill=k)) +
    guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=2))) +
    scale_y_continuous(limits=c(0, 1), n.breaks=6, labels=percent_format()) +
    labs(x="pseudotime", y="proportion of cells") +
    geom_col(position="fill", key_glyph="point") +
    scale_fill_manual(NULL, values=.pal) +
    ggtitle(.lab(wcs$sid, nrow(df))) +
    coord_equal(30, expand=FALSE) +
    theme_bw(6) + theme(
        panel.grid=element_blank(),
        legend.key=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.background=element_blank(),
        legend.key.size=unit(0, "lines"),
        legend.background=element_blank(),
        plot.title=element_text(hjust=0.5))

ggsave(args[[3]], units="cm", width=7, height=3.5)
