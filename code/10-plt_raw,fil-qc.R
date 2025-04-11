# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
raw <- readRDS(args[[1]])
fil <- readRDS(args[[2]])

# wrangling
ex <- !colnames(raw) %in% colnames(fil)
df <- mutate(
    data.frame(colData(raw), ex),
    `nCount_RNA/um2`=nCount_RNA/Area.um2,
    `nFeature_RNA/um2`=nFeature_RNA/Area.um2)
fd <- mutate(
    pivot_longer(df, matches("^n.*_RNA")),
    name=gsub("nCount_RNA", "counts", name),
    name=gsub("nFeature_RNA", "features", name))
mu <- fd |>
    group_by(name, ex) |>
    summarise_at("value", mean) |>
    mutate_at("value", round)

# plotting
gg <- ggplot(fd, aes(ex, value, fill=ex)) +
    scale_fill_manual(values=c("royalblue", "orange")) +
    scale_color_manual(values=c("royalblue", "orange")) +
    scale_y_continuous(limits=c(1, 5e4), trans="log10") +
    guides(fill=guide_legend("excluded", 
        override.aes=list(shape=21, stroke=NA, size=3))) +
    geom_boxplot(
        key_glyph="point", linewidth=0.2, 
        outlier.stroke=0.2, outlier.size=0.2) +
    geom_text(
        aes(y=1, col=ex, label=value), mu, 
        size=2, vjust=2, show.legend=FALSE) +
    ggtitle(bquote(bold(.(wcs$sid))~"(N ="~
            .(format(ncol(raw), big.mark=","))*";"~
            .(format(ncol(fil), big.mark=","))~"-"~
            .(round(100*mean(!ex), 2))*"%)")) +
    coord_cartesian(clip="off") +
    facet_grid(~name) +
    theme_bw(6) + theme(
        legend.position="bottom",
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.background=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5),
        legend.key.size=unit(0.5, "lines")) 

# saving
ggsave(args[[3]], gg, units="cm", width=6, height=4)
