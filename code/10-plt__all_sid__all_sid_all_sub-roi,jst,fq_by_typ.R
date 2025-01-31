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
df <- mapply(x=rep(args[[1]], each=3), y=args[[2]], \(x, y) {
    sce <- readRDS(x); ist <- readRDS(y)
    sub <- gsub(".*(epi|imm|str).*", "\\1", y)
    kid <- (k <- ist$clust)[match(colnames(sce), names(k))]
    data.frame(colData(sce)[c("sid", "typ")], sub, kid)
}, SIMPLIFY=FALSE) |> do.call(what=rbind) 

# wrangling
fd <- df |>
    mutate(typ=gsub(".*(REF|TVA|CRC).*", "\\1", typ)) |>
    mutate(typ=factor(typ, c("REF", "TVA", "CRC"))) |>
    filter(!is.na(typ), !is.na(kid)) |> 
    # count cells by cluster & region
    group_by(sub, kid, typ) |> tally() |>
    # get relative proportion by cluster
    mutate(p=n/sum(n)) |> ungroup()
    
# order by REF fraction
xo <- fd |>
    filter(typ == "REF") |> 
    arrange(p) |> pull(kid) |> unique()
fd <- mutate(fd, kid=factor(kid, xo))

# plotting
gg <- ggplot(fd, aes(p, kid, fill=typ)) + 
    facet_grid(rows="sub", scales="free_y", space="free") +
    geom_col(alpha=2/3, position="fill", width=1, key_glyph="point") +
    scale_fill_manual(values=c("limegreen", "royalblue", "tomato")) +
    scale_x_continuous(
        "proportion of cells", n.breaks=6, 
        labels=scales::label_percent()) +
    coord_cartesian(expand=FALSE) +
    theme_bw(6) + theme(
        legend.position="none",
        panel.grid=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        strip.background=element_rect(fill=NA))

# saving
ggsave(args[[3]], gg, units="cm", width=6, height=8)
