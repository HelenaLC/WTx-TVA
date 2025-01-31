# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

# loading
sid <- gsub(".*([0-9]{3}).*", "\\1", args[[1]])
ist <- lapply(setNames(args[[1]], sid), readRDS)

# wrangling
df <- lapply(ist, \(.) {
    ns <- table(k=.$clust)
    as.data.frame(ns, responseName="n")
}) |> bind_rows(.id="sid", )

# plotting
ps <- lapply(c(TRUE, FALSE), \(.) {
    if (.) {
        y <- "sid"; fill <- "k"
    } else {
        y <- "k"; fill <- "sid"
    }
    if (!.) {
        # reordering
        fd <- pivot_wider(df, 
            names_from=y, 
            values_from="n", 
            values_fill=0)
        mtx <- as.matrix(fd[, -1])
        rownames(mtx) <- fd[[1]]
        mtx <- prop.table(mtx, 1)
        hc <- hclust(dist(t(mtx)))
        yo <- colnames(mtx)[hc$order]
    } else {
        yo <- rev(sort(unique(df$sid)))
    }
    # rendering
    ggplot(df, aes(n, .data[[y]], fill=.data[[fill]])) +
        scale_fill_manual(NULL, values=unname(pals::trubetskoy())) +
        guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, col=NA, size=2))) +
        geom_col(position="fill", key_glyph="point", width=1, col="white", linewidth=0.1) +
        scale_x_continuous(labels=scales::percent_format()) +
        scale_y_discrete(limits=yo) +
        coord_equal(1/20, expand=FALSE) +
        theme_bw(6) + theme(
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            legend.key=element_blank(),
            plot.background=element_blank(),
            panel.grid.minor=element_blank(),
            legend.background=element_blank(),
            legend.key.size=unit(0, "lines"))
})

# saving
gg <- wrap_plots(ps, ncol=1, guides="collect") & 
    theme(legend.spacing=unit(0, "lines"))
n <- length(unique(df$sid))+length(unique(df$k))
ggsave(args[[2]], gg, units="cm", width=15, height=n/4)
