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
df <- lapply(ist, \(.)
    data.frame(kid=.$clust)) |>
    bind_rows(.id="sid")

# plotting
ps <- lapply(c(TRUE, FALSE), \(.) {
    if (.) {
        y <- "sid"
        f <- "kid"
        pal <- .pal_kid
    } else {
        y <- "kid"
        f <- "sid"
        pal <- .pal_sid
    }
    ns <- as.data.frame(table(df[[y]]))
    ns <- setNames(ns, c(y, "n"))
    gg <- .plt_fq(df, y, f, id="epi", h=TRUE) +
        geom_text(
            aes(.data[[y]], 1.2, label=format(n, big.mark=",")), 
            ns, inherit.aes=FALSE, hjust=1, size=1.5) +
        scale_fill_manual(NULL, values=pal) +
        guides(fill=guide_legend(override.aes=list(
            alpha=1, shape=21, stroke=0, size=1.5))) +
        theme(aspect.ratio=length(unique(df[[y]]))/20) +
        (if (!.) theme(plot.title=element_blank()))
    gg$coordinates$clip <- "off"
    gg
})

# saving
gg <- wrap_plots(ps, ncol=1, 
    guides="collect") & theme(
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.text.x=element_blank(),
        legend.spacing=unit(0, "lines"))
n <- length(unique(df$sid))+length(unique(df$kid))
ggsave(args[[2]], gg, units="cm", width=10, height=n/4)
