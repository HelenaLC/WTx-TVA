# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- lapply(args[[2]], readRDS)

# wrangling
sub <- gsub(".*(epi|imm|str).*", "\\1", basename(args[[2]]))
names(ist) <- names(sub) <- sub
df <- lapply(sub, \(.) {
    idx <- names(kid <- ist[[.]]$clust)
    data.frame(colData(sce)[idx, ], kid)
}) |> bind_rows(.id="sub") |>
    mutate(roi=gsub("^.*_", "", typ)) |>
    mutate(roi=factor(roi, names(.pal_roj))) |>
    mutate(roi=droplevels(roi))

# plotting
ps <- by(df, df$sub, \(fd) {
    gg <- .plt_fq(fd, "roi", "kid", hc=FALSE) +
        coord_equal(15, expand=FALSE) +
        .thm_fig_d("minimal", "f") + theme(
            legend.margin=margin(),
            legend.title=element_blank(),
            plot.title=element_blank(),
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    gg$guides$guides$fill$params$override.aes$size <- 1.5; gg
})

# saving
gg <- wrap_plots(ps, nrow=1)
ggsave(args[[3]], gg, units="cm", width=12, height=4)
