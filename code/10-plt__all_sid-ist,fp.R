# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(InSituType)
})

# loading
ist <- lapply(args[[1]], readRDS)

# plotting
pat <- ".*-(.*)\\.rds"
names(ist) <- gsub(pat, "\\1", args[[1]])
gg <- lapply(names(ist), \(.) .plt_fp(ist[[.]], .))
gg <- wrap_plots(gg, nrow=2) & theme(plot.margin=margin(r=2))

# saving
ggsave(args[[2]], gg, units="cm", width=12, height=7)
