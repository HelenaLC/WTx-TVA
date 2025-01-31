# args <- list(
#   list.files("outs", "^ist-", full.names=TRUE),
#   "plts/ist,mds.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(scater)
    library(scuttle)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
pat <- ".*-([0-9]+)\\.rds"
sid <- gsub(pat, "\\1", args[[1]])
rds <- setNames(args[[1]], sid)
es <- lapply(rds, \(.) {
    y <- readRDS(.)$profiles
    y <- normalizeCounts(y)
})

# restrict to genes present in all sections
gs <- Reduce(intersect, lapply(es, rownames))
es <- lapply(es, \(.) .[gs, ])

# fill in missing clusters
ks <- unique(unlist(lapply(es, colnames)))
es <- lapply(es, \(.) {
    nan <- setdiff(ks, colnames(.))
    if (!length(nan)) return(.[, ks])
    add <- replicate(length(nan), numeric(nrow(.)))
    colnames(add) <- nan
    cbind(., add)[, ks]
})

# pseudobulk-level MDS
mds <- calculateMDS(do.call(cbind, es))
colnames(mds) <- c("MDS1", "MDS2")
df <- data.frame(mds,
    kid=rep(ks, length(es)),
    sid=rep(sid, each=length(ks)))

# plotting
gg <- ggplot(df, aes(MDS1, MDS2, fill=kid)) +
    guides(fill=guide_legend(override.aes=list(
        alpha=1, size=1.5, stroke=0, col=NA))) +
    geom_point(alpha=2/3, size=2, shape=21, stroke=0.1) +
    scale_fill_manual(NULL, values=.pal) +
    scale_x_continuous(breaks=seq(-9, 9, 3), expand=c(0, 1)) +
    scale_y_continuous(breaks=seq(-9, 9, 3), expand=c(0, 1)) +
    coord_equal() +
    theme_bw(6) + theme(
        plot.margin=margin(),
        legend.key=element_blank(),
        plot.background=element_blank(),
        legend.key.size=unit(0, "lines"),
        panel.grid.minor=element_blank())

# saving
ggsave(args[[2]], gg, units="cm", width=6, height=3)
