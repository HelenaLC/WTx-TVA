# args <- list(
#   list.files("outs", "^ist-", full.names=TRUE),
#   "plts/ist,mgs.pdf")

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
df <- lapply(rds, \(.) {
    y <- normalizeCounts(readRDS(.)$profiles)
    data.frame(g=rownames(y), y, check.names=FALSE)
}) |> bind_rows(.id="sid")

# average across sections
mu <- df |>
    group_by(g) |>
    summarize(across(
        where(is.numeric),
        \(.) mean(., na.rm=TRUE))) |>
    pivot_longer(
        where(is.numeric),
        names_to="k", values_to="y")

# restrict to selected features
gs <- list(
    epi1=c("KRT8", "CEACAM5", "CEACAM6", "KRT80", "AXIN2", "OLFM4",
        "FABP1", "SLC26A2", "CA1", "KRT19", "KRT20", "LGR5", "KCNQ1", "ASCL2"),
    epi2=c("KRT23", "REG1A", "REG1B"),
    epi3=c("DEFA5", "DEFA6", "REG3A", "REG3G", "REG4", 
        "TSLP", "CDHR2", "ESPN", "ALPI"),
    epi4=c("KRT17"),
    mast=c("KIT", "IL1RL1", "CD69"),
    BC=c("CD79A", "CD40", "CIITA", "CD19", "MS4A1", "SELL", "IGHD"),
    PC_IgA=c("MZB1", "XBP1", "IGKC", "IGHA1", "IGHA2"),
    PC_IgG=c("IGHG1", "IGHG3"), 
    PC_IgM=c("IGHM"),
    TC=c("CD2", "CD3D", "CD8A", "CD28", "TIGIT", "GATA3", "NKG7", 
        "EOMES", "GZMK", "CXCR6", "CCL5", "ITGAL", "CD40LG", "CD4"),
    mye1=c(),
    mye2=c("TLR4", "CD14", "CD68", "APOE", "APOC1", "CTSD", "MMP9"),
    SMC=c("MYL9", "MYLK", "MYH11", "CCN1"),
    EC=c("CAV1", "PECAM1", "VWF", "CD34", "CD93", "EGFL7", "KLF2"),
    fib1=c("DCN", "ADAMDEC1", "C3", "MMP2", "COL1A1", "COL1A2", "COL3A1"),
    fib2=c("COL6A1", "COL6A2", "COL6A3", "F3", "PDGFRA", "FAP", "FN1"))
fd <- mu |>
    filter(g %in% unlist(gs)) |>
    mutate_at("g", factor, unlist(gs)) |>
    mutate_at("k", factor, rev(names(gs))) |>
    # z-scale across clusters
    group_by(g) |>
    mutate_at("y", .z)

# plotting
gg <- ggplot(fd, aes(g, k, fill=y)) +
    geom_tile() + coord_equal(4/3, expand=FALSE) +
    scale_fill_gradientn(
        "z-scaled\nmean expr.",
        colors=hcl.colors(9, "Vik"),
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
    theme_bw(6) + theme(
        plot.margin=margin(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        legend.key=element_blank(),
        plot.background=element_blank(),
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5))

# saving
ggsave(args[[2]], gg, units="cm", width=12, height=3.5)
