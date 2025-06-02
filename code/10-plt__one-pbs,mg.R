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
sce <- readRDS(args[[1]])

# restrict to selected features
gs <- list(
    epi1=c("CEACAM5", "CEACAM6", "PIGR", "KRT8", 
        "FABP1", "CA1", "CA2", "CA7", "PYY"),
    epi2=c("CLDN2", "LGR5", "CDX2", "KRT23"),        
    epi3=c("REG3A", "DEFA5", "DEFA6"),
    epi4=c("KRT19", "KRT20", "OLFM4", "KRT80", "DMBT1", "CTSE", "MSLN"),
    mast=c("KIT", "IL1RL1", "CD69", "HDC", "GATA2", "SLC18A2"),
    TC=c("CD2", "CD3D", "CD8A", "CD28", "TIGIT", "GATA3", "NKG7", "CXCR6", "CCL5", "ITGAL", "CD40LG"),
    BC=c("SELL", "CD79A", "CD19", "MS4A1", "CR2"),
    PC_IgA=c("IGHD", "MZB1", "XBP1", "IGKC", "IGHA1"),
    PC_IgG=c("IGHG1"),
    PC_IgM=c("IGHM"),
    "mye1"=c("CD4", "TLR4", "CD14", "CD68", "APOE", "APOC1", "CTSD"),
    "mye2"=c("S100A8", "S100A9", "CLEC9A", "MMP3"),
    SMC=c("MYH11", "MYL9", "MYLK", "TAGLN", "ACTG2"),
    EC=c("VWF", "PECAM1", "CD34", "CD93", "EGFL7", "KLF2"),
    fib1=c("C7", "THBS4", "SLIT2", "PI16", "SFRP1", "SFRP2", "COL14A1"),
    fib2=c(
        "COL1A1", "COL1A2", "COL3A1", "COL5A3", "COL6A1", "COL6A2",
        "CTHRC1", "COMP", "PDGFRA", "FAP", "FN1", "F3", "SULF1")) 

setdiff(unlist(gs), rownames(sce))
unlist(gs)[duplicated(unlist(gs))]

# wrangling
cd <- colData(sce)[c("kid", "sid")]
es <- t(logcounts(sce[unlist(gs), ]))
df <- data.frame(es, cd, check.names=FALSE)
fd <- df |>
    filter(sid == "all") |> 
    select(-sid) |>
    dplyr::rename(k=kid) |>
    pivot_longer(-k, names_to="g", values_to="y") |>
    group_by(g) |> mutate_at("y", .z) |>
    mutate_at("k", factor, rev(names(gs))) |>
    mutate_at("g", factor, unlist(gs))

# plotting
gg <- ggplot(fd, aes(g, k, fill=y)) +
    geom_tile() + coord_equal(4/3, expand=FALSE) +
    scale_fill_gradient2(
        "z-scaled\nmean expr.",
        low="dodgerblue2", high="tomato2",
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
    .thm_fig_c("minimal") + theme(
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5))

# saving
pdf(args[[2]], onefile=TRUE, width=15/2.54, height=4/2.54)
for (p in list(gg)) print(p); dev.off()
