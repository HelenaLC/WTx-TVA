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
    "enteroendocrine"=c("CHGA", "CHGB", "SCGN", "PYY", "SST"),
    "enterocyte_BEST4"=c("BEST4", "OTOP2", "CA7"),
    "enterocyte"=c("CA4", "CA1", "SLC26A3", "CEACAM7", "GPT", "PI3", "FABP2", "AQP3"),
    "goblet"=c("MUC2", "TFF3", "CLCA1", "ITLN1", "SPINK4"),
    "tuft"=c("HPGDS", "PLCG2", "DAPK1"),
    "stem"=c("LGR5", "OLFM4", "ASCL2", "RGMB", "TSPAN8", "DMBT1"),
    "TA"=c("PCNA", "MKI67", "TOP2A", "PCLAF"))

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
pdf(args[[2]], onefile=TRUE, width=12/2.54, height=2/2.54)
for (p in list(gg)) print(p); dev.off()
