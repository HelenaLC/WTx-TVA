# args <- list(
#   list.files("outs", "^lv2-", full.names=TRUE),
#   "plts/lv2,mgs.pdf")
# 
# # dependencies
# suppressPackageStartupMessages({
#     library(dplyr)
#     library(tidyr)
#     library(scater)
#     library(scuttle)
#     library(ggplot2)
#     library(SingleCellExperiment)
# })
# 
# # loading
# sid <- gsub(".*-([0-9]+),.*\\.rds", "\\1", args[[1]])
# sub <- gsub(".*(imm|epi|str).*", "\\1", args[[1]])
# 
# foo <- split(seq_along(args[[1]]), sub)
# dfs <- lapply(foo, \(i) {
#     lapply(i, \(j) {
#         y <- readRDS(args[[1]][j])
#         y <- normalizeCounts(y$profiles)
#         data.frame(y, g=rownames(y),
#             sub=sub[j], sid=sid[j],
#             check.names=FALSE)
#     }) |> bind_rows()
# })
# 
# # average across sections
# mu <- lapply(dfs, \(df) df |>
#     group_by(g) |>
#     summarize(across(
#         where(is.numeric),
#         \(.) mean(., na.rm=TRUE))) |>
#     pivot_longer(
#         where(is.numeric),
#         names_to="k", values_to="y"))

# restrict to selected features
gs <- list(
    imm=list(
        Bm=c("BCL7A", "CR2", "CD19", "CD40", "CD69", "CD79A", "CD83", "CIITA", "SELL", "IGHD"),
        Bn=c("BANK1", "MS4A1"),
        Ba=c("BIRC5"),
        PC_IgA=c("MZB1", "XBP1", "IGHA1"),
        PC_IgG=c("IGHG1"), 
        PC_IgM=c("IGHM"),
        Th=c("CD2", "CD3D", "CD4", "ICOS", "CCR4"),
        Tc=c("CD8A", "GZMA", "GZMK", "NKG7"),
        Tp=c("MKI67", "TOP2A"),
        mast=c("KIT", "IL1RL1"),
        cDC2_CCR7=c("CCR7"),
        cDC1=c("LAMP3", "CLEC10A", "ITGAX"),
        nphil=c("IL1B", "S100A8", "S100A9", "MMP9"),
        TRM=c("CD14", "CD68", "CD163", "APOE", "APOC1", "CTSD", "C1QA", "C1QB", "C1QC")
        # TC=c("CD2", "CD3D", "CD8A", "CD28", "TIGIT", "GATA3", "NKG7", 
        #     "EOMES", "GZMK", "CXCR6", "CCL5", "ITGAL", "CD40LG", "CD4"),
    ),
    epi=list(
        epi.EE=c("CHGA", "CHGB", "ENPP2", "SST", "PYY"),
        epi.goblet=c("MUC5AC"),
        epi.entero=c("SLC26A2", "CA1", "CA4", "CA7", "FABP1", "KRT19", "KRT20"),
        epi.stress=c("DUOX2", "DUOXA2", "CTSE", "KLK10"),
        `epi.goblet-progen`=c("REG4"),
        epi.stem=c("OLFM4", "KCNQ1"),
        epi.prolif=c("MKI67", "DNMT1"),
        epi.stem_LGR5.1=c("NOTUM", "AXIN2", "ZBTB4", "MCTP1", "ASCL2", "LGR5"),
        epi.stem_LGR5.2=c("CLDN2", "ODAM", "CCL15", "GHR"),
        epi.4=c("TMEM63C"),
        `epi.paneth-like`=c("REG1A", "REG1B", "REG3A", "DEFA5", "DEFA6"),
        epi.1=c("MMP7", "KRT23", "KRT17", "KRT80"),
        epi.6=c("DMXL2", "TUSC3"),
        epi.2=c("TSLP", "ESPN", "ALPI", "ATP12A"),
        epi.3=c("C3", "CISH"),
        epi.5=c("SPARC", "CALD1")
    ),
    str=list(
        SMC=c("MYL9", "MYLK", "MYH11", "ACTA2", "TAGLN"),
        BEC=c("CDH5", "FLT1", "CLDN5"),
        LEC=c("PECAM1", "VWF", "CD34", "CD93", "EGFL7", "KLF2", "LYVE1", "PROX1", "CCL21"),
        fib=c("PTGIS", "SFRP1", "SFRP2", "PI16"),
        fib_TLS=c("C3", "CLU", "LTA3", "CCL19", "CXCL13", "SPIB"),
        IAF_MMP=c("MMP1", "MMP3", "MMP10", "MMP19", "CXCL5", "CXCL8", "SMOX", "ASRGL1"),
        fib_ADAMDEC1=c("ADAMDEC1", "ADAM28", "BMP5", "FGFR4", "CCL11"),
        CAF_L1CAM=c("L1CAM", "ITGB8", "CCBE1"),
        CAF_CXCL14=c("COL4A5", "CXCL14", "F3", "PDGFRA"),
        CAF_CTHRC1=c("CTHRC1", "FN1", "COL10A1", "FAP"),
        CAF_COL4=c("PDGFA", "COL4", "COL5A3", "LZTS1")
    )
)
ps <- lapply(names(gs), \(.) {
fd <- mu[[.]] |>
    filter(g %in% unlist(gs[[.]])) |>
    mutate_at("g", factor, unlist(gs[[.]])) |>
    mutate_at("k", factor, rev(names(gs[[.]]))) |>
    # z-scale across clusters
    group_by(g) |>
    mutate_at("y", .z)

# plotting
ggplot(fd, aes(g, k, fill=y)) +
    geom_tile() + coord_equal(4/3, expand=FALSE) +
    scale_fill_gradientn(
        "z-scaled\nmean count",
        colors=rev(pals::brewer.rdbu(9)),
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
        axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5))
})

# saving
pdf(args[[2]], width=12/2.54, height=6.5/2.54, onefile=TRUE)
for (p in ps) print(p); dev.off()
