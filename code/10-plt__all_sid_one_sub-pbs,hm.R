args <- list(
    list.files("outs", "pbs.*epi", full.names=TRUE)
)
sce <- lapply(args[[1]], readRDS)

gs <- list(
    imm=list(
        Bm=c("CD79A", "BCL7A", "TCL1A", "CR2"),
        Bn=c("SELL", "CD19", "BANK1", "MS4A1", "FCER2", "FCRLA"),
        Ba=c("IGHD", "IER5", "BIRC5"),
        PC_IgA=c("XBP1", "MZB1", "IGHA1"),
        PC_IgG=c("IGHG1"), 
        PC_IgM=c("IGHM"),
        Th=c("CD2", "CD3D", "ICOS", "CCR4"),
        Tp=c("MKI67", "TOP2A"),
        Tc=c("CD8A", "GZMA", "GZMB", "GZMK", "NKG7"),
        mast=c("KIT", "IL1RL1", "GATA2", "SLC18A2"),
        nphil=c("IL1B", "S100A8", "S100A9"),
        cDC2_CCR7=c("CD40", "CD83", "CCR7", "BATF3", "CD1C", "CIITA"),
        cDC1=c("CLEC9A", "CLEC10A", "FCER1A", "ITGAX"),
        TRM=c("APOE", "APOC1", "CD4", "CD14", "CD68", "CD163", "CTSD", "C1QA", "C1QB", "C1QC", "MMP9")
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
        SMC=c("MYL9", "MYLK", "MYH11", "TAGLN", "ACTA2"),
        LEC=c("PROX1", "LYVE1", "CCL21", "CLDN5"),
        BEC=c("EGFL7", "CDH5", "KLF2", "FLT1", "PECAM1", "VWF", "CD34", "CD93"),
        fib=c("PTGIS", "SFRP1", "SFRP2", "PI16"),
        fib_TLS=c("CLU", "CCL19", "CXCL13"),
        IAF_MMP=c("MMP1", "MMP3", "MMP10", "CXCL5", "CXCL8", "SMOX"),
        fib_ADAMDEC1=c("ADAMDEC1", "ADAM28", "BMP5", "FGFR4", "CCL11"),
        CAF_L1CAM=c("L1CAM", "ITGB8"),
        CAF_CXCL14=c("CXCL14", "F3", "PDGFRA", "COL4A5", "COL4A6"),
        CAF_CTHRC1=c("CTHRC1", "FN1", "COL10A1", "FAP"),
        CAF_COL4=c("PDGFA", "COL5A3", "LZTS1", "COL4A1", "COL4A2")
    )
)

gs <- gs$epi
ks <- rev(names(gs))
gs <- unique(unlist(gs))

df <- lapply(sce, \(.) {
    es <- assay(.)
    data.frame(g=rownames(es), es, check.names=FALSE)
}) |>
    bind_rows(.id="s") |>
    filter(g %in% gs) |>
    group_by(g) |>
    summarise(across(-s, mean, na.rm=TRUE)) |>
    pivot_longer(-g, names_to="k", values_to="y") |>
    group_by(g) |>
    mutate_at("y", .z) |>
    mutate(g=factor(g, gs), k=factor(k, ks))

ggplot(df, aes(g, k, fill=y)) +
    coord_equal(expand=FALSE) +
    geom_tile() +
    scale_fill_gradientn(
        "z-scaled\nmean cpa",
        limits=c(-2.5, 2.5), n.breaks=5,
        colors=rev(pals::brewer.rdbu(9))) +
    .thm_fig_c("minimal") + theme(
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

w <- length(unique(df$g))+max(nchar(paste(df$g)))/2
h <- length(unique(df$k))+max(nchar(paste(df$k)))/2
ggsave("plts/pbs,hm,imm.pdf", units="cm", width=1+w/5, height=h/5)