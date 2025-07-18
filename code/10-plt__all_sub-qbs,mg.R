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
sub <- gsub(".*(epi|imm|str).*", "\\1", args[[1]])
sce <- lapply(setNames(args[[1]], sub), readRDS)

# restrict to selected features
mgs <- list(
    epi=list(
        "epi.EE"=c("PYY", "SST", "CHGA", "CHGB", "NEUROD1", "GCG", "INSL5"),
        "epi.goblet"=c("BEST2", "MUC2", "MUC4", "TFF3", "FCGBP", "ZG16", "ATOH1", "KLK1", "CLCA1", "ITLN1"),
        "epi.goblet-like_REG4"=c("MUC5AC", "REG4", "GRIP2", "RETNLB", "SYTL1", "PLXND1"),
        "epi.inflam_REG1A"=c("SAA1", "OLFM4", "DMBT1", "S100A9", "HLA-DRA", "HLA-DPA1", "CXCL1", "CXCL2", "CXCL3"),
        "epi.paneth-like_DEFA6"=c("REG1A", "REG1B", "REG3A", "DEFA5", "DEFA6", "DLL1", "ITLN2"),
        "epi.entero"=c("FABP1", "CA1", "CA2", "CA4", "GUCA2A", "MS4A12", "CEACAM7", "SLC26A3"),
        "epi.entero-like_ACSF2"=c("TTLL3", "PIGZ", "CCDC183", "SELENBP1", "ACSF2", "SHROOM1", "EXD3"),
        "epi.entero-like_EREG"=c(
            "EREG", "ANO9", "NCOA7", "TNNC2", "ARID3A", 
            "FITM2", "UCKL1", "SRPX2", "SPAG4", "THBS1", "HBEGF"),
        "epi.stem-like_LGR5"=c("LGR5", "FN1", "STRAP", "CDX2", "TGFBI", "CEL", "SMOC2", "ASCL2", "RNF43", "RGMB"),
        "epi.stem-like_CLU"=c("LEFTY1", "CLU", "PTGDR", "CD44", "SAMD5", "KIF13A", "TRPM5"),
        "epi.fetal-like_MMP7"=c(
            "NOTUM", "KRT23", "UBD", "FLNA", "FSTL3", "CAVIN1", "BST2", 
            "KLK6", "KLK10", "AKAP12", "TIMP2", "PPL", "S100A4", "SDC4", 
            "MMP7", "LAMA3", "LAMB3", "LAMC2", "PLAUR", "SLC2A1", "TNFAIP2", "RHOF"),
        "epi.fetal-like_EMP1"=c(
            "CXCL8", "MSLN", "DUSP5", "TM4SF1", "EMP1", "MUC17",
            "ANXA1", "DUOX2", "DUOXA2", "CD55", "ERO1A", "ECM1", "ENO2", 
            "PFKFB3", "EGLN3", "SERPINE2", "ERRFI1", "CTSE", "COL17A1"),
        "epi.prolif_TACC3"=c(
            "MKI67", "TOP2A", "MYO19", "NCAPG2", 
            "IQGAP3", "FANCA", "RAD54L", "TROAP", "TACC3",
            "STMN1", "SUMO3", "TUBA1B", "TUBB4B", 
            "H2AC4", "H1-3", "H1-5", "H3C7", "H2BC15", "CDK1",
            "RARRES1"),
        "epi.prolif_STMN1"=c(),
        "epi.prolif_H1"=c(),
        "epi.prolif"=c()
    ),
    imm=list(
        B.activated=c("MKI67", "TOP2A", "TUBA1B", "BCL7A", "CR2", "CD19", "CD69", "CD79A"),
        B.naive=c("BANK1", "MS4A1", "CD27", "CD38", "SELL"),
        PC_IgA1=c("MZB1", "XBP1", "SDC1", "JCHAIN"),
        PC_IgM=c("IGHA1", "IGHM", "IGHG1", "IGHD"),
        PC_IgG=c(), 
        PC_IgD=c(),
        Th=c("ITGAL", "TIGIT", "CD2", "CD3D", "CD28", "GATA3", "ICOS", "CCR4", "IL7R"),
        Tc=c("CD8A", "GZMA", "GZMK", "NKG7", "CXCR6", "CCL5"),
        mast=c("KIT", "IL1RL1", "HDC", "GATA2", "SLC18A2"),
        DC.lymphoid=c("CD40", "CD83", "CCR7", "IDO1", "FSCN1", "ID2", "IRF4", "IRF8", "LAMP3"),
        cDC=c("CLEC10A", "CD1C", "ITGAM", "CX3CR1", "ITGAX"),
        neutrophil=c("IL1B", "S100A8", "S100A9", "MMP9"),
        TAM=c("HLA-DPA1", "HLA-DRB1", "CD4", "C1QA", "C1QB", "C1QC", "CD14", "CD68", "APOE", "APOC1", "CTSD", "SPP1", "CCL2", "MARCO"),
        PVM=c("CD163", "CD209", "MS4A4E", "MS4A7", "STAB1", "SELENOP", "SIGLEC1", "MRC1", "LYVE1", "CCL18")
    ),
    str=list(
        pericytes=c("PDGFA", "COL5A3", "LZTS1", "NOTCH3", "RGS5", "CAV1", "GJC1", "TBX2", "ACAN", "ANO1", "NDUFA4L2", "COL4A2","COL18A1"),
        "SMC/fib.myo"=c("ACTA2", "MYL9", "MYLK", "MYH11", "TAGLN", "DES", "ACTG2", "LMOD1", "FHL1", "GREM1"),
        BEC=c(
            "CLDN5", "CDH5", "FLT1", "PECAM1", "VWF", "CD34", "CD93", "EGFL7", "KLF2", 
            "ESM1", "KDR", "RAMP3", "PODXL", "SLCO2A1", "NOTCH4", "SHANK3", "PLVAP"),
        LEC=c("LYVE1", "PROX1", "CCL21", "NTAN1", "RELN"),
        fib.TLS=c("CLU", "CCL19", "CXCL13", "LTB", "LTF"),
        CAF=c(),
        fib_C7=c("C7", "PDGFRL", "THBS4", "PID1", "SLIT2", "BOC", "KCNN3"),
        "fib_JUN/FOS"=c("NR4A1", "EGR1", "CCN1", "JUN", "FOS", "FOSB", "DUSP1"),
        fib_ADAMDEC1=c("ADAMDEC1", "ADAM28", "ABCA8", "HAAO", "BMP5", "FGFR4", "SCARA5"),
        fib_PI16=c("OGN", "PTGIS", "SFRP1", "SFRP2", "CFD", "GSN", "PI16", "MFAP5", "IGFBP6", "PLA2G2A"),
        IAF_MMP=c("MMP1", "MMP3", "MMP10","MMP19", "CXCL1", "CXCL5", "CXCL8", "SMOX", "AREG", "IL11", "HIF1A"),
        CAF_F3=c("PDGFRA", "COL4A5", "CXCL14", "F3", "SOX6", "WNT4", "WNT5A"),
        CAF.myo_FAP=c("MMP2",  "MMP11", "MMP14", "CTHRC1", "FN1", "FAP", "COMP", "SFRP4", "SULF1", "SULF2", "CDH11", "COL1A1", "COL1A2", "COL10A1", "COL11A1", "COL12A1", "LOXL2")
    )
)

lapply(mgs, \(.) unlist(.)[duplicated(unlist(.))])

# aesthetics
aes <- list(
    geom_tile(), 
    coord_equal(4/3, expand=FALSE),
    .thm_fig_c("minimal"), theme(
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5)),
    scale_fill_gradient2(
        "z-scaled\nmean expr.", 
        low="dodgerblue2", high="tomato2",
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)))

ps <- lapply(sub, \(.) {
    se <- sce[[.]]
    gs <- mgs[[.]]
    ns <- cumsum(sapply(gs, length))
    ns <- unique(ns); ns <- ns[-length(ns)]
    # joint
    x <- se[, se$sid == "all"]
    y <- t(logcounts(x[unlist(gs), ]))
    df <- data.frame(y, k=x$kid, check.names=FALSE)
    fd <- df |>
        pivot_longer(-k, names_to="g", values_to="y") |>
        group_by(g) |> mutate_at("y", .z) |>
        mutate_at("k", factor, rev(names(gs))) |>
        mutate_at("g", factor, unlist(gs))
    p <- ggplot(fd, aes(g, k, fill=y)) + aes +
        geom_vline(xintercept=ns+0.5, linewidth=0.2) +
        (if (. == "epi") scale_y_discrete(labels=\(.) gsub("^epi\\.", "", .)))
    # split
    x <- se[, se$sid != "all"]
    y <- t(logcounts(x[unlist(gs), ]))
    df <- data.frame(y, s=x$sid, k=x$kid, check.names=FALSE)
    df <- mutate(df, across(where(is.numeric), .z))
    q <- by(df, df$k, \(fd) {
        k <- fd$k[1]
        fd <- fd |>
            select(-k) |>
            pivot_longer(-s, names_to="g", values_to="y") |>
            mutate_at("g", factor, unlist(gs))
        # plotting
        hl <- c("plain", "bold")[1 + fd$g %in% gs[[k]]]
        ggplot(fd, aes(g, s, fill=y)) + aes + ggtitle(k) +
            geom_vline(xintercept=ns+0.5, linewidth=0.2) +
            scale_y_discrete(limits=\(.) rev(.)) + guides(fill="none") +
            theme(axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5, face=hl))
    })
    c(list(p), q)
}) |> Reduce(f=c)

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- (p <- ps[[.]])$data
    y <- ifelse(is.null(df$k), "s", "k")
    w <- length(unique(df$g))/10+max(nchar(paste(df[[y]])))/5
    h <- length(unique(df[[y]]))/10+max(nchar(paste(df$g)))/5
    pdf(tf[[.]], width=(2+w)/2.54, height=(0.5+h)/2.54)
    print(p); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[2]])
