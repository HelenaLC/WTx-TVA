# loading
sce <- readRDS(args[[1]])
library(SingleCellExperiment)

# aesthetics
aes <- scale_color_gradientn(
    limits=c(-2.5, 2.5), n.breaks=6,
    "z-scaled\nCNV score", colors=pals::jet())

# joint
sco <- colMeans(abs(assay(sce)))
p <- .plt_xy(sce, .z(unname(sco)), wcs$x) + aes

# split
chr <- rowData(sce)$chr
q <- lapply(unique(chr), \(.) {
    bin <- rownames(sce)[chr == .]
    sco <- colMeans(abs(assay(sce[bin, ])))
    .plt_xy(sce, .z(unname(sco))) +
    ggtitle(.lab(., length(bin))) + aes
})

# saving
df <- p$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=c(list(p), q), nrow=1, ncol=1, top=FALSE)
ggsave(args[[2]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
