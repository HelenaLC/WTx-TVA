# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(matrixStats)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[2]], readRDS)

# wrangling
pat <- ".*(epi|imm|str).*"
nms <- basename(args[[2]])
names(sce) <- gsub(pat, "\\1", nms)

df <- lapply(names(sce), \(sub) {
    cd <- colData(obj <- sce[[sub]])
    pcs <- reducedDim(obj, "PCA")[, seq_len(3)]
    pcs <- sweep(pcs, 1, rowMins(pcs), `-`)
    pcs <- sweep(pcs, 1, rowMaxs(pcs), `/`)
    df <- data.frame(cd, pcs, sub, check.names=FALSE)
    xy <- grep("global_mm$", names(df))
    names(df)[xy] <- c("x", "y")
    return(df)
})
    
# aesthetics
ds <- vapply(df, \(.) c(
    diff(range(.$x)),
    diff(range(.$y))),
    numeric(2))
dx <- max(ds[1, ])
dy <- max(ds[2, ])
pt <- min(ds)/100/4
    
# plotting
ps <- lapply(df, \(fd) ggplot(fd, aes(x, y, col=rgb(PC1, PC2, PC3))) + 
    geom_point_rast(shape=16, stroke=0, size=pt, raster.dpi=600) +
    ggtitle(.lab(paste0(wcs$sid, ": ", fd$sub[1]), nrow(fd))) +
    scale_color_identity() + .thm_xy)

# saving
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=0.5+dx/2, height=0.5+dy/2)
