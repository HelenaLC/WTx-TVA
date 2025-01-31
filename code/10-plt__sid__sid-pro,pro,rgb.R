# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(matrixStats)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])

# wrangling
pcs <- reducedDim(sce, "PCA")[, seq_len(3)]
pcs <- sweep(pcs, 1, rowMins(pcs), `-`)
pcs <- sweep(pcs, 1, rowMaxs(pcs), `/`)
df <- data.frame(colData(sce), pcs, check.names=FALSE)

# aesthetics
xy <- grep("global_mm$", names(df))
names(df)[xy] <- c("x", "y")
dx <- diff(range(df$x))
dy <- diff(range(df$y))
pt <- min(dx, dy)/100/4

# plotting
gg <- ggplot(df, aes(x, y, col=rgb(PC1, PC2, PC3))) + 
    geom_point_rast(shape=16, stroke=0, size=pt, raster.dpi=600) +
    ggtitle(.lab(wcs$sid, ncol(sce))) +
    scale_color_identity() + .thm_xy

# saving
ggsave(args[[3]], gg, units="cm", width=dx/2, height=0.5+dy/2)
