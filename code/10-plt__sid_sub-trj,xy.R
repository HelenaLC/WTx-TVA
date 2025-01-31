# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])

# wrangling
ts <- as.matrix(colData(sce)[grep("time", names(colData(sce)))])
cd <- colData(sce)[-grep("sling", names(colData(sce)))]
df <- data.frame(cd, t=.q(rowMeans(ts, na.rm=TRUE)))

# aesthetics
nc <- format(nrow(df), big.mark=",")
xy <- "Center(X|Y)_global_mm"
xy <- grep(xy, names(df))
names(df)[xy] <- c("x", "y")
dx <- diff(range(df$x))
dy <- diff(range(df$y))
pt <- min(dx, dy)/100

# plotting
gg <- ggplot(df[order(df$t), ], aes(x, y, col=t)) +
    ggtitle(bquote(bold(.(wcs$sid))~"(N ="~.(nc)*")")) +
    geom_point_rast(shape=16, stroke=0, size=pt, raster.dpi=600) +
    scale_color_gradientn(NULL, 
        colors=rev(hcl.colors(9, "Rocket")),
        breaks=c(0, 1), labels=c("early", "late")) +
    coord_equal(expand=FALSE) + 
    theme_void(6) + theme(
        legend.position="bottom",
        legend.key=element_blank(),
        plot.background=element_blank(),
        legend.background=element_blank(),
        plot.title=element_text(hjust=0.5),
        legend.key.width=unit(0.4, "lines"),
        legend.key.height=unit(0.2, "lines"))

# saving
ggsave(args[[2]], gg, units="cm", width=dx/2, height=1+dy/2)
