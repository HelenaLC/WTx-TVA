# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrastr)
    library(RColorBrewer)
    library(SingleCellExperiment)
})

# loading
ist <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# wrangling
xy <- "Center(X|Y)_global_mm"
xy <- grep(xy, names(colData(sce)))
names(colData(sce))[xy] <- c("x", "y")
cs <- match(colnames(sce), names(k <- ist$clust))
nk <- nlevels(k <- factor(k[cs]))
df <- data.frame(colData(sce), k)
nc <- format(ncol(sce), big.mark=",")
dx <- diff(range(df$x))
dy <- diff(range(df$y))
pt <- min(dx, dy)/100

# aesthetics
pal <- if (nk == 3) {
    c("gold", "cyan", "magenta")
} else if (nk > 12) {
    colorRampPalette(.pal)(nk)
} else .pal

# plotting
ggplot(df, aes(x, y, col=k)) +
    scale_color_manual(NULL, limits=levels(df$k), values=pal, drop=FALSE) +
    ggtitle(bquote(bold(.(wcs$sid))~"(N ="~.(nc)*")")) +
    guides(col=guide_legend(override.aes=list(size=2))) +
    geom_point_rast(shape=16, stroke=0, size=pt, raster.dpi=600) +
    coord_equal(expand=FALSE) + 
    theme_void(6) + theme(
        legend.key=element_blank(),
        plot.background=element_blank(),
        legend.background=element_blank(),
        legend.key.size=unit(0.5, "lines"),
        plot.title=element_text(hjust=0.5))

ps <- lapply(c(levels(df$k), NA), \(k) { 
    df$. <- grepl(sprintf("^%s$", k), df$k)
    nc <- format(sum(df$.), big.mark=",")
    df <- df[order(df$.), ]
    ggplot(df, aes(x, y, col=., size=.)) +
        geom_point_rast(shape=16, stroke=0, raster.dpi=600) +
        scale_color_manual(values=c("lavender", "blue")) +
        ggtitle(bquote(bold(.(k))~"(N ="~.(nc)*")")) +
        scale_size_manual(values=c(pt/2, pt)) +
        coord_equal(expand=FALSE) + 
        theme_void(6) + theme(
            legend.position="none",
            plot.background=element_blank(),
            plot.title=element_text(hjust=0.5))
})

# saving
pdf(args[[3]], onefile=TRUE, 
    width=(2+dx/2)/2.54, 
    height=(0.5+dy/2)/2.54)
for (p in c(list(p0), ps)) print(p); dev.off()
