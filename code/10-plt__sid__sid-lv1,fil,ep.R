# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(ggrastr)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
ist <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# wrangling
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
colnames(xy) <- c("x", "y")

j <- names(k <- ist$clust)
i <- match(colnames(sce), j)
xy <- xy[grep("^epi$", k[i]), ]

# analysis
ns <- nn2(xy, 
    r=(r <- 0.05), k=301, 
    searchtype="radius")

is <- ns$nn.idx[,-1]
is[is == 0] <- NA
range(rowSums(!is.na(is)))

as <- matrix(sce$Area.um2[c(is)], 
    nrow=nrow(is), ncol=ncol(is)) |>
    rowSums(na.rm=TRUE)
df <- data.frame(xy, z=as/(pi*((1e3*r)**2)))

# aesthetics
dx <- diff(range(df$x))
dy <- diff(range(df$y))
pt <- min(dx, dy)/1e3/2
nc <- format(nrow(df), big.mark=",")

# plotting
gg <- ggplot(
    df[order(df$z), ], aes(x, y, col=.z(z))) +
    geom_point_rast(shape=16, stroke=0, size=pt) +
    scale_color_gradient2(
        sprintf("epi. coverage\n(%sum radius)", 1e3*r), 
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2),
        low="royalblue3", mid="grey90", high="tomato3") +
    ggtitle(bquote(bold(.(wcs$sid))~"(N ="~.(nc)*")")) +
    coord_equal(expand=FALSE) + 
    theme_void(6) + theme(
        legend.position="bottom",
        legend.key=element_blank(),
        plot.background=element_blank(),
        legend.background=element_blank(),
        legend.key.width=unit(0.4, "lines"),
        legend.key.height=unit(0.2, "lines"),
        plot.title=element_text(hjust=0.5))

# saving
ggsave(args[[3]], gg, units="cm", width=dx/2, height=1+dy/2)
