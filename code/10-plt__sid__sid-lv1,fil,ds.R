# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
ist <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# wrangling
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
colnames(xy) <- c("x", "y")
ns <- nn2(xy, 
    k=(k <- 301), 
    r=(r <- 0.05), 
    searchtype="radius")
is <- ns$nn.idx[, -1]
range(rowSums(is > 0))
is[is == 0] <- NA
ks <- (ks <- ist$clust)[match(rownames(xy), names(ks))]
ks <- matrix(ks[c(is)], nrow(is), ncol(is))

# plotting
ps <- lapply(sort(unique(ist$clust)), \(sub) {
    fq <- .z(rowSums(ks == sub, na.rm=TRUE))
    .plt_xy(data.frame(xy), unname(fq), 
        id=paste0(wcs$x, "-", sub)) +
        scale_color_gradientn(
            sprintf("cellular density\n(%sum radius)", 1e3*r),
            colors=pals::jet(), limits=c(-2.5, 2.5), n.breaks=5)
})

# saving
df <- ps[[1]]$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=4+dx/2, height=0.5+dy/2)
