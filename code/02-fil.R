# dependencies
suppressPackageStartupMessages({
    library(scuttle)
    library(HDF5Array)
    library(matrixStats)
    library(SingleCellExperiment)
})

# loading
sce <- loadHDF5SummarizedExperiment(args[[1]])
pos <- read.csv(args[[2]])

# FOVs are 4256 x 4256 px; position
# coordinates  are top-left corners
w <- h <- 4256
x <- sce$CenterX_global_px
y <- sce$CenterY_global_px
i <- match(sce$fov, pos$FOV)
d <- cbind(
    l=x-pos$x_global_px[i], 
    t=pos$y_global_px[i]-y,
    r=pos$x_global_px[i]+w-x,
    b=y-pos$y_global_px[i]+h) 
ol1 <- rowAnys(d < (th_d <- 30))

# low counts overall
x <- sce$nCount_RNA
ol2 <- isOutlier(x, 3, "lower", TRUE)
th_n <- attr(ol2, "threshold")[1]

# low counts per area
x <- sce$nCount_RNA/sce$Area
ol3 <- isOutlier(x, 3, "lower", TRUE)
th_m <- attr(ol3, "threshold")[1]

# high background
x <- sce$nCount_negprobes
ol4 <- x > (th_np <- 10*mean(x))
x <- sce$nCount_falsecode
ol5 <- x > (th_fc <- 10*mean(x))

# filtering
(metadata(sce) <- list(
    th_np=th_np, th_fc=th_fc,
    th_d=th_d, th_n=th_n, th_m=th_m))
ol <- ol1 | ol2 | ol3 | ol4 | ol5
mean(ol1); mean(ol2); mean(ol3)
mean(ol4); mean(ol5)
sub <- sce[, !ol]

# summary
length(unique(sce$fov))
ncol(sce); ncol(sub); mean(ol)
diff(range(sce$CenterX_global_mm))
diff(range(sce$CenterY_global_mm))

# saving
base::saveRDS(sub, args[[3]])