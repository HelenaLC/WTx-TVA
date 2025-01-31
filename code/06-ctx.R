# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(imcRtools)
    library(BiocParallel)
    library(SingleCellExperiment)
})

# setup
th <- as.integer(args$ths)
bp <- MulticoreParam(th)
set.seed(241202)

# loading
sce <- mapply(x=args[[1]], y=args[[2]], \(x, y) {
    sce <- readRDS(x); ist <- readRDS(y)
    sce$kid <- (k <- ist$clust)[match(colnames(sce), names(k))]
    assays(sce) <- altExps(sce) <- list(); `rowData<-`(sce, value=NULL)
}, SIMPLIFY=FALSE)

# restrict to shared features & collapse
gs <- Reduce(intersect, lapply(sce, rownames))
sce <- do.call(cbind, lapply(sce, \(.) .[gs, ]))
table(sce$sid); dim(sce)

# wrangling
xy <- "Center(X|Y)_global_mm"
xy <- grep(xy, names(colData(sce)))
names(colData(sce))[xy] <- c("x", "y")

# build cellular KNN graph based on 
# Euclidean distances w/ thresholding
sce <- buildSpatialGraph(sce,
    img_id="sid", type="expansion", threshold=0.1, 
    coords=c("x", "y"),  name="knn", BPPARAM=bp)

# get subpopulation frequencies in each neighborhood
sce <-  aggregateNeighbors(sce, 
    colPairName="knn", name="fqs",
    aggregate_by="metadata", count_by="kid")
colPairs(sce) <- list() # free up some memo

# do k-means clustering on neighborhood compositions
fqs <- colData(sce)[["fqs"]]
rownames(fqs) <- colnames(sce)
fqs <- fqs[!is.na(fqs[, 1]), ]
ctx <- kmeans(fqs, centers=30)$clust
idx <- match(colnames(sce), rownames(fqs))
table(sce$ctx <- factor(ctx[idx], sort(unique(ctx))))

# saving
saveRDS(sce, args[[3]])
