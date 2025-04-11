# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# setup
set.seed(250410)

# loading
sce <- readRDS(args[[1]])
roi <- readRDS(args[[2]])
ist <- readRDS(args[[3]])

# root in reference-like region
ids <- rep("two", ncol(sce))
names(ids) <- colnames(sce)
ref <- grep("REF", roi$typ)
ref <- colnames(roi)[ref]
ref <- names(ids) %in% ref
ids[ref] <- "one"
table(ids)

kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
table(ids <- ifelse(grepl("EE|entero|goblet", kid), "one", "two"))

# analysis
sce <- slingshot(sce, 
    approx_points=100,
    clusterLabels=ids, 
    reducedDim="PCA", 
    start.clus="one")

# saving
base::saveRDS(sce, args[[4]])
