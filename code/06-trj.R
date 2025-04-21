# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(slingshot)
    library(SingleCellExperiment)
})

# setup
set.seed(250414)

# loading
sce <- readRDS(args[[1]])
roi <- readRDS(args[[2]])

# root in reference-like;
# end in most malignant
ids <- rep("foo", ncol(sce))
typ <- gsub("^.*_", "", roi$typ)
names(ids) <- colnames(sce)
names(typ) <- colnames(roi)
typ <- typ[names(ids)]
fin <- rev(names(.pal_roj))
fin <- fin[fin %in% typ][1]
fin <- names(which(typ == fin))
ref <- names(which(typ == "REF"))
ids[ref] <- "one"; ids[fin] <- "two"
table(ids); table(is.na(ids))

# analysis
sce <- slingshot(sce,
    clusterLabels=ids, 
    approx_points=100,
    reducedDim="PCA", 
    start.clus="one",
    end.clus="two")

# saving
base::saveRDS(sce, args[[3]])
