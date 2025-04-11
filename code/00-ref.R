# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(Seurat)
    library(scuttle)
    library(SingleCellExperiment)
})

# loading
dim(obj <- readRDS(args[[1]]))
sce <- as.SingleCellExperiment(obj)

# filtering
idx <- 
    sce$disease == "control" &
    sce$level_1_annot == "Epithelial" &
    sce$organ_groups == "Large_intestine"
sub <- sce[, idx]
table(idx)

# labelling
foo <- list(
    "TA"="TA",
    "tuft"="Tuft",
    "goblet"="Goblet",
    "stem"="Epithelial_stem",
    "enterocyte_BEST4"=c("BEST4_enterocyte_colonocyte"),
    "enterocyte"=c("Colonocyte", "Mature_colonocyte"),
    "enteroendocrine"=c(
        "Enteroendocrine", 
        "Enteroendocrine_progenitor"))
old <- as.character(sub$level_3_annot)
idx <- match(old, unlist(foo))
new <- rep.int(names(foo), sapply(foo, length))[idx]
table(old, new)

# selection
mgs <- findMarkers(sub, groups=new, direction="up")
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 100])
length(sel <- unique(unlist(top)))

# aggregation
sub <- logNormCounts(sub, transform="none")
pbs <- aggregateAcrossCells(sub, new, subset.row=sel,
    use.assay.type="normcounts", statistics="mean")

# saving
saveRDS(assay(pbs), args[[2]])