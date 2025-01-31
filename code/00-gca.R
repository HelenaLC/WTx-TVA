# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(scuttle)
    library(zellkonverter)
    library(SingleCellExperiment)
})

# retrieve data from
# https://www.gutcellatlas.org
# pan-GI cell atlas: healthy reference
options(timeout=600)
url <- "https://datasets.cellxgene.cziscience.com/97eaf591-644f-4c5b-ad31-ee730ce1e9e9.h5ad"
file.create(tf <- tempfile(fileext=".h5ad"))
download.file(url, tf, quiet=TRUE)
(sce <- readH5AD(tf))

# subset to cells of interest
epi <- sce[, sce$category == "Epithelial"]
ns <- table(droplevels(epi$cell_type), epi$tissue)
ns[, ts <- c("colon", "small intestine")]
sub <- epi[, epi$tissue %in% ts]
table(droplevels(sub$cell_type))

# simplify annotations
ee <- grep("enteroendocrine", 
    unique(sub$cell_type), 
    value=TRUE)
(lab <- c(ee,
    "stem cell"="stem",
    "progenitor cell"="progen",
    "epithelial cell"="BEST4+",
    "transit amplifying cell"="TA",
    "intestine goblet cell"="goblet",
    "enterocyte of colon"="enterocyte",
    "colon epithelial cell"="enterocyte",
    setNames(rep("EE", length(ee)), ee)))
table(sub$kid <- lab[match(sub$cell_type, names(lab))])

# (log-)normalization
assayNames(sub) <- "counts"
sub <- logNormCounts(sub, log=FALSE)
sub <- logNormCounts(sub, log=TRUE)

# feature selection
pro <- (rd <- rowData(sub))$feature_type == "protein_coding"
mgs <- findMarkers(sub[pro, ], sub$kid, add.summary=TRUE)
sel <- lapply(mgs, \(df) rownames(df)[df$Top <= 200])
sapply(sel, length); length(sel <- unique(unlist(sel)))

# reference profiles
pbs <- aggregateAcrossCells(sub[sel, ], ids=sub$kid)
sizeFactors(pbs) <- NULL; pbs <- logNormCounts(pbs, log=FALSE)
idx <- match(rownames(pbs), rownames(sub))
rownames(pbs) <- rowData(sub)$feature_name[idx]

# saving
rowData(sub)$sel <- rownames(sub) %in% sel
saveRDS(sub, "data/gca_sc.rds")
saveRDS(pbs, "data/gca_pb.rds")