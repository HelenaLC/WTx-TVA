# dependencies
suppressPackageStartupMessages({
    library(Matrix)
    library(scater)
    library(reticulate)
    library(BiocParallel)
    library(zellkonverter)
    library(SingleCellExperiment)
})

# setup
bin <- "/software/crgadm/software/Miniconda3/4.9.2/bin"
use_condaenv("commot", file.path(bin, "conda"), FALSE)
use_python(file.path(bin, "python"))
th <- as.integer(args$ths)
bp <- MulticoreParam(th)

# loading
sce <- readRDS(args[[1]])

# normalize by area
cpa <- normalizeCounts(sce,
    size.factors=sce$Area.um2,
    center.size.factors=FALSE)
assay(sce, "cpa") <- cpa

# wrangling
xy <- grep("global_mm$", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
reducedDim(sce, "spatial") <- xy

# run 'COMMOT'
ad <- import("anndata")
ct <- import("commot")
pd <- import("pandas")

# retrieve LR interaction from CellChatDB & filter for matches
db <- ct$pp$ligand_receptor_database(database="CellChat", species="human")
names(db) <- c("ligand", "receptor", "pathway", "type"); nrow(db)
db <- db[apply(db, 1, \(.) {
    rs <- strsplit(.["receptor"], "_")
    lr <- c(.["ligand"], unlist(rs))
    all(lr %in% rownames(sce))
}), ]; nrow(db)

# subset to features being considered
rs <- sapply(strsplit(db$receptor, "_"), .subset, 1)
nrow(sce <- sce[unique(c(db$ligand, rs)), ])

is <- split(colnames(sce), sce$fov)
sr <- bplapply(is, BPPARAM=bp, \(.) { 
    tryCatch(error=\(e) e, {
        # skip FOV when there are too few cells
        if (length(.) < 200) return(NULL) 
        ad <- SCE2AnnData(sce[, .], X_name="cpa")
        ct$tl$spatial_communication(ad,
            database_name="CellChatDB",
            # average cell is around 10x10um; here, we set
            # a distance threshold of 2 cells (20um=.02mm)
            dis_thr=0.02,  
            df_ligrec=db,
            heteromeric=TRUE,
            pathway_sum=TRUE,
            heteromeric_rule="min",
            heteromeric_delimiter="_")
        list(
            ad$obsm["commot-CellChatDB-sum-sender"], 
            ad$obsm["commot-CellChatDB-sum-receiver"])
    })
})
er <- sapply(sr, inherits, "error")
ex <- sapply(sr, is.null)
table(na <- er | ex)
sr <- sr[!na]

# wrangling
i <- colnames(sce)
s <- lapply(sr, \(.) .[[1]])
r <- lapply(sr, \(.) .[[2]])
s <- do.call(rbind, unname(s))
r <- do.call(rbind, unname(r))
res <- list(s=s[i, ], r=r[i, ])
res <- lapply(res, \(.) {
    mtx <- as(as.matrix(.), "dgCMatrix")
    nms <- gsub("^(s|r)-", "", colnames(.))
    `colnames<-`(mtx, nms)
})

# saving
saveRDS(res, args[[2]])
