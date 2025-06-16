# args <- list(
#     "outs/raw-232", "outs/fil-232.rds",
#     "outs/roi-232.rds", "outs/pro-232.rds",
#     "outs/ist-232.rds", "outs/lv1-232.rds",
#     "outs/sig-232.rds", "outs/ccc-232.rds",
#     "outs/cty-232.rds", "outs/trj-232.rds",
#     list.files("outs", "rep-232", full.names=TRUE),
#     list.files("outs", "jst-232", full.names=TRUE),
#     list.files("outs", "lv2-232", full.names=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(Matrix)  
    library(AUCell)  
    library(HDF5Array)  
    library(slingshot)  
    library(alabaster.sce)  
    library(zellkonverter)  
    library(SingleCellExperiment)  
})

# initialize output directory
args <- unlist(args)
dir <- tail(args, 1)
dir.create(dir)

# load up original data
raw <- grep("raw", args, value=TRUE)
sce <- loadHDF5SummarizedExperiment(raw)
assay(sce) <- as(assay(sce), "dgCMatrix")

# quality control
fil <- readRDS(grep("fil", args, value=TRUE))
metadata(sce) <- metadata(fil)
sce$fil <- colnames(sce) %in% colnames(fil)

# domains & regions of interest
roi <- readRDS(grep("roi", args, value=TRUE))
idx <- match(colnames(sce), colnames(roi))
cd <- colData(roi)[idx, c("typ", "roi")]
colData(sce) <- cbind(colData(sce), cd)

# signature scores
sig <- readRDS(grep("sig", args, value=TRUE))
idx <- match(colnames(sce), colnames(sig))
sig <- `colnames<-`(sig[, idx], colnames(sce))
assay(sig) <- as(assay(sig), "dgCMatrix")
altExp(sce, "AUCell") <- sig

# cell-cell communication
ccc <- readRDS(grep("ccc", args, value=TRUE))
names(ccc) <- c("sender", "receiver")
ccc <- lapply(ccc, \(.) {
    mtx <- as.matrix(.)
    mty <- t(mtx[!rowAlls(as.matrix(is.na(mtx))), ])
    idx <- match(colnames(sce), colnames(mty))
    mty <- as(mty[, idx], "dgCMatrix")
    `colnames<-`(mty, colnames(sce))
}) 
ccc <- SingleCellExperiment(ccc)
altExp(sce, "COMMOT") <- ccc

# spatial contexts
ctx <- readRDS(grep("cty", args, value=TRUE))
sce$ctx <- ctx$ctx[match(colnames(sce), ctx$cid)]
table(sce$ctx)

# epithelial trajectory
trj <- readRDS(grep("trj", args, value=TRUE))
idx <- match(colnames(sce), colnames(trj))
sce$trj <- slingAvgPseudotime(trj$slingshot)[idx]
table(is.na(sce$trj[sce$lv1 != "epi"]))

# principal components
pat <- ".*(pro|rep).*"
drs <- c(pro="PCA.fil", rep="PCA.epi")
for (rds in grep(pat, unlist(args), value=TRUE)) {
    pcs <- reducedDim(readRDS(rds), "PCA")
    mtx <- matrix(NA, ncol(sce), ncol(pcs))
    rownames(mtx) <- colnames(sce)
    colnames(mtx) <- colnames(pcs)
    mtx[rownames(pcs), ] <- pcs
    typ <- gsub(pat, "\\1", rds)
    reducedDim(sce, drs[typ]) <- mtx
}

# clustering & annotations
pat <- ".*(ist|jst|lv1|lv2).*"
ids <- \(.) paste(sort(unique(.)))
for (rds in grep(pat, unlist(args), value=TRUE)) {
    . <- gsub(pat, "\\1", rds)
    ist <- readRDS(rds)
    mtx <- ist$profiles
    kid <- ist$clust
    if (. %in% c("ist", "jst")) {
        # feature selection
        sel <- rownames(sce) %in% rownames(mtx)
        sub <- gsub(".*(epi|imm|str).*", "\\1", rds)
        idx <- ifelse(grepl("rds", sub), "all", sub)
        rowData(sce)[[idx]] <- sel
    }
    if (is.null(sce[[.]])) {
        # new entry
        idx <- match(colnames(sce), names(kid))
        sce[[.]] <- setNames(kid[idx], colnames(sce))
    } else {
        old <- ids(sce[[.]])
        # existing entry
        sce[[.]] <- setNames(paste(sce[[.]]), colnames(sce))
        sce[[.]][names(kid)] <- paste(kid)
        sce[[.]] <- factor(sce[[.]], c(old, ids(kid)))
    }
}
table(sce$ist, sce$lv1)
table(sce$lv2, sce$lv1)
sel <- c("all", "epi", "imm", "str")
sapply(sel, \(.) table(rowData(sce)[[.]]))

# saving
saveRDS(sce, file.path(dir, paste0(wcs$sid, ".rds")))
alabaster.sce::saveObject(sce, file.path(dir, wcs$sid))
zellkonverter::writeH5AD(sce, compression="gzip", file.path(dir, paste0(wcs$sid, ".h5ad")))
