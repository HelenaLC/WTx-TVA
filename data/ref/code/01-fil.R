# dependencies
suppressPackageStartupMessages({
    library(scuttle)
})

# loading
set.seed(240923)
sce <- readRDS(args[[1]])

# quality control
mt <- grep("^MT-", rownames(sce))
cd <- perCellQCMetrics(sce, subsets=list(mt=mt))
colData(sce) <- cbind(colData(sce), cd)

# drop cells with few detected features
ex1 <- isOutlier(cd$detected, log=TRUE, type="lower", nmads=2)
attr(ex1, "thresholds"); mean(ex1)

# drop cells with high mitochondrial content
ex2 <- isOutlier(cd$subsets_mt_percent, type="higher", nmads=3)
attr(ex2, "thresholds"); mean(ex2)

# keep features detected in at least 20 cells
table(gs <- rowSums(counts(sce) > 0) >= 20)

# subsetting
sub <- sce[gs, !(ex1 | ex2)]
tbl <- cbind(
    raw=raw <- dim(sce), 
    fil=fil <- dim(sub), 
    `%`=round(100*fil/raw, 2))
rownames(tbl) <- c("gs", "cs"); tbl

# saving
saveRDS(sub, args[[2]])
