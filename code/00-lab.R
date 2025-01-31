# dependencies
suppressPackageStartupMessages({
    library(jsonlite)
    library(SingleCellExperiment)
})

# loading
ist <- readRDS(args[[1]])
lab <- fromJSON(args[[2]])

# labelling
old <- unlist(lab)
int <- sapply(lab, length)
new <- rep.int(names(lab), int)
tbl <- data.frame(old, new)
jst <- .lab_ist(ist, tbl)
table(ist$clust, jst$clust)

# saving
saveRDS(jst, args[[3]])