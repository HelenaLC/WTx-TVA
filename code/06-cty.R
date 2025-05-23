# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- lapply(args[[1]], readRDS)
ist <- lapply(args[[2]], readRDS)

# wrangling
.sid <- \(.) gsub(".*([0-9]{3}).*", "\\1", .)
jst <- split(ist, .sid(args[[2]]))
names(sce) <- .sid(args[[1]])
set.seed(250505)

# neighborhoods
df <- mapply(x=sce, y=jst, \(x, y) {
    kid <- lapply(y, \(.) .$clust)
    sub <- rep.int(
        c("epi", "imm", "str"), 
        sapply(kid, length))
    x <- x[, names(kid <- unlist(kid))]
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(x)))
    xy <- as.matrix(colData(x)[xy])
    nn <- nn2(xy, searchtype="radius", r=0.05, k=301)
    is <- nn$nn.idx[, -1]
    is[is == 0] <- NA
    # blind analysis to...
    kid <- as.character(kid)
    kid[grep("^epi", kid)] <- "epi"
    # ...epithelial subpopulations
    id <- matrix(kid[is], nrow=ncol(x))
    id[is.na(id)] <- ""
    names(ks) <- ks <- unique(kid)
    ns <- sapply(ks, \(k) rowSums(id == k))
    fq <- prop.table(ns, 1)
    data.frame(sub, kid,
        cid=colnames(x), 
        sid=factor(x$sid),
        fq, check.names=FALSE)
}, SIMPLIFY=FALSE) |> bind_rows()
df[is.na(df)] <- 0

# clustering
fd <- select(df, where(is.numeric))
km <- kmeans(fd, centers=nk <- 10)$cluster
df$ctx <- factor(km, labels=paste0("N", seq_len(nk)))
with(df, table(ctx, sid))

# saving
for (. in seq_along(sce)) {
    fd <- df[df$sid == names(sce)[.], ]
    base::saveRDS(fd, args[[2+.]])
}
