rds <- list.files("meta/cnv", "ad$", full.names=TRUE)
names(rds) <- ids <- gsub(".*([0-9]{3}).*", "\\1", rds)

bin <- read.csv("meta/cnv/binning.csv")
bin <- bin[order(bin[[2]]), ]

gs <- read.csv("meta/cnv/mapping.csv.gz")
id <- gs[, 1]; gs <- gs[, -1]
bs <- paste0("bin_", gsub("^.*_", "", id))
gs <- setNames(asplit(unname(gs), 1), bs)

# 240 results include both, 
# lymph node section 241
# and colon section 242
for (. in setdiff(ids, "242")) {
    sce <- zellkonverter::readH5AD(rds[.])
    idx <- cut(seq(nrow(sce)), c(bin[[2]], nrow(sce)))
    rowData(sce)$chr <- bin[[1]][as.integer(idx)]
    rowData(sce)$ids <- gs[rownames(sce)]
    # retrieve corresponding cell metadata
    if (. != "240") {
        cd <- colData(readRDS(sprintf("outs/fil-%s.rds", .)))
        colData(sce) <- cd[colnames(sce), ]
        saveRDS(sce, sprintf("outs/cnv-%s.rds", .))
    } else {
        ist <- readRDS("outs/lv1-242.rds")$clust
        epi <- names(ist)[ist == "epi"]
        epi <- intersect(colnames(sce), epi)
        epj <- setdiff(colnames(sce), epi)
        idx <- list(`242`=epi, `241`=epj)
        for (. in names(idx)) {
            cd <- colData(readRDS(sprintf("outs/fil-%s.rds", .)))
            sub <- sce[, idx[[.]]]
            colData(sub) <- cd[colnames(sub), ]
            saveRDS(sub, sprintf("outs/cnv-%s.rds", .))
        }
    }
}
