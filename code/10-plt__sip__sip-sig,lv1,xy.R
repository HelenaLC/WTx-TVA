# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
sub <- rowData(sce)$sub
gs <- rep.int(names(sub), sapply(sub, length))
gs <- split(gs, unlist(sub))
sapply(gs, length)

# analysis
ks <- (ks <- ist$clust)[match(colnames(sce), names(ks))]
cs <- split(colnames(sce), ks)
cs <- cs[-2]; gs <- gs[-2] # skip imm
mu <- mapply(i=gs, j=cs, SIMPLIFY=FALSE, \(i, j) {
    # get nearest neighbors
    xy <- grep("global_mm", names(colData(sce)))
    xy <- as.matrix(colData(sce)[xy][j, ])
    nn <- nn2(xy, searchtype="radius", r=0.05, k=100)
    # smooth scores locally
    is <- nn$nn.idx[, -1]
    is[is == 0] <- NA
    es <- assay(sce[i, j])
    apply(es, 1, \(y) {
        z <- matrix(y[c(is)], nrow(is), ncol(is))
        rowMeans(z, na.rm=TRUE)
    }) |> `rownames<-`(value=j)
})

# plotting
ps <- mapply(y=mu, j=cs, SIMPLIFY=FALSE, \(y, j) { 
    val <- rep(NA, ncol(sce))
    names(val) <- colnames(sce)
    lapply(colnames(y), \(g) {
        val[rownames(y)] <- y[, g]
        .plt_xy(sce, .z(unname(val))) + 
            theme(legend.position="none") + 
            ggtitle(g) + scale_color_gradientn(NULL, 
            colors=pals::jet(), na.value="lightgrey") 
    })
}) |> Reduce(f=c)

# saving
df <- ps[[1]]$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
