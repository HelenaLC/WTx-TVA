# loading
cnv <- readRDS(args[[1]])
sce <- readRDS(grep("epi", args[[2]], value=TRUE))

# plotting
i <- match(colnames(sce), colnames(cnv))
y <- .z(unname(colMeans(assay(cnv))[i]))
p <- .plt_xy(sce, y, wcs$x, split=FALSE) + 
    guides(col=guide_colorbar("z-scaled\nCNV score"))

# saving
df <- p$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
ggsave(args[[3]], p, unit="cm", width=2+dx/2, height=0.5+dy/2)
