# loading
cnv <- readRDS(args[[1]])
sce <- readRDS(grep("epi", args[[2]], value=TRUE))

# wrangling
idx <- intersect(colnames(cnv), colnames(sce))
sce <- sce[, idx]; cnv <- cnv[, idx]
ids <- setNames(cnv$clu, idx)

# plotting
pal <- unname(pals::kelly(nlevels(cnv$clu)))
ps <- .plt_xy(sce, ids, wcs$x, split=TRUE, na=FALSE)
ps[[1]] <- ps[[1]] + scale_color_manual(values=pal) 
ps[[1]]$guides$guides$colour$params$override.aes$size <- 1

# saving
df <- ps[[1]]$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=4+dx/2, height=0.5+dy/2)
