# loading
sce <- readRDS(args[[1]])

# plotting
ids <- setNames(sce$clu, colnames(sce))
pal <- unname(pals::kelly(nlevels(ids)))
ps <- .plt_xy(sce, ids, wcs$x, split=TRUE, na=FALSE)
ps[[1]] <- ps[[1]] + scale_color_manual(NULL, values=pal) 
ps[[1]]$guides$guides$colour$params$override.aes$size <- 1

# saving
df <- ps[[1]]$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[2]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
