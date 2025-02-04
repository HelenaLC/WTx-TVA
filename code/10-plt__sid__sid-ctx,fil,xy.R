# loading
df <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# plotting
fd <- df[df$sid == wcs$sid, ]
ctx <- setNames(fd$ctx, fd$cid)
ps <- .plt_xy(sce, ctx, wcs$sid, split=FALSE)

# saving
df <- ps[[1]]$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
