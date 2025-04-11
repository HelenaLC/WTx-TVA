# loading
df <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# plotting
ctx <- setNames(df$ctx, df$cid)
gg <- .plt_xy(sce, ctx, wcs$sid, split=FALSE) +
    scale_color_manual("niche", values=.pal_ctx)

# saving
df <- gg$data; dx <- diff(range(df$x)); dy <- diff(range(df$y))
gs <- gridExtra::marrangeGrob(grobs=list(gg), nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=2+dx/2, height=0.5+dy/2)
