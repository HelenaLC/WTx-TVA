ps <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        df <- readRDS(x)
        sce <- readRDS(y)
        # plotting
        sid <- paste(sce$sid[1])
        ctx <- setNames(df$ctx, df$cid)
        .plt_xy(sce, ctx, sid, split=FALSE)
    })

# aesthetics
ps <- lapply(ps, \(.) {
    .$layers[[1]]$show.legend <- TRUE
    .$guides$guides$colour$params$override.aes$size <- 1
    . + scale_color_manual("niche", drop=FALSE, values=.pal_ctx)
})

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- ps[[.]]$data
    dx <- diff(range(df$x))
    dy <- diff(range(df$y))
    pdf(tf[[.]], 
        width=(2+dx/2)/2.54, 
        height=(0.5+dy/2)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
