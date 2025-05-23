ps <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        df <- readRDS(x)
        sce <- readRDS(y)
        # plotting
        sid <- paste(sce$sid[1])
        ctx <- setNames(df$ctx, df$cid)
        cty <- ctx; cty[!cty %in% c("N4", "N5")] <- NA
        p <- .plt_xy(sce, ctx, sid)
        q <- .plt_xy(sce, cty, sid, na=TRUE, split=FALSE)
        # aesthetics
        q$layers[[1]]$show.legend <- TRUE
        p[[1]]$layers[[1]]$show.legend <- TRUE
        p[[1]]$guides$guides$colour$params$override.aes$size <- 1
        pal <- scale_color_manual(
            "niche", drop=FALSE, 
            values=.pal_ctx, na.value="grey90")
        p[[1]] <- p[[1]] + pal; q <- q + pal
        c(p[1], list(q), p[-1])
    }) |> Reduce(f=c)

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