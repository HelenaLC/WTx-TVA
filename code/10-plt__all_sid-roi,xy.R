ps <- lapply(args[[1]], \(x) {
    # loading
    sce <- readRDS(x)
    # wrangling
    typ <- gsub("^.*_", "", sce$typ)
    typ <- factor(typ, names(.pal_roj))
    names(typ) <- colnames(sce)
    # plotting
    sid <- as.character(sce$sid[1])
    .plt_xy(sce, typ, sid, split=FALSE, na=TRUE) + 
        scale_color_manual("region", values=.pal_roj, 
            na.translate=FALSE, na.value="grey")
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
qpdf::pdf_combine(unlist(tf), output=args[[2]])
