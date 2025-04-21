ps <- lapply(args[[1]], \(x) {
    sid <- gsub(".*([0-9]{3}).*", "\\1", x)
    .plt_rgb(readRDS(x), sid)
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
qpdf::pdf_combine(unlist(tf), output=args[[2]])
