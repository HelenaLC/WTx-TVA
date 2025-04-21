# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

ps <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        raw <- readRDS(x)
        fil <- readRDS(y)

        # wrangling
        ex <- !colnames(raw) %in% colnames(fil)
        ex <- setNames(ex, colnames(raw))
        df <- data.frame(colData(raw), ex)
        df <- df[order(df$ex), ]

        # plotting
        sid <- paste(raw$sid[1])
        n <- ncol(raw)-ncol(fil)
        n <- format(n, big.mark=",")
        p <- round(100*mean(!ex), 2)
        .plt_xy(df, df$ex, split=FALSE) +
            theme(legend.position="none") +
            scale_color_manual(values=c("lavender", "purple")) +
            ggtitle(bquote(bold(.(sid))~"(N ="~.(n)*";"~.(p)*"%)"))
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