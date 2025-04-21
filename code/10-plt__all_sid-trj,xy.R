# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(slingshot)
    library(SingleCellExperiment)
})

ps <- lapply(args[[1]], \(x) {
    # loading
    sce <- readRDS(x)
    sid <- gsub(".*([0-9]{3}).*", "\\1", x)
    # wrangling
    ex <- grep("sling", names(cd <- colData(sce)))
    pt <- slingAvgPseudotime(sce$slingshot)
    df <- data.frame(cd[-ex], pt)
    # plotting
    .plt_xy(df, unname(.q(pt)), sid) +
        scale_color_gradientn(
            "q-scaled\npseudotime",
            colors=pals::jet(), limits=c(0, 1),
            n.breaks=2, labels=c("early", "late"))
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
