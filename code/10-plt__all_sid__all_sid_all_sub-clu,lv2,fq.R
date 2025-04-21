# loading
library(SingleCellExperiment())
ps <- mapply(
    x=args[[1]], y=grep("epi", args[[2]], value=TRUE),
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)$clust
        ist <- gsub("^epi\\.", "", ist)
        # wrangling
        sid <- gsub(".*([0-9]{3}).*", "\\1", x)
        idx <- intersect(colnames(sce), names(ist))
        df <- data.frame(ist=ist[idx], clu=sce[, idx]$clu)
        # plotting
        n <- length(unique(df$clu))
        m <- length(unique(df$ist))
        p <- .plt_fq(df, "clu", "ist", sid, h=TRUE) + 
            labs(y=NULL, fill="epi")
        q <- .plt_fq(df, "ist", "clu", h=TRUE) + 
            scale_fill_manual("cnv", values=pals::kelly()) +
            theme(plot.title=element_blank())
        aes <- list(shape=21, stroke=0, size=1)
        wrap_plots(p, q, ncol=1, heights=c(n, m),guides="collect") & 
            guides(fill=guide_legend(override.aes=aes)) &
            theme(
                legend.margin=margin(),
                legend.spacing=unit(0.2, "lines"),
                legend.justification=c(0, 0.5),
                axis.text.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_text(hjust=1))
    })

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    h <- sum(sapply(ps[[.]], \(.) length(unique(.$data[[2]]))))
    pdf(tf[[.]], width=9/2.54, height=(0.5+h/5)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
